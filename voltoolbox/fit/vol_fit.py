import math
import bisect
import datetime as dt
import numpy as np

from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
from itertools import groupby
from scipy.optimize import least_squares
from voltoolbox import (BusinessTimeMeasure, NormalizedSmileCurve, SplineSmileCurve, SmileVariationFilter,
                        bs_implied_volatility, longest_increasing_subsequence)
from voltoolbox.fit.option_quotes import OptionQuoteSlice, QuoteSlice, VolQuoteSlice, VolSlice
from voltoolbox.fit.fit_utils import act365_time, filter_quotes
from voltoolbox.fit.average_spline import AverageSpline


def _slice_implied_vol(quotes: QuoteSlice,
                       option_type: float,
                       forward: float,
                       time_to_maturity: float) -> VolSlice:
    xs = []
    vol_mids = []
    vol_errs = []
    for k, b, a in zip(quotes.strikes, quotes.bids, quotes.asks):
        vol_bid = bs_implied_volatility(forward, k, b, time_to_maturity, option_type)
        vol_ask = bs_implied_volatility(forward, k, a, time_to_maturity, option_type)
        if vol_bid > 0.0 and vol_ask > vol_bid:
            vol_mids.append(0.5 * (vol_ask + vol_bid))
            vol_errs.append(0.5 * (vol_ask - vol_bid))
            xs.append(math.log(k / forward))
    return VolSlice(tuple(xs),
                    tuple(vol_mids),
                    tuple(vol_errs))

def filter_vol_slice(vol_slice: VolSlice, time_to_maturity: float, 
                     power_payoff_existence: float,
                     std_dev_range: Tuple[float, float])-> VolSlice:
    sqrt_t = math.sqrt(time_to_maturity)

    xs = np.array(vol_slice.log_moneyness)
    mids = np.array(vol_slice.mids)
    errs = np.array(vol_slice.errs)

    inc_d_plus = longest_increasing_subsequence(xs / (sqrt_t * mids) - (power_payoff_existence - 0.5) * sqrt_t * mids)
    mids = mids[inc_d_plus]
    errs = errs[inc_d_plus]
    xs = xs[inc_d_plus]

    z_min, z_max = std_dev_range
    range_filter = (xs / (sqrt_t * mids) >= z_min) & (xs / (sqrt_t * mids) <= z_max)
    mids = mids[range_filter]
    errs = errs[range_filter]
    xs = xs[range_filter]

    return VolSlice(tuple(xs),
                    tuple(mids),
                    tuple(errs))


def prepare_vol_quotes(quote_slice: OptionQuoteSlice,
                       forward: float,
                       box_spread: float, 
                       pricing_dt: dt.datetime,
                       time_measure: BusinessTimeMeasure) -> VolQuoteSlice:
    assert(quote_slice.expiry > pricing_dt)

    box_discount = quote_slice.discount * math.exp(-box_spread * act365_time(pricing_dt, quote_slice.expiry))
    quote_slice = OptionQuoteSlice(quote_slice.symbol,
                                   quote_slice.expiry,
                                   box_discount,
                                   quote_slice.call,
                                   quote_slice.put)
    quote_slice = filter_quotes(quote_slice,
                                (forward * (1 - 0.0e-6), forward * (1 + 0.0e-6)),
                                discounted_output_quotes=True)

    time_to_maturity = time_measure.distance(pricing_dt, quote_slice.expiry)
    put_vols = _slice_implied_vol(quote_slice.put, -1.0, forward, time_to_maturity)
    call_vols = _slice_implied_vol(quote_slice.call, 1.0, forward, time_to_maturity)
    return VolQuoteSlice(quote_slice.symbol, quote_slice.expiry, time_to_maturity, forward, call_vols, put_vols)


def prepare_target_vol(vol_slice: VolQuoteSlice )-> VolSlice:

    put = filter_vol_slice(vol_slice.put, vol_slice.time_to_maturity, 2.0, (-5.0, 1.0))
    call = filter_vol_slice(vol_slice.call, vol_slice.time_to_maturity, -1.0, (-1.0, 5.0))

    vol_datas = list(zip(put.log_moneyness, put.mids, put.errs, ['p'] * len(put.mids)))
    vol_datas += list(zip(call.log_moneyness, call.mids, call.errs, ['c'] * len(call.mids)))
    vol_datas = sorted(vol_datas, key = lambda q: q[0])

    merged_datas = [vol_datas[0]]
    prev_q = vol_datas[0]
    for q in vol_datas[1:]:
        log_moneyness, mid, err, opt_type = q
        if abs(log_moneyness - prev_q[0]) > 1.0e-5:
            merged_datas.append((log_moneyness, mid, err))
        else:            
            if opt_type=='c':
                assert(prev_q[3]=='p')
                if log_moneyness < 0.0:
                    q0 = prev_q
                    q1 = q
                else:
                    q0 = q
                    q1 = prev_q                
            else :
                assert(prev_q[3]=='c')
                if log_moneyness < 0.0:
                    q0 = q
                    q1 = prev_q
                else:
                    q0 = prev_q
                    q1 = q
            
            diff = abs(q0[1] - q1[1])
            w = (1.0 - min(1.0, diff / q0[2])) * (1.0 - min(1.0, diff / q1[2]))
            w0 = 1.0 - 0.5 * w
            w1 = 1.0 - w0
            merged_datas.append((q0[0], w0 * q0[1] + w1 * q1[1], w0 * q0[2] + w1 * min(q0[2], q1[2])))
        prev_q = q

    xs = [q[0] for q in merged_datas]
    mids = [q[1] for q in merged_datas]
    errs = [q[2] for q in merged_datas]
    target_vols = VolSlice(xs, mids, errs)
    target_vols = filter_vol_slice(target_vols, vol_slice.time_to_maturity, 0.5, (-10.0, 10.0))
    return target_vols


@dataclass(frozen=True)
class TargetSlice():
    t: float
    zs: np.array
    mids: np.array
    errs: np.array

    @classmethod
    def create(cls, volquote_sl: VolQuoteSlice):

        vol_sl = prepare_target_vol(volquote_sl)

        sqrt_t = math.sqrt(volquote_sl.time_to_maturity)
        xs = np.array(vol_sl.log_moneyness)
        mid_vols = np.array(vol_sl.mids)
        err_vols = np.array(vol_sl.errs)

        bid_vols = mid_vols - err_vols
        bid_zs = xs / (sqrt_t * bid_vols)
        
        ask_vols = mid_vols + err_vols
        ask_zs = xs / (sqrt_t * ask_vols)

        mid_zs = 0.5 * (bid_zs + ask_zs)

        return cls(volquote_sl.time_to_maturity, 
                   mid_zs, mid_vols, err_vols)


def fit_atm_vol(target_sl: TargetSlice, 
                *, fit_width: float=0.6) -> Tuple[float, float]:
    target_zs = target_sl.zs
    target_mids = target_sl.mids
    target_errs = target_sl.errs

    for i in range(0, 3):
        atm_mask = (target_zs > -fit_width) & (target_zs < fit_width)
        target_zs = target_zs[atm_mask]
        
        if len(target_zs) > 4: # Enough strike
            target_mids = target_mids[atm_mask]
            target_errs = target_errs[atm_mask]
            weights = (np.maximum(0.0, 1.0 - (target_zs / fit_width))**2.0) * 3.0 / 4.0  # Epanechnikov kernel
            break
        else: # rescale with larger fit_width 
            fit_width *= 1.6
            target_zs = target_sl.zs

    avg_spline = AverageSpline()
    coeffs = avg_spline.fit_coeffs(target_zs, target_mids, target_errs / np.sqrt(weights), smoothing = 2.0e-7)

    atm_vol = avg_spline.compute_vals(np.array([0.0])).dot(coeffs)[0]

    fitted_targets = avg_spline.compute_avg_vals(target_zs).dot(coeffs)
    atm_vol_err = ((target_errs + np.abs(target_mids - fitted_targets)) * weights).sum() / weights.sum()

    return (atm_vol, atm_vol_err)

    
class VolCurveFitFunction:

    def __init__(self, ts, vols, errs, smoothing):
        assert len(ts)==len(vols)==len(errs)
        self.ts = np.array(ts)
        self.vols = vols
        self.errs = errs
        self.smoothing = smoothing

        ref_fwd_vars = [] 
        prev_t = 0.0
        prev_var = 0.0
        for t, v in zip(ts, vols):
            fwd_var = np.maximum(0.10 * v * (t - prev_t), (t * v**2 - prev_var))
            prev_t = t
            prev_var += fwd_var
            ref_fwd_vars.append(fwd_var)
        self.ref_fwd_vars = np.array(ref_fwd_vars)

    def fwd_vars(self, xs):
        return self.ref_fwd_vars * np.exp(xs)

    def vol_curve(self, xs):
        fwd_vars = self.fwd_vars(xs)
        vars = []
        current_var = 0.0
        for fwd_var in fwd_vars:
            current_var += fwd_var
            vars.append(current_var)
        return np.sqrt(np.array(vars) / self.ts)

    def fwd_vol_curve(self, xs):
        fwd_vars = self.fwd_vars(xs)
        dts =  np.ediff1d(np.append(np.array([0.0]), self.ts))
        return np.sqrt(fwd_vars / dts)

    def fit_residuals(self, xs):
        fit_vols = self.vol_curve(xs)
        scores = []
        for v, target_v, err in zip(fit_vols, self.vols, self.errs):
            scores.append((v - target_v) / err)

        #Curve smoothing penality term
        fwd_vols = self.fwd_vol_curve(xs)
        fwd_vol_smoothing = np.sqrt(self.smoothing) * np.ediff1d(fwd_vols) / np.ediff1d(self.ts)

        return np.append(np.array(scores), fwd_vol_smoothing)

    def __call__(self, xs):
        return self.fit_residuals(xs)


def fit_atm_vol_curve(target_slices :Dict[dt.datetime, TargetSlice],
                      *, smoothing :float=5.0e-3) -> Dict[dt.datetime, float]:
    ts = []
    atm_vol_mids = []
    atm_vol_errs = []
    for target_sl in target_slices.values():
        v, err = fit_atm_vol(target_sl)
        ts.append(target_sl.t)
        atm_vol_mids.append(v)
        atm_vol_errs.append(err)

    fit_func = VolCurveFitFunction(np.array(ts),
                                   np.array(atm_vol_mids),
                                   np.array(atm_vol_errs),
                                   smoothing)    
    x0 = np.array([0.0] * len(ts))
    res = least_squares(fit_func, x0)
    atm_vols = fit_func.vol_curve(res.x)
    return dict(zip(target_slices, atm_vols))


def sorted_fit_expiries(target_slices):
    slice_score_infos = {}
    for expi_dt, target_sl in target_slices.items():
        err_med = np.median(target_sl.errs)
        nb_strikes = len(target_sl.zs)

        bulk_zs = target_sl.zs[ np.abs(target_sl.zs) < 2.5]
        width = 0.5 *(abs(min(bulk_zs)) + max(bulk_zs))
        max_diff = np.quantile(np.diff(bulk_zs), 0.90)
        unif_diff = (max(bulk_zs) - min(bulk_zs)) / (len(bulk_zs) - 1)
        strike_uniformity =  max_diff / unif_diff

        slice_score_infos[expi_dt] = (width, err_med, nb_strikes, strike_uniformity)

    wing_score_threshold = np.median(np.array([w for w, _, _, _ in slice_score_infos.values()]))
    err_score_threshold = 2.0 * np.quantile(np.array([e for _, e, _, _ in slice_score_infos.values()]), 0.75)
    nb_strike_threshold = 0.75 * np.median(np.array([n for _, _, n, _ in slice_score_infos.values()]))
    strike_uniform_threshold = np.quantile(np.array([u for _, _, _, u in slice_score_infos.values()]), 0.5)

    slices_scores = {}
    for expi_dt, infos in slice_score_infos.items():
        width, err, nb_strikes, strike_uniformity = infos    
        # Score  1s build on three criteria : wing width, nb strikes, wvol spread(error) median
        # Weight on each criteria is chosen so that :
        #   1. wing width is more important that nb strike and vol spread
        #   2. nb strike is more important than vol spread
        score = 0.0
        if (width > wing_score_threshold):
            score += 2.0
        
        if (nb_strikes > nb_strike_threshold):
            score += 0.5

        if strike_uniformity < strike_uniform_threshold:
            score += 0.5

        if (err < err_score_threshold):
            score += 0.75

        slices_scores[expi_dt] = score

    #sorted by decreasing score
    res = []
    expis = sorted(slices_scores.items(), key=lambda item : -item[1])
    for sc, ts in groupby(expis, key=lambda item : -item[1]):
        res += reversed(sorted(((t, s) for t, s in ts)))

    return res


@dataclass(frozen=True)
class SmileBackbone:
    time_to_expiry: float
    atm_vol: float
    zs: np.array
    vol_ratios: np.array

    def normalized_curve(self):
        vols = self.atm_vol * self.vol_ratios
        return NormalizedSmileCurve(self.zs, vols)

    def sample_smile_curve(self, sampling_zs: np.array):
        vols = np.vectorize(self.normalized_curve().vol)(sampling_zs)
        xs = sampling_zs * vols * np.sqrt(self.time_to_expiry)
        return SplineSmileCurve(xs, vols)


class SurfaceBackbone:
    def __init__(self, deviations):
        self._deviations = deviations
        self.expiries = []
        self.slices = []

    @property
    def deviations(self):
        return self._deviations

    def add_pillar_slice(self, t: float, atm_vol: float, vol_ratios: np.array):
        if len(vol_ratios) != len(self._deviations):
                raise Exception('invalid slice format')

        self.slices.append(SmileBackbone(t, atm_vol, self._deviations, vol_ratios))
        self.slices = sorted(self.slices, key=lambda sl: sl.time_to_expiry)
        self.expiries = [sl.time_to_expiry for sl in self.slices]

    def bounding_pillar_slices(self, t) -> Tuple[Optional[SmileBackbone], Optional[SmileBackbone]]:
        if len(self.expiries)==0:
            return (None, None)
        
        idx = bisect.bisect_left(self.expiries, t)

        if idx==0:
            return (None, self.slices[idx])
        
        if idx==len(self.expiries):
            return (self.slices[idx-1], None)

        return (self.slices[idx-1], self.slices[idx])

    def slice(self, t, *, atm_vol=float('nan')) -> SmileBackbone:
        left_sl, right_sl = self.bounding_pillar_slices(t)

        if right_sl is None:
            ref_sl = self.slice(left_sl.time_to_expiry * 0.5)
            return _extrapolated_slice(t, ref_sl, left_sl, atm_vol=atm_vol)

        if left_sl is None:
            if np.isnan(atm_vol):
                atm_vol = right_sl.atm_vol
            return SmileBackbone(t, atm_vol, self._deviations, right_sl.vol_ratios) 

        return _interpolated_slice(t, left_sl, right_sl, atm_vol=atm_vol)

    def calendar_arb_slice(self, t):
        lower_sl, upper_sl = self.bounding_pillar_slices(t)

        if lower_sl is None:
            lower_sl = self.slice(0.10 * t)

        if upper_sl is None:
            upper_sl = self.slice(1.5 * t)

        return  CalendarSandwich(t, lower_sl, upper_sl)


def _interpolated_slice(t: float,
                       left_sl: SmileBackbone,
                       right_sl: SmileBackbone,
                       *, atm_vol=float('nan')):
    w = (t - left_sl.time_to_expiry) / (right_sl.time_to_expiry - left_sl.time_to_expiry)

    if np.isnan(atm_vol):
        atm_vol = np.sqrt(w * right_sl.atm_vol**2 + (1.0 - w) * left_sl.atm_vol**2)

    vrs = np.sqrt(w * right_sl.vol_ratios**2 + (1.0 - w) * left_sl.vol_ratios**2)
    return SmileBackbone(t, atm_vol,  left_sl.zs, vrs)


def _extrapolated_slice(t: float,
                        ref_sl: SmileBackbone,
                        left_sl: SmileBackbone,
                        *, atm_vol=float('nan')):
    assert(left_sl.time_to_expiry > ref_sl.time_to_expiry)
    fwd_var_ratios = (left_sl.vol_ratios**2 * left_sl.time_to_expiry - ref_sl.vol_ratios**2 * ref_sl.time_to_expiry) / ( left_sl.time_to_expiry - ref_sl.time_to_expiry)
    w = left_sl.time_to_expiry / t
    vrs = np.sqrt(left_sl.vol_ratios**2 * w + (1.0 - w) * fwd_var_ratios)
    
    if np.isnan(atm_vol):
        fwd_var_atm = (left_sl.atm_vol**2 * left_sl.time_to_expiry - ref_sl.atm_vol**2 * ref_sl.time_to_expiry) / ( left_sl.time_to_expiry - ref_sl.time_to_expiry)
        atm_vol = np.sqrt(left_sl.atm_vol**2 * w + (1.0 - w) * fwd_var_atm)

    return SmileBackbone(t, atm_vol,  left_sl.zs, vrs)


class CalendarSandwich:

    def __init__(self, 
                t: float, 
                lower_sl: SmileBackbone, 
                upper_sl: SmileBackbone) -> None:
        self.t = t
        
        self.lower_crv = lower_sl.normalized_curve()
        self.lower_t = lower_sl.time_to_expiry
        
        self.upper_crv = upper_sl.normalized_curve()
        self.upper_t = upper_sl.time_to_expiry

    def arbitrage_free_vol(self, z: float, vol: float):
        '''Return (desarb_vol, score, min_vol, max_vol) for calendar arbitrage.
           Score is an arbitrage indicator, if its value is in [-1, 1] there is no calendar arb.'''       
        min_var = self.lower_crv.vol(z)**2 * self.lower_t
        min_vol = np.sqrt(min_var / self.t)
        
        max_var = self.upper_crv.vol(z)**2 * self.upper_t 
        max_vol = np.sqrt(max_var / self.t)

        w = (self.t - self.lower_t) / (self.upper_t - self.lower_t) 
        mid_var = w * max_var + (1.0 - w) * min_var

        var = vol**2 * self.t
        if var > mid_var:
            score = (var - mid_var) / (max_var - mid_var)
        else:
            score = (var - mid_var) / (mid_var - min_var)

        # map score into [-1, 1]
        if (score > 0.5):
            desarb_score = 1.0 - np.exp(-4.0 * (score - 0.5)) / (1.0 + np.exp(-4.0 * (score - 0.5)))
        elif (score < -0.5):
            desarb_score = -(1.0 - np.exp(-4.0 * (-score - 0.5)) / (1.0 + np.exp(-4.0 * (-score - 0.5))))
        else:
            desarb_score = score

        if desarb_score > 0.0:
            desarb_var = mid_var + desarb_score * (max_var - mid_var)
        else:
            desarb_var = mid_var + desarb_score * (mid_var - min_var)

        desarb_vol = np.sqrt(desarb_var / self.t)
                
        return desarb_vol, score, min_vol, max_vol  


def fit_target_slices(target_slices: Dict[dt.datetime, TargetSlice],
                      *, fit_zs) -> SurfaceBackbone:
    # TODO clean-up all magic constants
    atm_vols = fit_atm_vol_curve(target_slices, smoothing=1.0e-3)

    fitted_surf = SurfaceBackbone(fit_zs)
    expis = sorted_fit_expiries(target_slices)
    desarb_targets = {}
    for expi_dt, score in expis:
        target_sl = target_slices[expi_dt]
        atm_vol = atm_vols[expi_dt]
        
        if len(fitted_surf.expiries) > 0:
            calendar_sandwich = fitted_surf.calendar_arb_slice(target_sl.t)
            desarb_mids = []
            for z, v in zip(target_sl.zs, target_sl.mids):
                desarb_v, score, min_vol, max_vol =calendar_sandwich.arbitrage_free_vol(z, v)
                desarb_mids.append(desarb_v)
            desarb_targets[expi_dt] = desarb_mids

            target_sl = TargetSlice(target_sl.t,
                                    target_sl.zs,
                                    np.array(desarb_mids),
                                    target_sl.errs)


        # BUILD A PRIOR
        if len(fitted_surf.expiries)==0:
            avg_spline = AverageSpline([3.25, 2.75, 2.25, 1.75, 1.0, 0.5, 0.25],
                                    [0.25, 0.5, 1.0, 1.5, 2.0, 2.5])
            coeffs = avg_spline.fit_coeffs(target_sl.zs, target_sl.mids, target_sl.errs, smoothing=1e-7)
            basis = avg_spline.compute_avg_vals(fit_zs)
            prior_smile = SmileBackbone(target_sl.t, atm_vol, fit_zs,  basis.dot(coeffs) / atm_vol)
        else:
            prior_smile = fitted_surf.slice(target_sl.t)
            prior_smile = SmileBackbone(target_sl.t, atm_vol, fit_zs, prior_smile.vol_ratios)

        prior_crv = prior_smile.normalized_curve()
        prior_vols = np.array([prior_crv.vol(z) for z in target_sl.zs])
        prior_vol_diffs = target_sl.mids - prior_vols

        # FILTER MID - PRIOR
        noise_devs = np.maximum(5.0 / 10000.0, 0.5 * target_sl.errs)
        atm_dev = 0.01
        atm_skew_dev = 0.02
        z_ref = 3.0
        dvol_filter = SmileVariationFilter(list(target_sl.zs), list(prior_vol_diffs), 
                                        list(noise_devs),
                                        atm_dev,
                                        atm_skew_dev,
                                        z_ref)
        filtered_vol_diffs = np.array([dvol_filter.dvol_with_error(z)[0] for z in prior_smile.zs])
        fitted_vols = prior_smile.atm_vol * prior_smile.vol_ratios + filtered_vol_diffs 
        fitted_atm_vol = prior_smile.atm_vol + dvol_filter.dvol_with_error(0.0)[0]

        fitted_surf.add_pillar_slice(target_sl.t, fitted_atm_vol, fitted_vols / fitted_atm_vol)

    return fitted_surf
