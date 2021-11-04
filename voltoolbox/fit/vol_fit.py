import math
import bisect
import datetime as dt
import numpy as np

from typing import Tuple, Dict, Optional
from dataclasses import dataclass
from itertools import groupby
from scipy.optimize import least_squares
from voltoolbox import (BusinessTimeMeasure, NormalizedSmileCurve, 
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
    put_vols = filter_vol_slice(put_vols, time_to_maturity, 2.0, (-5.0, 1.0))

    call_vols = _slice_implied_vol(quote_slice.call, 1.0, forward, time_to_maturity)
    call_vols = filter_vol_slice(call_vols, time_to_maturity, -1.0, (-1.0, 5.0))

    return VolQuoteSlice(quote_slice.symbol, quote_slice.expiry, time_to_maturity, forward, call_vols, put_vols)


def prepare_target_vol(vol_slice: VolQuoteSlice )-> VolSlice:
    put = vol_slice.put
    vol_datas = list(zip(put.log_moneyness, put.mids, put.errs, ['p'] * len(put.mids)))
    call = vol_slice.call
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
    def create(cls, vol_sl: VolQuoteSlice):
        target_vols = prepare_target_vol(vol_sl)
        mids = np.array(target_vols.mids)
        zs = np.array(target_vols.log_moneyness) / (math.sqrt(vol_sl.time_to_maturity) * mids)
        return cls(vol_sl.time_to_maturity, zs, mids, np.array(target_vols.errs))


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
        max_diff = np.diff(bulk_zs).max()
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

    def smile_curve(self):
        vols = self.atm_vol * self.vol_ratios
        return NormalizedSmileCurve(self.zs, vols)


class SurfaceBackbone:
    def __init__(self):
        self.deviations = None
        self.expiries = []
        self.slices = []

    def add_slice(self, smile: SmileBackbone) -> Tuple[Optional[SmileBackbone], Optional[SmileBackbone]]:
        if self.deviations is None: 
            self.deviations = smile.zs
        else:
            if not np.array_equal(self.deviations, smile.zs):
                raise Exception('invalid slice')
        self.slices.append(smile)
        self.slices = sorted(self.slices, key=lambda sl: sl.time_to_expiry)
        self.expiries = [sl.time_to_expiry for sl in self.slices]

    def bounding_slices(self, t):
        idx = bisect.bisect_left(self.expiries, t)

        if idx==0:
            return (None, self.slices[idx])
        
        if idx==len(self.expiries):
            return (self.slices[idx-1], None)

        return (self.slices[idx-1], self.slices[idx])

    def interpolated_slice(self, t, *, atm_vol=float('nan')) -> SmileBackbone:
        left_sl, right_sl = self.bounding_slices(t)

        if right_sl is None:
            raise Exception('no extrapolation allowed')

        if left_sl is None:
            if np.isnan(atm_vol):
                atm_vol = right_sl.atm_vol
            return SmileBackbone(t, atm_vol, self.deviations, right_sl.vol_ratios) 

        w = (t - left_sl.time_to_expiry) / (right_sl.time_to_expiry - left_sl.time_to_expiry)

        if np.isnan(atm_vol):
            atm_vol = np.sqrt(w * right_sl.atm_vol**2 + (1.0 - w) * left_sl.atm_vol**2)
        
        vrs = np.sqrt(w * right_sl.vol_ratios**2 + (1.0 - w) * left_sl.vol_ratios**2)
        return SmileBackbone(t, atm_vol, self.deviations, vrs)
