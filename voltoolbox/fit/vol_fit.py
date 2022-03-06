import math
import datetime as dt
import numpy as np

from typing import Tuple, Dict
from dataclasses import dataclass
from itertools import groupby
from scipy.optimize import least_squares
from voltoolbox import SmileVariationFilter, longest_increasing_subsequence
from .option_quotes import VolQuoteSlice, VolSlice
from .average_spline import AverageSpline
from .vol_surface import SurfaceBackbone, SmileBackbone


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
                assert prev_q[3]=='p'
                if log_moneyness < 0.0:
                    q0 = prev_q
                    q1 = q
                else:
                    q0 = q
                    q1 = prev_q             
            else :
                assert prev_q[3]=='c'
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

    for _ in range(0, 3):
        atm_mask = (target_zs > -fit_width) & (target_zs < fit_width)
        target_zs = target_zs[atm_mask]

        if len(target_zs) > 4: # Enough strike
            target_mids = target_mids[atm_mask]
            target_errs = target_errs[atm_mask]
            weights = (np.maximum(0.0, 1.0 - (target_zs / fit_width))**2.0) * 3.0 / 4.0  # Epanechnikov kernel
            break
        else:  # rescale with larger fit_width
            fit_width *= 1.6
            target_zs = target_sl.zs

    avg_spline = AverageSpline.average_spline()
    coeffs = avg_spline.fit_coeffs(target_zs, target_mids, target_errs / np.sqrt(weights), smoothing = 2.0e-7)

    atm_vol = avg_spline.sample_basis(np.array([0.0])).dot(coeffs)[0]

    fitted_targets = avg_spline.sample_basis(target_zs).dot(coeffs)
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

        # Curve smoothing penality term
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
        # Score is build on three criteria : wing width, nb strikes, vol spread
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
    for _, ts in groupby(expis, key=lambda item : -item[1]):
        res += reversed(sorted(((t, s) for t, s in ts)))

    return res


def fit_target_slices(target_slices: Dict[dt.datetime, TargetSlice],
                      *, fit_zs) -> SurfaceBackbone:
    # TODO clean-up all magic constants
    atm_vols = fit_atm_vol_curve(target_slices, smoothing=1.0e-3)

    fitted_surf = SurfaceBackbone(fit_zs)
    expis = sorted_fit_expiries(target_slices)
    desarb_targets = {}
    for expi_dt, _ in expis:
        target_sl = target_slices[expi_dt]
        atm_vol = atm_vols[expi_dt]

        calendar_sandwich = None 
        if len(fitted_surf.expiries) > 0:
            calendar_sandwich = fitted_surf.calendar_arb_slice(target_sl.t)

        if calendar_sandwich:    
            desarb_mids, _, _, _ = calendar_sandwich.arbitrage_free_vol(target_sl.zs, target_sl.mids)            
        else:
            desarb_mids = target_sl.mids
        desarb_targets[expi_dt] = desarb_mids

        # BUILD A PRIOR
        if len(fitted_surf.expiries)==0:
            avg_spline = AverageSpline.average_spline([1.0, 1.75, 3.0],
                                                      [1.0, 2.0])

            coeffs = avg_spline.fit_coeffs(target_sl.zs, desarb_mids, target_sl.errs, smoothing=1e-7)
            basis = avg_spline.sample_basis(fit_zs)
            prior_smile = SmileBackbone(target_sl.t, atm_vol, fit_zs,  basis.dot(coeffs) / atm_vol)
        else:
            prior_smile = fitted_surf.slice(target_sl.t, atm_vol = atm_vol)
        prior_crv = prior_smile.normalized_curve()
        prior_vols = np.vectorize(prior_crv.vol)(target_sl.zs)

        # FILTER MID - PRIOR
        prior_vol_diffs = desarb_mids - prior_vols
        noise_devs = np.maximum(5.0 / 10000.0, 0.5 * target_sl.errs)
        ATM_DEV = 0.01
        ATM_SKEW_RATIO = 2.0
        WING_SKEW_RATIO = 0.5
        Z_REF = 2.5
        dvol_filter = SmileVariationFilter(target_sl.zs,
                                           prior_vol_diffs,
                                           noise_devs,
                                           ATM_DEV,
                                           ATM_SKEW_RATIO,
                                           WING_SKEW_RATIO,
                                           Z_REF)
        filtered_vol_diffs = np.vectorize(dvol_filter.dvol)(prior_smile.zs)
        filtered_atm_vol_diff = dvol_filter.dvol(0.0)

        fitted_vols = np.vectorize(prior_crv.vol)(prior_smile.zs) + filtered_vol_diffs
        fitted_atm_vol = prior_crv.vol(0.0) + filtered_atm_vol_diff

        if calendar_sandwich:
            fitted_vols, _, _, _ = calendar_sandwich.arbitrage_free_vol(prior_smile.zs, fitted_vols)
            atf_vs, _, _, _ = calendar_sandwich.arbitrage_free_vol(np.array([0.0]), np.array([fitted_atm_vol]))
            fitted_atm_vol = atf_vs[0]

        fitted_surf.add_pillar_slice(target_sl.t, fitted_atm_vol, fitted_vols / fitted_atm_vol)

    return fitted_surf
