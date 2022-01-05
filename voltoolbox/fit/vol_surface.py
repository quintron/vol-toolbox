import bisect
import numpy as np

from typing import Tuple, Optional
from dataclasses import dataclass
from voltoolbox import SplineCurve, SmileSplineCurve


@dataclass(frozen=True)
class SmileBackbone:
    time_to_expiry: float
    atm_vol: float
    zs: np.array
    vol_ratios: np.array

    def normalized_curve(self):
        vols = self.atm_vol * self.vol_ratios
        return SplineCurve(self.zs, vols)

    def sample_smile_curve(self, sampling_zs: np.array=None):
        if sampling_zs is None:
            mid_pts = 0.5 * (self.zs[1:] + self.zs[:-1])
            sampling_zs = np.append(self.zs, mid_pts)

        vols = np.vectorize(self.normalized_curve().vol)(sampling_zs)
        xs = sampling_zs * vols * np.sqrt(self.time_to_expiry)
        return SmileSplineCurve(xs, vols)


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
        
        lower_crv = lower_sl.normalized_curve()
        self.lower_vol_fn = np.vectorize(lower_crv.vol)
        self.lower_t = lower_sl.time_to_expiry

        upper_crv = upper_sl.normalized_curve()
        self.upper_vol_fn = np.vectorize(upper_crv.vol)
        self.upper_t = upper_sl.time_to_expiry

    def arbitrage_free_vol(self, zs, vols):
        '''Return (desarb_vol, score, min_vol, max_vol) for calendar arbitrage.
           Score is an arbitrage indicator, if its value is in [-1, 1] there is no calendar arb.'''       
        assert len(zs) == len(vols)
        min_var = self.lower_vol_fn(zs)**2 * self.lower_t
        min_vol = np.sqrt(min_var / self.t)

        max_var = self.upper_vol_fn(zs)**2 * self.upper_t
        max_vol = np.sqrt(max_var / self.t)

        need_fix = min_vol >= max_vol
        if need_fix.any():
            avg_vol = 0.5 * (min_vol + max_vol)
            delta_vol = 0.5 * np.abs(max_vol - min_vol)
            EPS = 0.10
            fixed_min_vol = np.select([need_fix, ~need_fix],
                                    [avg_vol - EPS * delta_vol, min_vol])
            fixed_max_vol = np.select([need_fix, ~need_fix],
                                    [avg_vol + EPS * delta_vol, max_vol])
            min_vol = fixed_min_vol
            max_vol = fixed_max_vol
            min_var = min_vol**2 * self.t
            max_var = max_vol**2 * self.t

        w = (self.t - self.lower_t) / (self.upper_t - self.lower_t) 
        mid_var = w * max_var + (1.0 - w) * min_var

        var = vols**2 * self.t
        score = np.select([var > mid_var], 
                          [(var - mid_var) / (max_var - mid_var)],
                          (var - mid_var) / (mid_var - min_var))

        # map score into [-1, 1]
        desarb_score = np.select([score > 0.5, score <= -0.5],
                                 [1.0 - np.exp(-4.0 * (score - 0.5)) / (1.0 + np.exp(-4.0 * (score - 0.5))), 
                                  -(1.0 - np.exp(-4.0 * (-score - 0.5)) / (1.0 + np.exp(-4.0 * (-score - 0.5))))],
                                 score)

        desarb_var = np.select([desarb_score > 0.0],
                               [mid_var + desarb_score * (max_var - mid_var)],
                                mid_var + desarb_score * (mid_var - min_var))

        desarb_vol = np.sqrt(desarb_var / self.t)
        return desarb_vol, score, min_vol, max_vol  
