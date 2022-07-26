import numpy as np
from scipy.interpolate import interp1d
from voltoolbox.utils import act365_time


def _rate_interpolation(ts, rates, kind='cubic'):
    total_yield = interp1d(np.append([0.0], ts),
                           np.append([0.0], ts * rates), kind=kind)
    def rate_fn(t):
        if t <= 0.0:
            return 0.0
        else:
            return total_yield(t) / t
    return rate_fn


class DiscountRateCurve:

    def __init__(self,
                 pricing_date,
                 maturities,
                 rates,
                 *, interpol_type: str='cubic'):
        '''
        interpol_type: 'cubic' or 'linear
          '''

        self.pricing_date = pricing_date
        ts = np.array([act365_time(pricing_date, m) for m in maturities])
        self.rate_fn = _rate_interpolation(ts, rates, kind = interpol_type)

    def discount_rate(self, d):
        return self.rate_fn(act365_time(self.pricing_date, d))

    def discount_factor(self, d):
        t = act365_time(self.pricing_date, d)
        return np.exp(- t * self.rate_fn(t))
