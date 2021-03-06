import math
import numpy as np
import datetime as dt
from typing import Dict, Optional, Tuple

from voltoolbox import bs_implied_volatility
from voltoolbox.fit.option_quotes import OptionSnapshot, OptionQuoteSlice, QuoteSlice
from fit_utils import act365_time, filter_quotes


def prepare_quotes_for_fit(quote_slice: OptionQuoteSlice,
                           pricing_date: dt.datetime,
                           spot_prior: float,
                           *, yield_threshold=0.10
                           ) -> OptionQuoteSlice:
    """Return discounted option quotes, outlier have been removed.
    """

    time_to_mat = act365_time(pricing_date, quote_slice.expiry)
    ratio_min = min(math.exp(- time_to_mat * yield_threshold), (1.0 - yield_threshold))
    ratio_max = math.exp(time_to_mat * yield_threshold)
    
    return filter_quotes(quote_slice, (spot_prior * ratio_min, spot_prior * ratio_max), discounted_output_quotes = True)


class OptionKernelRegression:

    def __init__(self, opt_slice: QuoteSlice):
        self.opt_slice = opt_slice

        call_ks = np.array(opt_slice.call.strikes)
        call_bids = np.array(opt_slice.call.bids)
        call_asks = np.array(opt_slice.call.asks)
        put_bids = np.array(opt_slice.put.bids)
        put_asks = np.array(opt_slice.put.asks)
        
        self.ks = np.array(opt_slice.call.strikes + opt_slice.put.strikes)
        self.put_target = np.concatenate(
                            (0.5 * (call_bids + call_asks) + call_ks, 
                             0.5 * (put_bids + put_asks)
                            ), axis=0)        
        self.put_err = np.concatenate((0.5 * (call_asks - call_bids), 
                                       0.5 * (put_asks - put_bids)
                                     ), axis=0)

        self.smoothing = 1.0e-10 #Numerical parameter 

    def call_put_premium(self, k: float, kernel_width: float):

        call_len = len(self.opt_slice.call.strikes)
        put_len = len(self.opt_slice.put.strikes)
        call_ones = np.array([1.0] * call_len + [0.0] * put_len)
        put_ones = np.array([0.0] * call_len + [1.0] * put_len)

        def convex(z):
            a, b, c = (-1.5, 0.0, 1.5)

            res = (max(0, z - a)**3 - max(0, z - b)**3) / (b - a)
            res -= (max(0, z - b)**3 - max(0, z - c)**3) / (c - b)
            res -= (b - a)**2
                
            res += (max(0, c - z)**3 - max(0, b - z)**3) / (c - b)
            res -= (max(0, b - z)**3 - max(0, a - z)**3) / (b - a)
            res -= (c-b)**2
            return res / 12.0        

        basis = np.column_stack([call_ones, 
                                 put_ones,
                                 self.ks - k, 
                                 np.vectorize(convex)((self.ks / k - 1.0) / kernel_width)])

        ws = np.sqrt(np.exp(-0.5 * ((self.ks / k - 1.0) / kernel_width) **2) / self.put_err)
        b_ws = np.column_stack([ws] * 4) 

        var = np.matmul((basis * b_ws).T, basis * b_ws)
        var_regul = self.smoothing * var.trace() / float(var.shape[0]) * np.identity(var.shape[0])
        var_regul[0, 0] = 0.0
        var_regul[1, 1] = 0.0
        var += var_regul

        cov = (basis * b_ws).T.dot(self.put_target * ws)
        coeffs = np.linalg.solve(var, cov)
        return (coeffs[0] - k), coeffs[1]


def fit_forward_curve(quotes: OptionSnapshot, 
                      box_spread: float) -> Dict[dt.datetime, float]:

    pricing_dt = quotes.time_stamp

    previous_fit_fwd = quotes.ref_spot
    forwards = {}
    for quote_sl in quotes.slices:        

        if (quote_sl.expiry < pricing_dt):
            continue

        t = act365_time(pricing_dt, quote_sl.expiry)
        box_discount = quote_sl.discount * math.exp(-box_spread * t)
        quote_sl = OptionQuoteSlice(quote_sl.symbol,
                                    quote_sl.expiry,
                                    box_discount,
                                    quote_sl.call,
                                    quote_sl.put)
        quote_sl = prepare_quotes_for_fit(quote_sl,
                                          pricing_dt,
                                          previous_fit_fwd,
                                          yield_threshold=0.01)
        
        # Raw estimation of atm spot volatility
        DEFAULT_VOL = 0.15
        opt_kernel = OptionKernelRegression(quote_sl)
        c, p = opt_kernel.call_put_premium(quotes.ref_spot, 0.5 * DEFAULT_VOL * np.sqrt(t))
        raw_fwd = c - p + quotes.ref_spot
        raw_vol = 0.5 * (bs_implied_volatility(raw_fwd, quotes.ref_spot, p, t, -1.0)
                         + bs_implied_volatility(raw_fwd, quotes.ref_spot, c, t, 1.0))

        dev = raw_vol * np.sqrt(t)
        KERNEL_NB_DEV_WIDTH = 0.10
        kernel_width = KERNEL_NB_DEV_WIDTH * dev
        ks = raw_fwd * np.exp(dev * np.linspace(-0.5, 0.5, 11))
        estimated_fwds = []
        for k in ks:
            c, p = opt_kernel.call_put_premium(k, kernel_width)
            estimated_fwds.append(c - p + k)
        forward = np.median(estimated_fwds)

        forwards[quote_sl.expiry] = forward
        previous_fit_fwd = forward

    return forwards
