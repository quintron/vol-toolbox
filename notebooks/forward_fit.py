import math
import numpy as np
import datetime as dt
from voltoolbox import longest_increasing_subsequence
from voltoolbox.fit.option_quotes import OptionQuoteSlice, QuoteSlice


def act365_time(t0: dt.datetime, t1: dt.datetime):
    delta = t1 - t0
    return (delta.days + delta.seconds / 86400.0) / 365.0


def filter_quotes(quote_slice: OptionQuoteSlice, 
                  pricing_date: dt.datetime,
                  discount : float,
                  spot_prior: float,
                  *, yield_threshold=0.10):

    time_to_mat = act365_time(pricing_date, quote_slice.expiry)
    ratio_min = min(math.exp(- time_to_mat * yield_threshold), (1.0 - yield_threshold))
    ratio_max = math.exp(time_to_mat * yield_threshold)
    
    call_sl = quote_slice.call
    call_ks = np.array(call_sl.strikes)
    call_bid = np.array(call_sl.bids)
    call_ask = np.array(call_sl.asks)
    call_mid = 0.5 * (call_ask + call_bid)

    put_plus_fwd = call_mid / discount + call_ks 
    filter_ = put_plus_fwd > spot_prior * ratio_min
    inc_subseq = longest_increasing_subsequence(put_plus_fwd[filter_])
    filt_call_sl = QuoteSlice(tuple(call_ks[filter_][inc_subseq]),
                              tuple(call_bid[filter_][inc_subseq]),
                              tuple(call_ask[filter_][inc_subseq]))


    put_sl = quote_slice.put
    put_ks = np.array(put_sl.strikes)
    put_bid = np.array(put_sl.bids)
    put_ask = np.array(put_sl.asks)
    put_mid = 0.5 * (put_ask + put_bid)
    
    filter_ = put_mid / discount - put_ks > - spot_prior * ratio_max
    inc_subseq = longest_increasing_subsequence(put_mid[filter_])
    filt_put_sl = QuoteSlice(tuple(put_ks[filter_][inc_subseq]),
                             tuple(put_bid[filter_][inc_subseq]),
                             tuple(put_ask[filter_][inc_subseq]))
    
    return OptionQuoteSlice(quote_slice.symbol, 
                            quote_slice.expiry,
                            filt_call_sl, 
                            filt_put_sl)
