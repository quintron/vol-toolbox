import datetime as dt
import numpy as np

from typing import Tuple
from voltoolbox import longest_increasing_subsequence
from voltoolbox.fit.option_quotes import OptionQuoteSlice, QuoteSlice


def act365_time(t0: dt.datetime, t1: dt.datetime):
    delta = t1 - t0
    return (delta.days + delta.seconds / 86400.0) / 365.0


def filter_quotes(quote_slice: OptionQuoteSlice,
                  fwd_bracket: Tuple[float, float],
                  *, discounted_output_quotes: bool=True
                 ) -> OptionQuoteSlice:
    """Return discounted option quotes, outlier have been removed.
    """
    fwd_min, fwd_max = fwd_bracket
    assert(fwd_max >= fwd_min)

    discount = quote_slice.discount
    out_discount = 1.0 if discounted_output_quotes else discount

    # Prepare call slice
    call_sl = quote_slice.call
    call_ks = np.array(call_sl.strikes)
    call_bid = np.array(call_sl.bids) / discount
    call_ask = np.array(call_sl.asks) / discount

    call_sanity_filter = (call_bid > 0.0) & (call_ask > call_bid)
    call_ks = call_ks[call_sanity_filter]
    call_bid = call_bid[call_sanity_filter]
    call_ask = call_ask[call_sanity_filter]
    call_mid = 0.5 * (call_ask + call_bid)

    put_intrinsic_filter = call_mid > np.maximum(0.0, fwd_min - call_ks)
    inc_subseq = longest_increasing_subsequence(-call_mid[put_intrinsic_filter])
    filt_call_sl = QuoteSlice(tuple(call_ks[put_intrinsic_filter][inc_subseq]),
                              tuple(out_discount * call_bid[put_intrinsic_filter][inc_subseq]),
                              tuple(out_discount * call_ask[put_intrinsic_filter][inc_subseq]))

    # Prepare put slice
    put_sl = quote_slice.put
    put_ks = np.array(put_sl.strikes)
    put_bid = np.array(put_sl.bids) / discount
    put_ask = np.array(put_sl.asks) / discount

    put_sanity_filter = (put_bid > 0.0) & (put_ask > put_bid)
    put_ks = put_ks[put_sanity_filter]
    put_bid = put_bid[put_sanity_filter]
    put_ask = put_ask[put_sanity_filter]

    put_mid = 0.5 * (put_ask + put_bid)

    call_intrinsic_filter = put_mid > np.maximum(0.0, put_ks - fwd_max)
    inc_subseq = longest_increasing_subsequence(put_mid[call_intrinsic_filter])
    filt_put_sl = QuoteSlice(tuple(put_ks[call_intrinsic_filter][inc_subseq]),
                             tuple(out_discount * put_bid[call_intrinsic_filter][inc_subseq]),
                             tuple(out_discount * put_ask[call_intrinsic_filter][inc_subseq]))

    return OptionQuoteSlice(quote_slice.symbol, 
                            quote_slice.expiry,
                            discount,
                            filt_call_sl, 
                            filt_put_sl)