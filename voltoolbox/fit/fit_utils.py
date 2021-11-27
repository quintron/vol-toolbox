import math
import datetime as dt
import numpy as np

from typing import Tuple
from voltoolbox import (BusinessTimeMeasure, 
                        bs_implied_volatility, longest_increasing_subsequence)
from voltoolbox.fit.option_quotes import OptionQuoteSlice, QuoteSlice, VolQuoteSlice, VolSlice


def prepare_vol_quotes(quote_slice: OptionQuoteSlice,
                       forward: float,
                       box_spread: float, 
                       pricing_dt: dt.datetime,
                       time_measure: BusinessTimeMeasure) -> VolQuoteSlice:
    """Return (implied) vol quote slice from option quote slice"""
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
