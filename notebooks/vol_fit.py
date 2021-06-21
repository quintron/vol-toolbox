import math
import datetime as dt
import numpy as np

from voltoolbox import BusinessTimeMeasure, bs_implied_volatility
from voltoolbox.fit.option_quotes import OptionQuoteSlice, VolQuoteSlice, VolSlice
from fit_utils import act365_time, filter_quotes


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

    put = quote_slice.put
    put_xs = []
    put_mids = []
    put_errs = []
    for k, b, a in zip(put.strikes, put.bids, put.asks):
        vol_bid = bs_implied_volatility(forward, k, b, time_to_maturity, -1.0)
        vol_ask = bs_implied_volatility(forward, k, a, time_to_maturity, -1.0)
        if vol_bid > 0.0 and vol_ask > vol_bid:
            put_mids.append(0.5 * (vol_ask + vol_bid))
            put_errs.append(0.5 * (vol_ask - vol_bid))
            put_xs.append(math.log(k / forward))
    put_vols = VolSlice(tuple(put_xs),
                        tuple(put_mids),
                        tuple(put_errs))

    call = quote_slice.call
    call_xs = []
    call_mids = []
    call_errs = []
    for k, b, a in zip(call.strikes, call.bids, call.asks):
        vol_bid = bs_implied_volatility(forward, k, b, time_to_maturity, 1.0)
        vol_ask = bs_implied_volatility(forward, k, a, time_to_maturity, 1.0)
        if vol_bid > 0.0 and vol_ask > vol_bid:
            call_mids.append(0.5 * (vol_ask + vol_bid))
            call_errs.append(0.5 * (vol_ask - vol_bid))
            call_xs.append(math.log(k / forward))
    call_vols = VolSlice(tuple(call_xs),
                        tuple(call_mids),
                        tuple(call_errs))

    return VolQuoteSlice(quote_slice.symbol, quote_slice.expiry, time_to_maturity, forward, call_vols, put_vols)
