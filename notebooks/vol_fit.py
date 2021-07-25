import math
import datetime as dt
from typing import Union, Tuple
import numpy as np

from voltoolbox import BusinessTimeMeasure, bs_implied_volatility, longest_increasing_subsequence
from voltoolbox.fit.option_quotes import OptionQuoteSlice, QuoteSlice, VolQuoteSlice, VolSlice
from fit_utils import act365_time, filter_quotes


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
        
        if abs(q[0] - prev_q[0]) > 1.0e-5:
            merged_datas.append((q[0], q[1], q[2]))
        else:            
            if q[3]=='c':
                assert(prev_q[3]=='p')
                if q[0] < 0.0:
                    q0 = prev_q
                    q1 = q
                else:
                    q0 = q
                    q1 = prev_q                
            else :
                assert(prev_q[3]=='c')
                if q[0] < 0.0:
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
