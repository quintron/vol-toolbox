import datetime as dt
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from voltoolbox.utils import act365_time
from voltoolbox.fit.option_quotes import OptionSnapshot, OptionQuoteSlice


def forward_bidask(fwd_bids: pd.Series,
                   fwd_asks: pd.Series) -> Tuple[float, float]:

    def _clean_up_data(fwd_bids, fwd_asks, q: float):
        bid_mask = fwd_bids < np.quantile(fwd_asks, q)
        ask_mask = fwd_asks > np.quantile(fwd_bids, 1.0 - q)
        return fwd_bids.loc[bid_mask], fwd_asks.loc[ask_mask]

    def _adjust_box_spread(fwd_bids, fwd_asks):
        reg_res = stats.theilslopes(np.concatenate((fwd_asks.values, fwd_bids.values), axis=0),
                                    np.concatenate((fwd_asks.index, fwd_bids.index), axis=0), 0.90)
        box_spd_discount =(1.0 - reg_res[0])
        fwd_bids = (fwd_bids - fwd_bids.index) / box_spd_discount + fwd_bids.index
        fwd_asks = (fwd_asks - fwd_asks.index) / box_spd_discount + fwd_asks.index
        return fwd_bids, fwd_asks

    Q_THRESH = (0.90, 0.75, 0.5, 0.25, 0.10)
    for q in Q_THRESH:
        fwd_bids, fwd_asks = _clean_up_data(fwd_bids, fwd_asks, q)
        fwd_bids, fwd_asks = _adjust_box_spread(fwd_bids, fwd_asks)

    fwd_best_bid = np.quantile(fwd_bids, 1.0 - Q_THRESH[-1])
    fwd_best_ask = np.quantile(fwd_asks, Q_THRESH[-1])
    return fwd_best_bid, fwd_best_ask


def fit_forward(option_quote_sl: OptionQuoteSlice,
                pricing_dt: dt.datetime,
                box_spread_guess: float) -> Tuple[float, float]:
    time_to_mat = act365_time(pricing_dt, option_quote_sl.expiry)
    discount = option_quote_sl.discount * np.exp(-box_spread_guess * time_to_mat)

    call = option_quote_sl.call
    call_asks = pd.Series(call.asks, call.strikes, dtype=float)
    call_bids = pd.Series(call.bids, call.strikes, dtype=float)

    put = option_quote_sl.put
    put_asks = pd.Series(put.asks, put.strikes, dtype=float)
    put_bids = pd.Series(put.bids, put.strikes, dtype=float)

    quote_df = pd.DataFrame({'call_bids' : call_bids,
                             'call_asks' : call_asks,
                             'put_bids' : put_bids,
                             'put_asks' : put_asks})
    valid_quote_mask = quote_df.put_bids > 0.0
    valid_quote_mask &= quote_df.put_bids <= quote_df.put_asks
    valid_quote_mask &= quote_df.call_bids > 0.0
    valid_quote_mask &= quote_df.call_bids <= quote_df.call_asks
    quote_df = quote_df.loc[valid_quote_mask]

    # fwd without box spread discount adjustement
    fwd_asks = (quote_df.call_asks - quote_df.put_bids) / discount
    fwd_asks += fwd_asks.index
    fwd_bids = (quote_df.call_bids - quote_df.put_asks) / discount
    fwd_bids += fwd_bids.index

    return forward_bidask(fwd_bids, fwd_asks)


def fit_forward_curve(quotes: OptionSnapshot,
                      box_spread_guess: float) -> Dict[dt.datetime, Tuple[float, float]]:

    pricing_dt = quotes.time_stamp
    forwards = {}
    for quote_sl in quotes.slices:

        if quote_sl.expiry < pricing_dt:
            continue

        forwards[quote_sl.expiry] = fit_forward(quote_sl, pricing_dt, box_spread_guess)

    return forwards
