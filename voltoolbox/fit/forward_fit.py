import datetime as dt
import numpy as np
import pandas as pd
from typing import Dict, Tuple
from scipy import optimize
from voltoolbox.fit.option_quotes import OptionSnapshot, OptionQuoteSlice
from voltoolbox.fit.fit_utils import act365_time


ForwardBidAsk = Tuple[float, Tuple[float, float]]


def fit_forward_from_bidask(fwd_bids: np.array, 
                            fwd_asks: np.array) -> ForwardBidAsk:
    fwd_ask_med = fwd_asks.median()
    fwd_bid_med = fwd_bids.median()

    fwd_asks = fwd_asks[fwd_asks > fwd_bid_med]  # outlier clean-up
    fwd_bids = fwd_bids[fwd_bids <  fwd_ask_med]  # outlier clean-up

    fwd_ask = fwd_asks.quantile(0.25)
    fwd_bid = fwd_bids.quantile(0.75)

    fwd_asks = fwd_asks[fwd_asks < fwd_ask + 10.0 * (fwd_ask - fwd_bid)] # outlier clean-up
    fwd_bids = fwd_bids[fwd_bids > fwd_bid - 10.0 * (fwd_ask - fwd_bid)] # outlier clean-up

    ref_spd =  max(0.5 * np.log(fwd_ask / fwd_bid), 1.0 / 10000.0)
    def logit_loglikehood_deriv(f: float):
        fwd_bid_probas = 1.0 / (1.0 + np.exp(-np.log(fwd_bids / f) / ref_spd))
        fwd_ask_probas = 1.0 / (1.0 + np.exp(-np.log(fwd_asks / f) / ref_spd))
        return (-fwd_bid_probas).sum() + (1.0 - fwd_ask_probas).sum()
    fitted_fwd = optimize.brenth(logit_loglikehood_deriv, fwd_bid, fwd_ask)

    return fitted_fwd, (fwd_bid, fwd_ask)


def fit_forward(option_quote_sl: OptionQuoteSlice, 
                pricing_dt: dt.datetime,
                box_spread: float) -> ForwardBidAsk:
    time_to_mat = act365_time(pricing_dt, option_quote_sl.expiry)
    discount = option_quote_sl.discount * np.exp(-box_spread * time_to_mat)

    call = option_quote_sl.call
    call_asks = pd.Series(call.asks, call.strikes)
    call_bids = pd.Series(call.bids, call.strikes)

    put = option_quote_sl.put
    put_asks = pd.Series(put.asks, put.strikes)
    put_bids = pd.Series(put.bids, put.strikes)

    fwd_asks = (call_asks - put_bids) / discount
    fwd_asks += fwd_asks.index
    fwd_asks = fwd_asks.dropna()

    fwd_bids = (call_bids - put_asks) / discount
    fwd_bids += fwd_bids.index 
    fwd_bids = fwd_bids.dropna()

    return fit_forward_from_bidask(fwd_bids, fwd_asks)


def fit_forward_curve(quotes: OptionSnapshot,
                      box_spread: float) -> Dict[dt.datetime, ForwardBidAsk]:

    pricing_dt = quotes.time_stamp
    forwards = {}
    for quote_sl in quotes.slices:

        if (quote_sl.expiry < pricing_dt):
            continue

        forwards[quote_sl.expiry] = fit_forward(quote_sl, pricing_dt, box_spread)

    return forwards
