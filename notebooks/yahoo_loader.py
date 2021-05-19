import requests
import pytz
import datetime as dt
import numpy as np

from bs4 import BeautifulSoup
from pandas_datareader.yahoo.options import Options
from voltoolbox.fit.option_quotes  import OptionSnapshot, OptionQuoteSlice, QuoteSlice


_SYMBOLS = {
    'spx' : '%5ESPX',
    'rut' : '%5ERUT',
    'goog': 'GOOG'
}


PM_SETTLEMENT = (dt.time(15, 0), pytz.timezone('US/Eastern'))

AM_SETTLEMENT = (dt.time(9, 30), pytz.timezone('US/Eastern'))

AM_SETTLEMENT_SYMBOLS = ('spx', )


def option_settlement_datetime(expiry: dt.date, symbol: str):  
    if str.lower(symbol) in AM_SETTLEMENT_SYMBOLS:
        exp_time, tz = AM_SETTLEMENT
    else:
        exp_time, tz = PM_SETTLEMENT

    exp_dt = dt.datetime.combine(expiry, exp_time, tz)
    return exp_dt.astimezone(pytz.UTC)


def _datetime64_to_datetime(dt64):
    ts = (dt64 - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
    return pytz.UTC.localize(dt.datetime.utcfromtimestamp(ts))


def snap_options(symbol) -> OptionSnapshot:
    yahoo_symb = _SYMBOLS.get(str.lower(symbol), None)
    if yahoo_symb is None:
        raise Exception(f'Unknown symbol {symbol}')

    opt = Options(yahoo_symb)
    option_yahoo_df = opt.get_all_data()
    # filter useful information
    option_yahoo_df = option_yahoo_df[['Bid', 'Ask', 'Vol', 'Open_Int', 'Underlying_Price', 'IV', 'Root', 'Quote_Time']]
    # filter useless quotes
    option_yahoo_df = option_yahoo_df[(option_yahoo_df.Bid > 0.0) | (option_yahoo_df.Ask > 0.0)]
    option_table = option_yahoo_df.groupby(['Expiry','Root', 'Type', 'Strike']).first()

    slices = []
    slices_keys = sorted({(exp, root) for exp, root, _, _ in option_table.index})
    for exp, root in slices_keys:

        calls = option_table.loc[(exp, root, 'call')]
        call_slice = QuoteSlice(tuple(calls.index.values),
                                tuple(calls.Bid.values),
                                tuple(calls.Ask.values))

        puts = option_table.loc[(exp, root, 'put')]
        put_slice = QuoteSlice(tuple(puts.index.values),
                               tuple(puts.Bid.values),
                               tuple(puts.Ask.values))

        expiry = option_settlement_datetime(exp.date(), root)
        slices.append(OptionQuoteSlice(root, expiry, call_slice, put_slice))

    ref_spot = float(option_table.Underlying_Price.median())
    time_stamp = _datetime64_to_datetime(max(t for t in option_table.Quote_Time.values))
    return OptionSnapshot(time_stamp, ref_spot, slices)


_TREASURY_URL = "http://www.treasury.gov/resource-center/data-chart-center/interest-rates/Pages/TextView.aspx?data=yield"


def snap_ois_rate_curve():

    r = requests.get(_TREASURY_URL)
    soup = BeautifulSoup(r.text, 'html.parser')

    table = soup.find("table", attrs={'class' : 't-chart'})
    rows = table.find_all('tr')
    lastrow = len(rows)-1
    cells = rows[lastrow].find_all("td")
    date = cells[0].get_text()
    m1 = float(cells[1].get_text())
    m3 = float(cells[2].get_text())
    m6 = float(cells[3].get_text())
    y1 = float(cells[4].get_text())
    y2 = float(cells[5].get_text())
    y3 = float(cells[6].get_text())
    y5 = float(cells[7].get_text())
    y7 = float(cells[8].get_text())
    y10 = float(cells[9].get_text())
    y20 = float(cells[10].get_text())
    y30 = float(cells[11].get_text())

    years = (1 / 12, 3 / 12, 6 / 12, 12 / 12, 24 / 12, 36 / 12, 60 / 12, 84 / 12, 120 / 12, 240 / 12, 360 / 12)
    rates = (m1 / 100, m3 / 100, m6 / 100, y1 / 100, y2 / 100, y3 / 100, y5 / 100, y7 / 100, y10 / 100, y20 / 100, y30 / 100)
    return dict(zip(years, rates))
