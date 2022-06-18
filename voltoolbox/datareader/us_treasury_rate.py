import datetime as dt
import pandas as pd


_TREASURY_URL = "https://home.treasury.gov/resource-center/data-chart-center/interest-rates/TextView?type=daily_treasury_yield_curve&field_tdr_date_value_month={MONTH}"


def scrap_us_rate_curve(d: dt.date) -> pd.Series:
    url = _TREASURY_URL.format(MONTH=d.strftime('%Y%m'))
    rate_table = pd.read_html(url)[0]
    rate_table.set_index('Date', inplace=True)
    rate_curve = rate_table.loc[d.strftime('%m/%d/%Y')]
    rate_curve = rate_curve[['1 Mo', '2 Mo', '3 Mo', '6 Mo', '1 Yr', '2 Yr', '3 Yr', '5 Yr', '7 Yr', '10 Yr', '20 Yr', '30 Yr']]
    rate_curve.index = ['1m', '2m', '3m', '6m', '1y', '2y', '3y', '5y', '7y', '10y', '20y', '30y']
    rate_curve /= 100.0
    rate_curve.name = d
    return rate_curve


def scrap_last_us_rate_curve(d: dt.date) -> pd.Series:
    res = None
    shift = 0
    while res is None and shift < 10:
        try:
            res = scrap_us_rate_curve(d + dt.timedelta(days=-shift))
        except Exception:
            res=None
        shift += 1
    if res is None:
        raise Exception(f'Unable to find a curve for {d}')
    return res
