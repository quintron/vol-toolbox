import requests
import pandas as pd
from voltoolbox.datareader.scrap_utils import random_header


_OPTIONS_BASE_URL = "https://query1.finance.yahoo.com/" "v7/finance/options/"


_OPTION_FIELDS = (
    ('expiration', 'expiration', int),
    ('strike', 'strike', float),
    ('symbol', 'contractSymbol', str),
    ('contract_size', 'contractSize', str),
    ('currency', 'currency', str),
    ('bid', 'bid', float),
    ('ask', 'ask', float),
    ('last_price', 'lastPrice', float),
    ('last_trade_timestamp', 'lastTradeDate', int),
    ('volume', 'volume', int),
    ('open_interest', 'openInterest', int),
)


class YahooClient:

    def __init__(self):
        self.session = requests.Session()

    def _load_url(self, url: str):
        res = self.session.get(url, headers=random_header())
        return res.json()

    def snap_option_quotes(self, symb: str):
        symb_base_url = f'{_OPTIONS_BASE_URL}/{symb}'
        base_res = self._load_url(symb_base_url)
        expiration_timestamps = base_res["optionChain"]["result"][0]["expirationDates"]

        quotes_datas = []
        for ts in expiration_timestamps:
            expi_url = f'{symb_base_url}?date={ts}'
            expi_res = self._load_url(expi_url)

            for option in expi_res["optionChain"]["result"][0]["options"]:
                for typ in ('calls', 'puts'):
                    for q in option[typ]:
                        quote = [{'calls': 'call', 'puts': 'put'}[typ]]  # option_type
                        for _, rkey, ntype in _OPTION_FIELDS:
                            try:
                                val = ntype(q[rkey])
                            except (KeyError, ValueError):
                                val = {int : -1, float: float('nan'), str: ''}[ntype]
                            quote.append(val)
                        quotes_datas.append(quote)
        quote_df = pd.DataFrame(data=quotes_datas,
                                columns = ['option_type'] + [t[0] for t in _OPTION_FIELDS])

        for f in ('expiration', 'last_trade_timestamp'):
            quote_df[f] = pd.to_datetime(quote_df[f], unit='s', origin='unix')

        for f in ('option_type', 'contract_size', 'currency'):
            quote_df[f] = quote_df[f].astype('category')

        quote_df['symbol'] = quote_df['symbol'].astype(str)

        return quote_df
