import time
import pytz
import datetime as dt
import pandas as pd
import numpy as np
from voltoolbox.datareader.yahoo import YahooClient


_YAHOO_CLIENT = YahooClient()


_YAHOO_SYMBOLS = {
    'spx' : '%5ESPX',
    'rut' : '%5ERUT',
    'goog': 'GOOG'
}


def snap_and_write_options(symb: str) -> None:
    print(f'snap {symb}')
    try:
        obs_time = pd.Timestamp.now(tz='utc')

        yahoo_symb = _YAHOO_SYMBOLS.get(str.lower(symb), None)
        if yahoo_symb is None:
            raise Exception(f'Unknown symbol {symb}')        
        option_df = _YAHOO_CLIENT.snap_option_quotes(yahoo_symb)
        option_df['observation_time'] = obs_time
        option_df.set_index(['observation_time', 'option_type', 'expiration', 'strike'], inplace=True)

        file_time_stamp = obs_time.strftime('%Y%m%d_%H%M')
        file_name = f'option_{symb}_{file_time_stamp}.csv'
        option_df.to_csv(file_name)

        print(f'    {file_name}')
    except Exception as e:
        print('snap failed :')
        print(e)


TRADING_SESSION = (dt.time(9, 30), dt.time(16, 0), pytz.timezone('US/Eastern'))


def is_trading():
    start, end, tz = TRADING_SESSION
    now = pytz.UTC.localize(dt.datetime.utcnow())
    now = now.astimezone(tz)

    is_weekend = now.weekday() in (5, 6)
    if is_weekend:
        return False

    return now.time() >= start and now.time() <= end


SNAP_FREQUENCY = 60 * 15 # 15 minutes


RAN_GEN = np.random.default_rng()


while True:
    if is_trading():
        snap_and_write_options('SPX')

    # randomize frequency to trick yahoo server
    wait_time = int(SNAP_FREQUENCY * RAN_GEN.uniform(0.8, 1.2, 1)[0])
    time.sleep(wait_time)
