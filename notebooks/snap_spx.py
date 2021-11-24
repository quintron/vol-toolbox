import time
import pytz
import datetime as dt
from yahoo_loader import snap_options


def snap_and_write_options(symb: str) -> None:
    print(f'snap {symb}')
    try:
        option_quotes = snap_options(symb)

        time_stamp = option_quotes.time_stamp.strftime('%Y%m%d_%H%M')
        file_name = f'vol_{symb}_{time_stamp}.json'
        with open(file_name, 'w') as f:
            f.write(option_quotes.to_json())

        print(f'{file_name}')
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

while True:
    if is_trading():
        snap_and_write_options('SPX')
    time.sleep(SNAP_FREQUENCY)
