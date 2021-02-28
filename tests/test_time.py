import datetime as dt
import voltoolbox
from voltoolbox import check_datetime, check_date 


def validate_datetime(t):
    assert check_datetime(t) == t
    assert check_date(t.date()) == t.date()


def test_datetime():
    validate_datetime(dt.datetime(2020, 2, 29, 23, 59, 59))
    validate_datetime(dt.datetime(2021, 2, 26, 18, 0, 15))
