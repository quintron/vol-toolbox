import datetime as dt
import voltoolbox
from voltoolbox import check_datetime


def valid_datetime(t):
    assert check_datetime(t) == t


def test_datetime():
    valid_datetime(dt.datetime(2020, 2, 29, 23, 59, 59))
    valid_datetime(dt.datetime(2021, 2, 26, 18, 0, 15))
