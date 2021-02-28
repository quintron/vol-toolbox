import datetime as dt
from voltoolbox.calendar import nyse_calendar, nyse_holidays


def validate_count_open_days(cal, start, end):
    nb_opens = cal.count_open_days(start, end)

    count = 0
    d = start
    while d < end:
        if not cal.is_closed(d):
            count = count + 1
        d = d + dt.timedelta(days=1)

    assert nb_opens == count


def test_calendar_count_open_days():
    cal = nyse_calendar(from_year=2020, to_year=2030)
    for year in (2020, 2021, 2022, 2023, 2024, 2025, 2026):
        for month in (1, 5, 10, 12):
            start = dt.date(year, month, 1)
            for end_days in (7, 30, 250, 365):
                end = start + dt.timedelta(days=end_days)
                validate_count_open_days(cal, start, end)

def test_nyse_holidays_regression():
    holidays_2020 = nyse_holidays(from_year=2020, to_year=2020)
    ref_holidays_2020 = (dt.date(2020, 1, 1),
                         dt.date(2020, 1, 20),
                         dt.date(2020, 2, 17),
                         dt.date(2020, 4, 10),
                         dt.date(2020, 5, 25),
                         dt.date(2020, 7, 3),
                         dt.date(2020, 9, 7),
                         dt.date(2020, 11, 26),
                         dt.date(2020, 12, 25))
    for h, h_ref in zip(holidays_2020, ref_holidays_2020):
        assert h == h_ref

    holidays_2021 = nyse_holidays(from_year=2021, to_year=2021)
    ref_holidays_2021 = (dt.date(2021, 1, 1),
                         dt.date(2021, 1, 18),
                         dt.date(2021, 2, 15),
                         dt.date(2021, 4, 2),
                         dt.date(2021, 5, 31),
                         dt.date(2021, 7, 5),
                         dt.date(2021, 9, 6),
                         dt.date(2021, 11, 25),
                         dt.date(2021, 12, 24))
    for h, h_ref in zip(holidays_2021, ref_holidays_2021):
        assert h == h_ref
