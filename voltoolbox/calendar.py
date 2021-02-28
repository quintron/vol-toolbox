import datetime as dt
from voltoolbox import Calendar

def _weekend_adjusted_date(d, adjust_sunday_only = False):
    """ return date adjusted to monday if sunday, or friday if saturday""" 
    wd = d.weekday()
    if wd==6:
        return d + dt.timedelta(days=1)
    elif wd==5 and not adjust_sunday_only:
        return d + dt.timedelta(days=-1)
    else:
        return d

def _first_january_holiday(year, adjust_sunday_only = False):
    """1st January adjusted to monday if sunday, or friday if saturday"""
    return _weekend_adjusted_date(dt.date(year, 1, 1)) 


def _christmas_holiday(year):
    """25th December adjusted to monday if sunday, or friday if saturday"""
    return _weekend_adjusted_date(dt.date(year, 12, 25))


def _independance_day_holiday(year):
    """4th July adjusted to monday if sunday, or friday if saturday"""
    return _weekend_adjusted_date(dt.date(year, 7, 4))


def _martin_luther_king_holiday(year):
    """Martin Luther King's birthday (third Monday in January)"""
    if year < 1983:
        return None
    for day in range(15, 22):
        d = dt.date(year, 1, day)
        if d.weekday()==0:
            return d
    raise Exception('Should never get there')


def _president_day(year):
    if year < 1971:
        return _weekend_adjusted_date(dt.date(year, 2, 22))
    else:
        for day in range(15, 22):
            d = dt.date(year, 2, day)
            if d.weekday()==0:
                return d
        raise Exception('Should never get there')


def _memorial_day(year):
    if year < 1971:
        return _weekend_adjusted_date(dt.date(year, 5, 30)) 
    else:
        # last Monday in May
        for day in range(25, 32):
            d = dt.date(year, 5, day)
            if d.weekday()==0:
                return d
        raise Exception('Should never get there')


def _labor_day(year):
    """First monday in September"""
    for day in range(1, 8):
        d = dt.date(year, 9, day)
        if d.weekday()==0:
            return d
    raise Exception('Should never get there')


def _thanksgiving_day(year):
    for day in range(22, 29):
        d = dt.date(year, 11, day)
        if d.weekday()==3:
            return d
    raise Exception('Should never get there')


def easter_monday(y):
    EasterMonday = [
             98,  90, 103,  95, 114, 106,  91, 111, 102,        # 1901-1909
             87, 107,  99,  83, 103,  95, 115,  99,  91, 111,   # 1910-1919
             96,  87, 107,  92, 112, 103,  95, 108, 100,  91,   # 1920-1929
            111,  96,  88, 107,  92, 112, 104,  88, 108, 100,   # 1930-1939
             85, 104,  96, 116, 101,  92, 112,  97,  89, 108,   # 1940-1949
            100,  85, 105,  96, 109, 101,  93, 112,  97,  89,   # 1950-1959
            109,  93, 113, 105,  90, 109, 101,  86, 106,  97,   # 1960-1969
             89, 102,  94, 113, 105,  90, 110, 101,  86, 106,   # 1970-1979
             98, 110, 102,  94, 114,  98,  90, 110,  95,  86,   # 1980-1989
            106,  91, 111, 102,  94, 107,  99,  90, 103,  95,   # 1990-1999
            115, 106,  91, 111, 103,  87, 107,  99,  84, 103,   # 2000-2009
             95, 115, 100,  91, 111,  96,  88, 107,  92, 112,   # 2010-2019
            104,  95, 108, 100,  92, 111,  96,  88, 108,  92,   # 2020-2029
            112, 104,  89, 108, 100,  85, 105,  96, 116, 101,   # 2030-2039
             93, 112,  97,  89, 109, 100,  85, 105,  97, 109,   # 2040-2049
            101,  93, 113,  97,  89, 109,  94, 113, 105,  90,   # 2050-2059
            110, 101,  86, 106,  98,  89, 102,  94, 114, 105,   # 2060-2069
             90, 110, 102,  86, 106,  98, 111, 102,  94, 114,   # 2070-2079
             99,  90, 110,  95,  87, 106,  91, 111, 103,  94,   # 2080-2089
            107,  99,  91, 103,  95, 115, 107,  91, 111, 103,   # 2090-2099
             88, 108, 100,  85, 105,  96, 109, 101,  93, 112,   # 2100-2109
             97,  89, 109,  93, 113, 105,  90, 109, 101,  86,   # 2110-2119
            106,  97,  89, 102,  94, 113, 105,  90, 110, 101,   # 2120-2129
             86, 106,  98, 110, 102,  94, 114,  98,  90, 110,   # 2130-2139
             95,  86, 106,  91, 111, 102,  94, 107,  99,  90,   # 2140-2149
            103,  95, 115, 106,  91, 111, 103,  87, 107,  99,   # 2150-2159
             84, 103,  95, 115, 100,  91, 111,  96,  88, 107,   # 2160-2169
             92, 112, 104,  95, 108, 100,  92, 111,  96,  88,   # 2170-2179
            108,  92, 112, 104,  89, 108, 100,  85, 105,  96,   # 2180-2189
            116, 101,  93, 112,  97,  89, 109, 100,  85, 105    # 2190-2199
        ]
    return EasterMonday[y - 1901]

def _good_friday_day(year):
    em = easter_monday(year)
    return dt.date(year -1 , 12, 31) + dt.timedelta(days=em - 3)


def nyse_holidays(from_year=2000, to_year=2100):
    holidays = []
    for year in range(from_year, to_year + 1):
        holidays +=[_first_january_holiday(year, True),
                    _president_day(year),
                    _good_friday_day(year),
                    _memorial_day(year),                    
                    _independance_day_holiday(year),
                    _labor_day(year),
                    _thanksgiving_day(year),
                    _christmas_holiday(year)
                   ]
        if year >= 1998:
            holidays.append(_martin_luther_king_holiday(year))

        # President election day
        if year <= 1968 or (year <= 1980 and year % 4 == 0):
            for day in range(1, 8):
                d = dt.date(year, 11, day)
                if d.weekday()==3:
                    holidays.append(d)
                    break

    return sorted(holidays)


def nyse_calendar(from_year=2000, to_year=2100):
    return Calendar(nyse_holidays(from_year=from_year, to_year=to_year))


if __name__ == "__main__":
    for d in nyse_holidays(2020, 2010):
        print(d)

