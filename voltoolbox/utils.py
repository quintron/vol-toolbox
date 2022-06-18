import re
import pandas as pd


def act365_time(t0, t1):
    delta = t1 - t0
    return (delta.days + delta.seconds / 86400.0) / 365.0


NB_REG=re.compile('\d+')
UNIT_REG=re.compile('[hHdDmMyY]')
DUR_REG = re.compile(f'^{NB_REG.pattern}{UNIT_REG.pattern}$')


def parse_offset(s):
    if DUR_REG.search(s) :
        nb_unit = int(NB_REG.match(s).group(0))
        unit_desc = str.lower(UNIT_REG.findall(s)[0])
        if unit_desc=='h':
            return pd.DateOffset(hours=nb_unit)
        elif unit_desc=='d':
            return pd.DateOffset(days=nb_unit)
        elif unit_desc=='m':
            return pd.DateOffset(months=nb_unit)
        elif unit_desc=='y':
            return pd.DateOffset(years=nb_unit)
    raise ValueError(f'Unable to parse offset :{s}')
