import datetime as dt


def act365_time(t0: dt.datetime, t1: dt.datetime):
    delta = t1 - t0
    return (delta.days + delta.seconds / 86400.0) / 365.0
