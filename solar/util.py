from typing import Union, Iterable
import numpy as np
from datetime import datetime, timedelta
from dateutil import tz

import timezonefinder


def vec_timedelta(**kwargs):
    """
    Just like datetimes timedelta method, but will return numpy arrays of
    timedelta objects for vectorized inputs. kwargs are exactly as those
    supported by timedelta, and all argument inputs may be either singular
    int, float, or ndarrays with dtype int or float.

    for example:

        td_array = vec_timedelta(days=np.array(range(1,366))
        returns a numpy array of timedelta objects for between 1 and 365
        days
    """

    # all arguments may be vectors, but must be the same shape!
    if len(kwargs) > 1:
        if not all([v1.shape == v2.shape for v1 in kwargs.values() for v2 in kwargs.values()]):
            raise Exception(
                "If vectorized arguments are used, they must all be "
                "identical shapes, got [{}]".format(
                    ",".join(
                        ["{}=({})".format(k, len(v)) for k, v in kwargs.items()]
                    )))

    # argument type conditioning
    is_vec = False
    for k, v in kwargs.items():
        if isinstance(v, np.ndarray):
            is_vec = True
            kwargs[k] = kwargs[k].astype(dtype=float)

    # produce vectorized version of timedelta.
    vf = np.vectorize(timedelta)

    # return the appropriate one
    if is_vec:
        return vf(**kwargs)
    else:
        return timedelta(**kwargs)


def lat_lon_to_timezone(lat: float, lon: float) -> str:
    """
    Looks up the time zone of a lat-lon.
    """
    tf = timezonefinder.TimezoneFinder()
    tz_str = tf.certain_timezone_at(lat=lat, lng=lon)
    return tz_str


def localize_datetimes_by_lat_lon(
        dt: Union[datetime, Iterable[datetime]],
        lat: float,
        lon: float)\
        -> Union[datetime, Iterable[datetime]]:
    """
    Assigns a time zone to a datetime or array-like of datetimes based on a Point coordinate.
    """

    lat_lon_tz_str = lat_lon_to_timezone(lat=lat, lon=lon)
    lat_lon_tz = tz.gettz(lat_lon_tz_str)

    def _do_one(dt: datetime):
        utc_aware_datetime = dt.replace(tzinfo=tz.tzutc())
        local_datetime = utc_aware_datetime.astimezone(lat_lon_tz)
        return local_datetime

    if hasattr(dt, "__iter__"):
        return [_do_one(ld) for ld in dt]
    else:
        return _do_one(dt)
