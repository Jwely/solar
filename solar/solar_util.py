from typing import Union, List, Iterable
from datetime import datetime, timedelta, timezone
from collections.abc import Iterable as IterableType

import timezonefinder
from dateutil import tz
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from shapely.geometry import Point

from solar.solar_cls import SolarCalculator


def local_datetime_array(
        lat_lon_point: Point,
        start_date: datetime,
        end_date: datetime,
        time_step: timedelta):
    """
    Helper utility to instantiate a SolarCalculator for a single point in space, but many
    points in time.
    """

    # we want to return the full time series betwen the start_date and end_date boundaries
    # at local time, but the inputs should be in UTC. what we want is the exact UTC time
    # that marks the beginning and end of the day at the specific coordinate we are worried about.
    def _offset_utc(date, point):
        local_date = localize_utc_datetimes_by_lat_lon(date, point)
        # here we pretend localized date is in UTC again so we can compute the difference.
        offset = date.replace(tzinfo=timezone.utc) - local_date.replace(tzinfo=timezone.utc)
        return date + offset

    start_date = _offset_utc(start_date, lat_lon_point)
    end_date = _offset_utc(end_date, lat_lon_point)

    #  build a set of date samples. Now all in the right UTC range to perfectly cover the local date range.
    span = end_date - start_date
    num_steps = int(span / time_step)
    dts = [start_date + (time_step * i) for i in range(num_steps)]
    return dts


def lat_lon_to_timezone(point: Point) -> tz.tz.tzfile:
    """
    Looks up the time zone of a lat-lon and returns a timezone instance
    """
    lat = point.y
    lon = point.x

    tf = timezonefinder.TimezoneFinder()
    tz_str = tf.certain_timezone_at(lat=lat, lng=lon)
    return tz.gettz(tz_str)


def localize_utc_datetimes_by_lat_lon(
        x: Union[datetime, Iterable[datetime]],
        point: Point) -> Union[datetime, Iterable[datetime]]:
    """
    Input a timezone-naive datetime THAT IS IN UTC, and a corresponding
    coordinate (lon,lat) where that time is desired to be converted to local time.

    This would be identical to writing
        my_date: datetime = ???
        point: Point = ???
        local_tz = lat_lon_to_timezone(point)
        utc_datetime = my_date.replace(tzinfo=timezone.utc)
        local_datetime = utc_datetime.astimezone(local-tz)

    This function will return the same datetime, but localized
    :param x: single datetime,
    :param point:
    :return:
    """

    lat_lon_tz = lat_lon_to_timezone(point)

    def _localize_datetime_by_lat_lon(_x):
        if _x.tzinfo is not None and _x.tzinfo.utcoffset(_x) is not None:
            # this datetime is tzaware.
            if str(_x.tzinfo) == str(timezone.utc):
                utc_datetime = _x
            else:
                raise AttributeError(f"datetime {_x} is already timezone-aware and is {_x.tzinfo}, not UTC!")
        else:
            utc_datetime = _x.replace(tzinfo=timezone.utc)

        local_datetime = utc_datetime.astimezone(lat_lon_tz)
        return local_datetime

    if isinstance(x, IterableType):
        return [_localize_datetime_by_lat_lon(dt) for dt in x]
    else:
        return _localize_datetime_by_lat_lon(x)


def solar_calculator_daily_riemann_sum_(
        lat_lon_point: Point,
        start_date: datetime,
        end_date: datetime,
        method_str: str,
        time_step: timedelta = None) -> pd.DataFrame:
    """
    Calculates a daily riemann sum for any one method of the solar calculator class.
    Designed for "surface_par", but can easily work for any other

    :param lat_lon_point: lat_lon point, as Point(lon, lat) NOTE x,y convention!
    :param start_date: first utc date to integrate daily local PAR for.
    :param end_date: last utc date to integrate daily local PAR for.
    :param method_str: the solar class method to integrate over the time specified.
    :param time_step: the step size of the riemann sum integral. defaults to 20 minutes.
    :return:
    """
    # TODO: the version of this function to accumulate broadband

    if time_step is None:
        # 20 minutes shows < 0.05% daily sum variance from 20 second sampling in riemann integral.
        # increasing it from here rapidly increases approximation error.
        time_step = timedelta(minutes=20)

    dts = local_datetime_array(
        lat_lon_point=lat_lon_point,
        start_date=start_date,
        end_date=end_date,
        time_step=time_step)

    # low mem is desireable in this case since we're grabbing one variable then dumping it.
    sc = SolarCalculator(lat=lat_lon_point.y, lon=lat_lon_point.x, dt=dts, low_mem=True)

    # take a riemann sum "rectangle integrals"
    method = getattr(sc, method_str)    # get the method of the SolarCalculator (like surface_par)
    par_series = method()
    par_areas = par_series * time_step.total_seconds()  # integrated time dimension.

    # localize the date array back to local time
    dts = localize_utc_datetimes_by_lat_lon(dts, point=lat_lon_point)

    # aggregate to daily, check, and return the results.
    df = pd.DataFrame()
    df["datetime"] = dts
    df[method_str] = par_areas
    df["date"] = df["datetime"].apply(lambda x: datetime(x.year, x.month, x.day))

    # take the integral and the number of points in each.
    grps = df[["date", method_str]].groupby(by="date")
    daily = pd.DataFrame({
        method_str: grps[method_str].sum(),
        # "count": grps["dpi"].count()   # count was used to ensure
    })
    return daily.reset_index()


def clear_day_daily_par_integral(
        lat_lon_point: Point,
        start_date: datetime,
        end_date: datetime,
        time_step: timedelta = None) -> pd.DataFrame:
    """
    Wrapper for the solar calculator to compute surface_par at a sufficient time sampling
    to approximately integrate function over an entire day.

    daily_par_integral has units of joules/meter, representing the total photosynthetically
    active energy to fall on a patch of earth on the given date if the weather is totally clear.
    (no amospheric effect)


    Example:
    ```
        fairbanks_ak = Point(-147.7, 64.8)          # just south of the arctic circle (66.5 degrees)
        buenos_aires_arg = Point(-58.43, -34.6)     # southern hemisphere
        fredericksburg_va = Point(-77.5, 38.3)      # norhtern hemisphere

        dpari = clear_day_daily_par_integral(
            lat_lon_point=fredericksburg_va,
            start_date=datetime(2000, 1, 1),
            end_date=datetime(2020, 1, 1))

        # it takes much longer to plot the data than it does to compute it.
        sns.lineplot("date", "surface_par", data=dpari)
        plt.show()
    ```

    :param lat_lon_point: lat_lon point, as Point(lon, lat) NOTE x,y convention!
    :param start_date: first utc date to integrate daily local PAR for.
    :param end_date: last utc date to integrate daily local PAR for.
    :param time_step: the step size of the riemann sum integral. defaults to 20 minutes.

    :returns:
        A dataframe with columns ["date", "surface_par"] for Daily PAR Integral
    """

    df = solar_calculator_daily_riemann_sum_(
        lat_lon_point=lat_lon_point,
        start_date=start_date,
        end_date=end_date,
        method_str="surface_par",
        time_step=time_step)
    return df


