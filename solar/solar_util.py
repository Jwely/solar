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

    # we want to return the full time series between the start_date and end_date boundaries
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
        slope: float,
        aspect: float,
        method_str: str,
        time_step: timedelta = None) -> pd.DataFrame:
    """
    Calculates a daily riemann sum for any one method of the solar calculator class.
    Designed for taking daily incident radiation totals.

    :param lat_lon_point: lat_lon point, as Point(lon, lat)
    :param start_date: first utc date to integrate daily local PAR for.
    :param end_date: last utc date to integrate daily local PAR for.
    :param method_str: the solar class method to integrate over the time specified.
    :param time_step: the step size of the riemann sum integral. defaults to 20 minutes.
    :return:
    """
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
    sc = SolarCalculator(lat=lat_lon_point.y, lon=lat_lon_point.x, dt=dts, slope=slope, aspect=aspect, low_mem=True)

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
    })
    return daily.reset_index()


def clear_day_daily_par_integral(
        lat_lon_point: Point,
        start_date: datetime,
        end_date: datetime,
        slope: float,
        aspect: float,
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
        slope=slope,
        aspect=aspect,
        method_str="surface_par",
        time_step=time_step)
    return df


def solar_energy_stc():
    """
    Replicates a scenario equal to the "standard testing condition" for solar power.
    which is an iiradiance of 1000W/m2 and a sun-facing (aspect = azimuth) surface
    and a solar elevation of 41.81 degrees. This is designed to represent equinoxes
    in the USA. (My opinion is that this is a shit standardization)
    """

    # This scenario created by visual guess and check iteration.

    # started with washington DC coordinate and moved north to get elevation correct.
    # date of the autumn equinox, in the morning at 10:02:40 am local time.
    # azimuth and aspect equivalent, standard 37 degree slope.
    # elevation ~41.81 degrees.
    sc = SolarCalculator(
        lat=40,
        lon=-77.0369,
        dt=datetime(2020, 9, 22, 15, 2, 40),
        slope=37,
        aspect=138.65987758,
        _air_mass_method="kasten-young",
        _pollution_condition="average")
    sc.summarize(verbose=True)


def _main():
    from typing import Dict, List
    from matplotlib import pyplot as plt
    import seaborn as sns
    import numpy as np
    import pandas as pd
    from pprint import pprint

    # input params
    cp = Point(-77.499069, 38.291597)
    slope = 25
    aspect = 175
    start_date = datetime(2019, 1, 1, 0, 10)
    end_date = datetime(2020, 1, 1, 0, 0)
    time_step = timedelta(minutes=10)

    dates = local_datetime_array(
        cp,
        start_date=start_date,
        end_date=end_date,
        time_step=time_step)

    sc = SolarCalculator(
        lat=cp.y,
        lon=cp.x,
        dt=dates,
        slope=slope,
        aspect=aspect,
        _pollution_condition="average")
    dates = localize_utc_datetimes_by_lat_lon(sc.dt(), point=cp)

    # set up the dataframe
    df = pd.DataFrame()
    df["utc"] = sc.dt()
    df["datetime"] = dates
    df["tzoffset"] = df["datetime"].apply(lambda x: x.utcoffset().total_seconds() / (24 * 60 * 60))
    df["date"] = df["datetime"].apply(lambda x: datetime(x.year, x.month, x.day))

    # set up reimann sum of incident energy
    df["surf_norm_boa_irradiance"] = sc.surf_norm_boa_irradiance() * time_step.total_seconds()
    df["zenith"] = sc.zenith()
    df["elevation"] = sc.elevation()
    df["sunrise"] = (sc.sunrise() + df["tzoffset"]) * 24
    df["sunset"] = (sc.sunset() + df["tzoffset"]) * 24
    df["solar_noon"] = (sc.solar_noon() + df["tzoffset"]) * 24
    df["sunlight"] = sc.sunlight() * 24
    df["sun_surf_angle"] = sc.sun_surf_angle()

    grps = df.groupby(by="date")

    daily = pd.DataFrame({
        "surf_norm_irradiance (kwh/m2*day)": grps["surf_norm_boa_irradiance"].sum() / 3.6e6,
        "max elevation (deg)": grps["elevation"].max(),
        "sunrise": grps["sunrise"].first(),
        "solar_noon": grps["solar_noon"].first(),
        "sunset": grps["sunset"].first(),
        "sunlight (hours)": grps["sunlight"].mean(),
        "min sun_surf_angle (deg)": grps["sun_surf_angle"].min(),
    }).reset_index()

    # panel based conversions
    panel_efficiency = 0.195
    conversion_efficiency = 0.925
    n_panels = 2
    panel_area = 1.67
    cloud_fraction = 0.40
    cloud_eff = 0.25
    cloud_factor = (1 - cloud_fraction) + (cloud_fraction * cloud_eff)

    daily["recoverable power (kwh/day)"] = daily["surf_norm_irradiance (kwh/m2*day)"] * \
        (cloud_factor * conversion_efficiency * panel_efficiency * n_panels * panel_area)

    fix, ax = plt.subplots(3, 1)

    sns.lineplot(x="date", y="sunset", data=daily, ax=ax[0], c="blue", label="sunset")
    sns.lineplot(x="date", y="solar_noon", data=daily, ax=ax[0], c="red", label="solar_noon")
    sns.lineplot(x="date", y="sunrise", data=daily, ax=ax[0], c="orange", label="sunrise")
    sns.lineplot(x="date", y="sunlight (hours)", data=daily, ax=ax[0], c='grey', label="sunlight")
    ax[0].set_ylabel("hours or hour of day")
    ax[0].set_yticks(range(0, 25, 2))
    ax[0].set_yticklabels([f"{x}:00" for x in range(0,25,2)])
    ax[0].grid()

    sns.lineplot(x="date", y="max elevation (deg)", data=daily, ax=ax[1], label="max elevation")
    sns.lineplot(x="date", y="min sun_surf_angle (deg)", data=daily, ax=ax[1], label="min sun_surf_angle")
    ax[1].grid()

    sns.lineplot(x="date", y="surf_norm_irradiance (kwh/m2*day)", data=daily,
                 ax=ax[2], c='blue', label="surface_norm_irradiance (kwh/m2*day)")
    sns.lineplot(x="date", y="recoverable power (kwh/day)", data=daily,
                 ax=ax[2], c='purple', label="recoverable power (kwh/day)")

    ax[2].text(x=dates[0], y=-2,
               s=f"panel eff={panel_efficiency}\n"
                 f"conversion eff={conversion_efficiency}\n"
                 f"area={n_panels*panel_area} m2\n"
                 f"cloud factor={cloud_factor}\n"
                 f"roof_slope={slope}\n"
                 f"aspect={aspect}")
    ax[2].grid()

    plt.show()

    daily.to_csv("solar_report.csv")
    print(daily.mean())
    print(daily.sum())


if __name__ == "__main__":
    #solar_energy_stc()
    _main()


