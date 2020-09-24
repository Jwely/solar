from datetime import datetime, timedelta
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from shapely.geometry import Point

from solar.calculator import SolarCalculator
from solar.util import local_datetime_array, localize_utc_datetimes_by_lat_lon


def basic_plots():
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

    sns.lineplot(x="date", y="sunset", data=daily, ax=ax[0], label="sunset")
    sns.lineplot(x="date", y="solar_noon", data=daily, ax=ax[0], label="solar_noon")
    sns.lineplot(x="date", y="sunrise", data=daily, ax=ax[0], label="sunrise")
    sns.lineplot(x="date", y="sunlight (hours)", data=daily, ax=ax[0], label="sunlight")
    ax[0].set_ylabel("hours or hour of day")
    ax[0].set_yticks(range(0, 25, 2))
    ax[0].set_yticklabels([f"{x}:00" for x in range(0,25,2)])
    ax[0].grid()

    sns.lineplot(x="date", y="max elevation (deg)", data=daily, ax=ax[1], label="max elevation")
    sns.lineplot(x="date", y="min sun_surf_angle (deg)", data=daily, ax=ax[1], label="min sun_surf_angle")
    ax[1].grid()

    sns.lineplot(x="date", y="surf_norm_irradiance (kwh/m2*day)", data=daily,
                 ax=ax[2], label="surface_norm_irradiance (kwh/m2*day)")
    sns.lineplot(x="date", y="recoverable power (kwh/day)", data=daily,
                 ax=ax[2], label="recoverable power (kwh/day)")

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
    basic_plots()