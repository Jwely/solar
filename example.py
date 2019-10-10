from datetime import datetime, timedelta
from solar.calculator import SolarCalculator
from solar.util import localize_datetimes_by_lat_lon


# An example for precisely calculating sunrise time for a date and coordinate.
def get_sunrise(utc_date: datetime, lat: float, lon: float):

    sc = SolarCalculator(lat=lat, lon=lon, dt=utc_date)
    sunrise = sc.sunrise_time()

    # convert to datetime and localize
    sunrise_utc = datetime(utc_date.year, utc_date.month, utc_date.day) + sunrise
    sunrise_local = localize_datetimes_by_lat_lon(sunrise_utc, lat=lat, lon=lon)
    return sunrise_local


# calculating the precise time of sunset and sunrise
lat = 38.303693
lon = -77.549559
s = get_sunrise(datetime.utcnow(), lat=lat, lon=lon)
print(s)