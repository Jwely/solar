from datetime import datetime
from shapely.geometry import Point
from solar.util import lat_lon_to_timezone

REFERENCE_LAT = 38.9072
REFERENCE_LON = -77.0369
tz = lat_lon_to_timezone(Point(REFERENCE_LON, REFERENCE_LAT))
REFERENCE_DATE = datetime(2020, 1, 1, tzinfo=tz)