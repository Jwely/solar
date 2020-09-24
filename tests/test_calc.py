import unittest
import numpy as np
from datetime import datetime, timedelta

from solar import SolarCalculator
from shapely.geometry import Point
from solar.util import lat_lon_to_timezone

# scalar values (somewhere in washington DC)
TEST_LAT = 38.9072
TEST_LON = -77.0369
tz = lat_lon_to_timezone(Point(TEST_LON, TEST_LAT))
TEST_DATE = datetime(2020, 1, 1, tzinfo=tz)

# vector tests
_test_lats = [TEST_LAT + i / 10 for i in range(-2, 3)]
_test_lons = [TEST_LON + i / 10 for i in range(-2, 3)]
TEST_LATS_GRID, TEST_LONS_GRID = np.meshgrid(_test_lats, _test_lons)
TEST_DATES = [TEST_DATE + timedelta(hours=hour) for hour in range(24)]


class TestCalculator(unittest.TestCase):

    def test_single_point_and_time(self):
        """
        test solar calculator for single point on earths surface at single time.
        Results must be manually verified!
        """

        s = SolarCalculator(lat=TEST_LAT, lon=TEST_LON, dt=TEST_DATE)
        # todo

    def test_multipoint(self):
        """
        test solar calculator for grid of points on earth surface at single time.
        Results must be manually verified!
        """
        s = SolarCalculator(lat=TEST_LATS_GRID, lon=TEST_LONS_GRID, dt=TEST_DATE)
        # todo

    def test_multitime(self):
        """
        test solar calculator for single point on earths surface at multiple times.
        Results must be manually verified!
        """

        dt = [TEST_DATE +
              timedelta(hours=hour)
              for hour in range(24)]

        s = SolarCalculator(lat=TEST_LAT, lon=TEST_LON, dt=dt)
        # todo

    def test_multitime_multipoint(self):
        """
        test solar calculator for multiple points on earths surface at multiple times.
        Results must be manually verified!
        """
        # example times on todays date at - GMT (tested in east cost USA during DST)
        s = SolarCalculator(lat=TEST_LATS_GRID, lon=TEST_LONS_GRID, dt=TEST_DATES)
        # todo

    def test_high_compute(self):
        # TODO: more of a performance test, don't want to run something hefty on every test run.
        """ test computations for very large vector space """
        # do houry going back 50 years
        now = datetime.utcnow()
        dt = [TEST_DATE +
              timedelta(hours=hour)
              for hour in range(24*365*50)]

        lat = [38.820450 + i / 10 for i in range(-10, 10)]
        lon = [-77.050552 + i / 10 for i in range(-10, 10)]
        lat, lon = np.meshgrid(lat, lon)

        s = SolarCalculator(lat=lat, lon=lon, dt=dt, low_mem=True)
        # todo
