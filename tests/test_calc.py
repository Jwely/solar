import unittest
import numpy as np
from datetime import datetime, timedelta

from solar import SolarCalculator
from tests import REFERENCE_LAT, REFERENCE_LON, REFERENCE_DATE


class TestCalculator(unittest.TestCase):

    def test_single_point_and_time(self):
        """
        test solar calculator for single point on earths surface at single time.
        Results must be manually verified!
        """

        s = SolarCalculator(lat=REFERENCE_LAT, lon=REFERENCE_LON, dt=REFERENCE_DATE)
        # todo

    def test_multipoint(self):
        """
        test solar calculator for grid of points on earth surface at single time.
        Results must be manually verified!
        """

        dt = REFERENCE_DATE

        lat = [REFERENCE_LAT + i / 10 for i in range(-2, 1)]
        lon = [REFERENCE_LON + i / 10 for i in range(-2, 1)]
        lat, lon = np.meshgrid(lat, lon)

        s = SolarCalculator(lat=lat, lon=lon, dt=dt)
        # todo

    def test_multitime(self):
        """
        test solar calculator for single point on earths surface at multiple times.
        Results must be manually verified!
        """

        dt = [REFERENCE_DATE +
              timedelta(hours=hour)
              for hour in range(24)]

        s = SolarCalculator(lat=REFERENCE_LAT, lon=REFERENCE_LON, dt=dt)
        # todo

    def test_multitime_multipoint(self):
        """
        test solar calculator for multiple points on earths surface at multiple times.
        Results must be manually verified!
        """
        # example times on todays date at - GMT (tested in east cost USA during DST)
        dt = [REFERENCE_DATE +
              timedelta(hours=hour)
              for hour in range(24)]

        lat = [REFERENCE_LAT + i / 10 for i in range(-2, 3)]
        lon = [REFERENCE_LON + i / 10 for i in range(-2, 3)]
        lat, lon = np.meshgrid(lat, lon)

        #lat = np.repeat(lat[:, :, np.newaxis], len(dt), axis=2)
        #lon = np.repeat(lon[:, :, np.newaxis], len(dt), axis=2)

        s = SolarCalculator(lat=lat, lon=lon, dt=dt)
        # todo

    def test_high_compute(self):
        # TODO: more of a performance test, don't want to run something hefty on every test run.
        """ test computations for very large vector space """
        # do houry going back 50 years
        now = datetime.utcnow()
        dt = [REFERENCE_DATE +
              timedelta(hours=hour)
              for hour in range(24*365*50)]

        lat = [38.820450 + i / 10 for i in range(-10, 10)]
        lon = [-77.050552 + i / 10 for i in range(-10, 10)]
        lat, lon = np.meshgrid(lat, lon)

        lat = np.repeat(lat[:, :, np.newaxis], len(dt), axis=2)
        lon = np.repeat(lon[:, :, np.newaxis], len(dt), axis=2)

        s = SolarCalculator(lat=lat, lon=lon, dt=dt, low_mem=True)
        # todo
