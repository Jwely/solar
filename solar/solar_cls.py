import numpy as np
from datetime import timedelta
from typing import Tuple, Dict
import pandas as pd

from solar import solar_fn as sfn
from solar import FlexNum, FlexDate


def _vec_timedelta(**kwargs):
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


def store_value(store: Dict):  # decorator factory
    """
    Decorator to intercept class method calls. If this method has already been
    called, it's result will be stored in the store dictionary. On subsequent
    calls the same result will be returned.

    A special access count key is used to store a dictionary of the number of
    times each of the keys has been accessed. A curious tidbit for performance
    profiling.

    NOTE: the store_value is in a separate dictionary, not directly
    against class attributes.

    TODO WARNING: This obviously comes with great assumptions of input param immutability.
        Any attributed stored in the store_value dict is inherently inconvenient to access
        and is *unlikely* to be mutated by accident. However, measures must be taken to prevent
        mutation of __init__ params.

    """
    access_counts = "__access_counts"
    if access_counts not in store.keys():
        store[access_counts] = dict()

    def decorator(func):
        def wrapper(*args, **kwargs):
            vs_name = f"{func.__name__}"
            if vs_name in store.keys():
                # return existing value
                result = store[vs_name]
                store[access_counts][vs_name] += 1
            else:
                # compute , store, and return
                result = func(*args, **kwargs)
                store[vs_name] = result
                store[access_counts][vs_name] = 1
            return result
        return wrapper
    return decorator


class SolarCalculator(object):
    """
    Handles common solar calculations commonly needed for energy balance equation work,
    from evapotranspirative and photosynthesis modelling to photovoltaic performance
    profiling.

    Intended to scale moderately and flexibly for high dimensionality problems. Any
    combination scalar, 1d, or 2d space parameters and scalar or 1d time parameters may
    be used: for up to 3 dimensions of simultaneous computation.

    Time dimensionality:
        scalar:
            a single reference datetime input 'dt'.
        1d:
            a list of datetime values covering multiple times. smaller time increments allow
            integration of variables to daily values commonly desired.

    space dimensionality:
        scalar:
            single point on the earths surface. lat, lon, slope, aspect are all scalars.
        1d:
            performant handling of list-like inputs. lat, lon, slope and aspect are all
            1d numpy arrays of the same length.
        2d:
            meshgrid handling well suited for combination of final outputs with raster
            data structures. For example, resampling a DEM to match a spectral data source
            like landsat or sentinel and extracting the lat/lon pairs of that DEM allows
            full specification of the lat, lon, slope, and aspect inputs for the full raster space.

    TODO: for some complex partial solar tracking functions, a 3d input for slope and aspect
            may be desired when using 2 spatial dims and 1 time dim. Presently 2d spatial
            inputs are cast to 3d arrays by copying into the 3rd dimension.

    Table of variable descriptions

        lat                     decimal degrees latitude (float OR numpy array)
        lon                     decimal degrees longitude (float OR numpy array)
        reference_datetime      datetime at UTC
        slope                   slope of land at lat,lon (0 is flat, 90 is vertical)
        aspect                  aspect of land at lat,lon (0 = N, 90 = E, 180 = S, 270 = W)
        -------------------------------------------------------------------------------
        lat                     latitude                                    (space)
        lon                     longitude                                   (space)
        rdt                     reference datetime object (date_time_obj)   (time)
        ajd                     absolute julian day                         (time)
        ajc                     absolute julian century                     (time)
        geomean_long            geometric mean longitude of the sun         (time)
        geomean_anom            geometric mean longitude anomaly of the sun (time)
        earth_eccent            eccentricity of earths orbit                (time)
        sun_eq_of_center        the suns equation of center                 (time)
        true_long               true longitude of the sun                   (time)
        true_anom               true longitude anomaly of the sun           (time)
        app_long                the suns apparent longitude                 (time)
        oblique_mean_ellipse    earth oblique mean ellipse                  (time)
        oblique_corr            correction to earths oblique ellipse        (time)
        right_ascension         suns right ascension angle                  (time)
        declination             solar declination angle                     (time)
        equation_of_time        equation of time (minutes)                  (time)
        hour_angle_sunrise      the hour angle at sunrise                   (space, time)
        solar_noon              LST of solar noon                           (space, time)
        sunrise                 LST of sunrise time                         (space, time)
        sunset                  LST of sunset time                          (space, time)
        sunlight                LST fractional days of sunlight             (space, time)
        true_solar              LST for true solar time                     (space, time)
        hour_angle              total hour angle                            (space, time)
        zenith                  zenith angle                                (space, time)
        elevation               elevation angle                             (space, time)
        azimuth                 azimuthal angle                             (space, time)
        earth_distance          earths distance to sun in meters            (time)
        sun_norm_irradiance     incident solar energy at earth distance     (time)
        sun_vec_cartesian       cartesian vector pointing at sun            (space, time)
        surf_vec_cartesian      cartesian vector normal to surface          (space)
        sun_surf_angle          angle between sun_vec and surf_vec          (space, time)
        surf_norm_irradiance    irradiance on surface with surf_vec         (space, time)
        surf_norm_par           photosynthetically active surf_norm rad     (space, time)

    Units unless otherwise labeled.
        angle    = degrees
        distance = meters
        energy   = joules
        time     = datetime, timedeltas, or absolute julian time. ALL TIME is in UTC time zone!
        power    = watts
    """

    __result_store = dict()

    def __init__(self,
                 lat: FlexNum,
                 lon: FlexNum,
                 dt: FlexDate,
                 slope: FlexNum = None,
                 aspect: FlexNum = None,
                 _air_mass_method: str = None,
                 _pollution_condition: str = None):
        """
        Accepts vectorized positions for lat and lon, (use meshgrid) and/or list of
        multiple datetimes. Computations are fully vectorized

        'lat' and 'lon' must have the same dimensions:
            1d arrays: like a list of coordinate pairs
            2d arrays: should form a meshgrid.

        'dt' may be a single time or a 1d list of times.

        'slope' and 'aspect' must either be None, or have the same dimensions as lat/lon

        :param lat: decimal degrees latitude
        :param lon: decimal degrees longitude
        :param dt: reference datetime or a list of datetimes
        :param slope: decimal degrees of slope (0 = flat, 90 for vertical surface)
        :param aspect: decimal degrees of aspect (0 = N, 90 = E, 180 = S, 270 = W)

        """

        # check inputs for dimensional compatibility
        lat, lon, dt, vt, vs, slope, aspect, space_shape, time_shape = \
            self._validate_init(lat, lon, dt, slope, aspect)

        self.__time_is_vec = vt
        self.__space_is_vec = vs
        self.space_shape = space_shape
        self.time_shape = time_shape

        self._dt = dt    # list of datetimes.
        self._lat = lat  # degrees
        self._lon = lon  # degrees
        self._lat_rad = np.radians(lat)
        self._lon_rad = np.radians(lon)
        self._slope = slope     # degrees
        self._aspect = aspect   # degrees
        self._ajd, self._ajc = sfn.get_absolute_julian_day_and_century(dt)

        # solar power configuration params
        self._air_mass_method = _air_mass_method
        self._pollution_condition = _pollution_condition


# = input checks and properties of the class
    @staticmethod
    def _validate_init(lat: FlexNum, lon: FlexNum, dt: FlexDate, slope: FlexNum, aspect: FlexNum) \
            -> Tuple[FlexNum, FlexNum, FlexDate, bool, bool, FlexNum, FlexNum, Tuple, Tuple]:
        """
        Ensures that input lat, lon, and dt arguments are dimensionally compatible.

        latitude and longitude should be scalars or form a meshgrid.
        """

        # does this solver vectorize calculations over time?
        if isinstance(dt, list):
            vec_time = True
            time_shape = len(dt),
        else:
            vec_time = False
            time_shape = tuple()

        if isinstance(lat, np.ndarray) and isinstance(lon, np.ndarray):
            vec_space = True
            if lat.shape != lon.shape:
                raise AttributeError("Non scalar (lat,lon) inputs must "
                                     "form a meshgrid such as: `lat, lon = np.meshgrid(lat, lon)`")
            else:
                space_shape = lat.shape

            if slope is None and aspect is None:
                slope = np.zeros(space_shape)
                aspect = np.zeros(space_shape)

            if slope.shape != space_shape or aspect.shape != space_shape:
                raise AttributeError("'slope', 'aspect', 'lat', and 'lon' must all have the same shape!")

            # if both space and time are vectorized, add 3rd time dimension to lat/lon
            # TODO; this is inefficient: might be consequential for large spatial arrays.
            if vec_time:
                if len(space_shape) == 2:
                    # make spatial arrays 3 dimensions for compatibility with dt
                    lat = np.repeat(lat[:, :, np.newaxis], len(dt), axis=2)
                    lon = np.repeat(lon[:, :, np.newaxis], len(dt), axis=2)
                    slope = np.repeat(slope[:, :, np.newaxis], len(dt), axis=2)
                    aspect = np.repeat(aspect[:, :, np.newaxis], len(dt), axis=2)
        else:
            vec_space = False
            space_shape = tuple()

        return lat, lon, dt, vec_time, vec_space, slope, aspect, space_shape, time_shape

    @property
    def _result_store(self):
        return self.__result_store

    @property
    def time_is_vec(self):
        return self.__time_is_vec

    @property
    def space_is_vec(self):
        return self.__space_is_vec

# = set on init
    def lat(self):
        """ degrees lat as user input """
        return self._lat

    def lon(self):
        """ degrees lon as user input """
        return self._lon

    def lat_rad(self):
        """ lat in radians """
        return self._lat_rad

    def lon_rad(self):
        """ lon in radians """
        return self._lon_rad

    def slope(self):
        return self._slope

    def aspect(self):
        return self._aspect

    def ajc(self):
        return self._ajc

    def ajd(self):
        return self._ajd

    def dt(self):
        return self._dt

# = Use value store to allow on the fly compute
    @store_value(__result_store)
    def geomean_long(self):
        return sfn.get_geomean_long(
            ajc=self.ajc())

    @store_value(__result_store)
    def geomean_anom(self):
        return sfn.get_geomean_anom(
            ajc=self.ajc())

    @store_value(__result_store)
    def earth_eccent(self):
        return sfn.get_earth_eccent(
            ajc=self.ajc())

    @store_value(__result_store)
    def sun_eq_of_center(self):
        return sfn.get_sun_eq_of_center(
            ajc=self.ajc(),
            geomean_anom=self.geomean_anom())

    @store_value(__result_store)
    def true_long(self):
        return sfn.get_true_long(
            geomean_long=self.geomean_long(),
            sun_eq_of_center=self.sun_eq_of_center())

    @store_value(__result_store)
    def true_anom(self):
        return sfn.get_true_anom(
            geomean_anom=self.geomean_anom(),
            sun_eq_of_center=self.sun_eq_of_center())

    @store_value(__result_store)
    def rad_vector(self):
        return sfn.get_rad_vector(
            earth_eccent=self.earth_eccent(),
            true_anom=self.true_anom())

    @store_value(__result_store)
    def app_long(self):
        return sfn.get_app_long(
            true_long=self.true_long(),
            ajc=self.ajc())

    @store_value(__result_store)
    def oblique_mean_ellipse(self):
        return sfn.get_oblique_mean_ellipse(
            ajc=self.ajc())

    @store_value(__result_store)
    def oblique_corr(self):
        return sfn.get_oblique_corr(
            ajc=self.ajc(),
            oblique_mean_ellipse=self.oblique_mean_ellipse())

    @store_value(__result_store)
    def right_ascension(self):
        return sfn.get_right_ascension(
            app_long=self.app_long(),
            oblique_corr=self.oblique_corr())

    @store_value(__result_store)
    def declination(self):
        return sfn.get_declination(
            app_long=self.app_long(),
            oblique_corr=self.oblique_corr())

    @store_value(__result_store)
    def equation_of_time(self):
        return sfn.get_equation_of_time(
            oblique_corr=self.oblique_corr(),
            geomean_long=self.geomean_long(),
            geomean_anom=self.geomean_anom(),
            earth_eccent=self.earth_eccent())

    @store_value(__result_store)
    def hour_angle_sunrise(self):
        return sfn.get_hour_angle_sunrise(
            declination=self.declination(),
            lat_r=self.lat_rad())

    @store_value(__result_store)
    def solar_noon(self):
        return sfn.get_solar_noon(
            lon=self.lon(),
            equation_of_time=self.equation_of_time())

    @store_value(__result_store)
    def sunrise(self):
        return sfn.get_sunrise(
            solar_noon=self.solar_noon(),
            hour_angle_sunrise=self.hour_angle_sunrise())

    @store_value(__result_store)
    def sunset(self):
        return sfn.get_sunset(
            solar_noon=self.solar_noon(),
            hour_angle_sunrise=self.hour_angle_sunrise())

    @store_value(__result_store)
    def sunlight(self):
        return sfn.get_sunlight(
            hour_angle_sunrise=self.hour_angle_sunrise())

    @store_value(__result_store)
    def true_solar(self):
        return sfn.get_true_solar(
            lon=self.lon(),
            equation_of_time=self.equation_of_time(),
            reference_datetime=self.dt())

    @store_value(__result_store)
    def hour_angle(self):
        return sfn.get_hour_angle(
            true_solar=self.true_solar())

    @store_value(__result_store)
    def zenith(self):
        return sfn.get_zenith(
            declination=self.declination(),
            hour_angle=self.hour_angle(),
            lat_r=self.lat_rad())

    @store_value(__result_store)
    def air_mass(self):
        return sfn.get_air_mass(
            zenith=self.zenith(),
            method=self._air_mass_method)

    @store_value(__result_store)
    def elevation(self):
        return sfn.get_elevation(
            zenith=self.zenith())

    @store_value(__result_store)
    def elevation_noatmo(self):
        return sfn.get_elevation_noatmo(
            zenith=self.zenith())

    @store_value(__result_store)
    def azimuth(self):
        return sfn.get_azimuth(
            lat_r=self.lat_rad(),
            declination=self.declination(),
            hour_angle=self.hour_angle(),
            zenith=self.zenith())

    @store_value(__result_store)
    def earth_distance(self):
        return sfn.get_earth_distance(
            rad_vector=self.rad_vector())

    @store_value(__result_store)
    def sun_norm_toa_irradiance(self):
        return sfn.get_sun_norm_toa_irradiance(
            earth_distance=self.earth_distance())

    @store_value(__result_store)
    def sun_norm_boa_irradiance(self):
        return sfn.air_mass_correction_(
            _toa_irradiance=self.sun_norm_toa_irradiance(),
            air_mass=self.air_mass(),
            pollution_condition=self._pollution_condition)

    @store_value(__result_store)
    def sun_vec_cartesian(self):
        return sfn.get_sun_vec_cartesian(
            azimuth=self.azimuth(),
            elevation=self.elevation())

    @store_value(__result_store)
    def surf_vec_cartesian(self):
        return sfn.get_surf_vec_cartesian(
            aspect=self.aspect(),
            slope=self.slope())

    @store_value(__result_store)
    def sun_surf_angle(self):
        return sfn.get_sun_surf_angle(
            sun_vec=self.sun_vec_cartesian(),
            surface_vec=self.surf_vec_cartesian())

    @store_value(__result_store)
    def surf_norm_toa_irradiance(self):
        return sfn.get_surf_norm_toa_irradiance(
            sun_surf_angle=self.sun_surf_angle(),
            sun_norm_irradiance=self.sun_norm_toa_irradiance())

    @store_value(__result_store)
    def surf_norm_boa_irradiance(self):
        return sfn.air_mass_correction_(
            _toa_irradiance=self.surf_norm_toa_irradiance(),
            air_mass=self.air_mass(),
            pollution_condition=self._pollution_condition)

    @store_value(__result_store)
    def surf_norm_toa_par(self):
        return sfn.get_surf_norm_toa_par(self.surf_norm_toa_irradiance())

# = formatted
    @store_value(__result_store)
    def solar_noon_time(self):
        return _vec_timedelta(days=self.solar_noon())

    @store_value(__result_store)
    def sunlight_time(self):
        return _vec_timedelta(days=self.sunlight())

    @store_value(__result_store)
    def sunrise_time(self):
        return _vec_timedelta(days=self.sunrise())

    @store_value(__result_store)
    def sunset_time(self):
        return _vec_timedelta(days=self.sunset())

# = handlers
    def summarize(self, verbose=False):
        compute_list = [f for f in dir(self)
                        if callable(getattr(self, f))
                        and not f.startswith("_")
                        and f != "summarize"]
        for f in compute_list:
            result = getattr(self, f).__call__()

            if verbose:
                if isinstance(result, np.ndarray):

                    # TODO: some variables are always cartesian vectors, having a ending dimension of 3
                    #   and it is inappropriate to average across the last dimension.
                    #   How to identify? better tracking of variable dimensionality? perhaps:
                    #   shape = (x, y, t, c) ->
                    #   space_dims = (True, True, False, False),
                    #   time_dims = (False, False, True, False)

                    # TODO hacky exception for cartesian vector variables
                    if "vec_cartesian" in f:
                        mean_axis = tuple(i for i in range(len(result.shape)-1))
                    else:
                        mean_axis = tuple(i for i in range(len(result.shape)))
                    result_mean = result.mean(axis=mean_axis)

                    # for vectorized types, verbose summary prints means across space and time dimensions
                    print(f"{f:27s} {str(result.shape):20s}{str(result_mean):27s}")

                else:
                    if isinstance(result, list):
                        dim = len(result)
                    else:
                        dim = 1,
                    print(f"{f:27s} {str(dim):20s}{str(result):27s}")


