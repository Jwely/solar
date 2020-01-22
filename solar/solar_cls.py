import numpy as np
from datetime import timedelta

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


class SolarCalculator(object):
    """
    Chains solar equations as needed, saving intermediates as attributes of this class.

    Table of variable descriptions

        lat                     decimal degrees latitude (float OR numpy array)
        lon                     decimal degrees longitude (float OR numpy array)
        reference_datetime      datetime at UTC
        slope                   slope of land at lat,lon for solar energy calculations
        aspect                  aspect of land at lat,lon for solar energy calculations
        -------------------------------------------------------------------------------
        lat                     latitude                                    (array)
        lon                     longitude                                   (array)
        tz                      time zone                                   (scalar)
        rdt                     reference datetime object (date_time_obj)   (scalar)
        ajd                     absolute julian day                         (scalar)
        ajc                     absolute julian century                     (scalar)
        geomean_long            geometric mean longitude of the sun         (scalar)
        geomean_anom            geometric mean longitude anomaly of the sun (scalar)
        earth_eccent            eccentricity of earths orbit                (scalar)
        sun_eq_of_center        the suns equation of center                 (scalar)
        true_long               true longitude of the sun                   (scalar)
        true_anom               true longitude anomaly of the sun           (scalar)
        app_long                the suns apparent longitude                 (scalar)
        oblique_mean_ellipse    earth oblique mean ellipse                  (scalar)
        oblique_corr            correction to earths oblique ellipse        (scalar)
        right_ascension         suns right ascension angle                  (scalar)
        declination             solar declination angle                     (scalar)
        equation_of_time        equation of time (minutes)                  (scalar)
        hour_angle_sunrise      the hour angle at sunrise                   (array)
        solar_noon              LST of solar noon                           (array)
        sunrise                 LST of sunrise time                         (array)
        sunset                  LST of sunset time                          (array)
        sunlight                LST fractional days of sunlight             (array)
        true_solar              LST for true solar time                     (array)
        hour_angle              total hour angle                            (array)
        zenith                  zenith angle                                (array)
        elevation               elevation angle                             (array)
        azimuth                 azimuthal angle                             (array)
        norm_surf_irradiance    incident energy on earths surface           (array)
        surface_par             photosynthetically active radiation         (array)
        rad_vector              radiation vector (distance in AU)           (scalar)
        earth_distance          earths distance to sun in meters            (scalar)
        norm_irradiance         incident solar energy at earth distance     (scalar)

    Units unless otherwise labeled.
        angle    = degrees
        distance = meters
        energy   = watts or joules
        time     = mostly in datetime objects or absolute julian time. ALL TIME is in UTC time zone!
    """

    # TODO: More precise shape descriptions, space and time variant support has
    #     # TODO: made the expected behavior slightly more complex.
    def __init__(self, lat: FlexNum, lon: FlexNum, dt: FlexDate,
                 slope: FlexNum = None, aspect: FlexNum = None, low_mem=False):
        """
        Accepts vectorized positions for lat and lon, (use meshgrid) and/or list of
        multiple datetimes. Computations are fully vectorized

        :param lat:     decimal degrees latitude
        :param lon:     decimal degrees longitude
        :param dt:  reference datetime or a list of datetimes
        :param low_mem: if False, intermediate calculations will not be preserved. Repeated calls
                        will recompute all values including intermediates every time.
                        Use for very large vector spaces.
        """

        if slope is not None:
            raise NotImplementedError(f"'slope' argument not yet supported")

        if aspect is not None:
            raise NotImplementedError(f"'aspect' argument not yet supported")

        self.low_mem = low_mem  # if True, do not save intermediate computations

        # check inputs for dimensional compatibility
        lat, lon, dt, vt, vs = self._coerce_dims(lat, lon, dt)

        self._time_is_vec = vt
        self._space_is_vec = vs

        self._dt = dt
        self._lat = lat
        self._lon = lon
        self._lat_rad = np.radians(lat)
        self._lon_rad = np.radians(lon)

        self._ajd, self._ajc = sfn.get_absolute_julian_day_and_century(dt)

# = input checks and properties of the class
    @staticmethod
    def _coerce_dims(lat: FlexNum, lon: FlexNum, dt: FlexDate):
        """
        Ensures that input lat, lon, and dt arguments are dimensionally compatible.

        latitude and longitude should be scalars or form a meshgrid.

        :param lat:
        :param lon:
        :param dt:
        :return:
        """

        # does this solver vectorize calculations over time?
        if isinstance(dt, list):
            vec_time = True
        else:
            vec_time = False

        if isinstance(lat, np.ndarray) and isinstance(lon, np.ndarray):

            vec_space = True

            # if lat and lon shapes are not equal
            if lat.shape != lon.shape or lat.shape == 1:
                lat, lon = np.meshgrid(lat, lon)

            # if both space and time are vectorized, add 3rd time dimension to lat/lon
            space_shape = lat.shape
            if vec_time:
                if len(space_shape) == 2:
                    # make lat and lon 3 dimensions for compatibility with dt
                    lat = np.repeat(lat[:, :, np.newaxis], len(dt), axis=2)
                    lon = np.repeat(lon[:, :, np.newaxis], len(dt), axis=2)
        else:
            vec_space = False

        return lat, lon, dt, vec_time, vec_space

    @property
    def time_is_vec(self):
        return self._time_is_vec

    @property
    def space_is_vec(self):
        return self._space_is_vec

# = function decorators and meta handlers
    def _ref(self, function, attribute, *args, **kwargs):
        """
        decorator to check the given attribute for the result of a function first,
        if it"s None, then execute the function with the input arguments and store
        the result in the attribute (if store=True).

        for example: the function `geomean_long` is writen as



        ```
            def geomean_long(self):
                return self.reference_hidden(get_geomean_long, "_geomean_long")(self._ajc)
        ```
        which is equivalent to

        ```
            def geomean_long(self):
                result = self._geomean_long
                if result is None:
                    result = get_geomean_long(self._ajc)
                    self._geomean_long = result
                return result
        ```
        """

        def wrapped():
            if hasattr(self, attribute):
                result = getattr(self, attribute)
            else:
                result = None

            if result is None:
                result = function(*args, **kwargs)

                if self.low_mem:
                    setattr(self, attribute, result)
            return result

        return wrapped()

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

    def ajc(self):
        return self._ajc

    def ajd(self):
        return self._ajd

    def dt(self):
        return self._dt

# = computed on the fly
    def geomean_long(self):
        return self._ref(sfn.get_geomean_long, "_geomean_long",
                         self.ajc())

    def geomean_anom(self):
        return self._ref(sfn.get_geomean_anom, "_geomean_anom",
                         self.ajc())

    def earth_eccent(self):
        return self._ref(sfn.get_earth_eccent, "_earth_eccent",
                         self.ajc())

    def sun_eq_of_center(self):
        return self._ref(sfn.get_sun_eq_of_center, "_sun_eq_of_center",
                         self.ajc(), self.geomean_anom())

    def true_long(self):
        return self._ref(sfn.get_true_long, "_true_long",
                         self.geomean_long(), self.sun_eq_of_center())

    def true_anom(self):
        return self._ref(sfn.get_true_anom, "_true_anom",
                         self.geomean_anom(), self.sun_eq_of_center())

    def rad_vector(self):
        return self._ref(sfn.get_rad_vector, "_rad_vector",
                         self.earth_eccent(), self.true_anom())

    def app_long(self):
        return self._ref(sfn.get_app_long, "_app_long",
                         self.true_long(), self.ajc())

    def oblique_mean_ellipse(self):
        return self._ref(sfn.get_oblique_mean_ellipse, "_oblique_mean_ellipse",
                         self.ajc())

    def oblique_corr(self):
        return self._ref(sfn.get_oblique_corr, "_oblique_corr",
                         self.ajc(), self.oblique_mean_ellipse())

    def right_ascension(self):
        return self._ref(sfn.get_right_ascension, "_right_ascension",
                         self.app_long(), self.oblique_corr())

    def declination(self):
        return self._ref(sfn.get_declination, "_declination",
                         self.app_long(), self.oblique_corr())

    def equation_of_time(self):
        return self._ref(sfn.get_equation_of_time, "_equation_of_time",
                         self.oblique_corr(), self.geomean_long(), self.geomean_anom(), self.earth_eccent())

    def hour_angle_sunrise(self):
        return self._ref(sfn.get_hour_angle_sunrise, "_hour_angle_sunrise",
                         self.declination(), self.lat_rad())

    def solar_noon(self):
        return self._ref(sfn.get_solar_noon, "_solar_noon",
                         self.lon(), self.equation_of_time())

    def sunrise(self):
        return self._ref(sfn.get_sunrise, "_sunrise",
                         self.solar_noon(), self.hour_angle_sunrise())

    def sunset(self):
        return self._ref(sfn.get_sunset, "_sunset",
                         self.solar_noon(), self.hour_angle_sunrise())

    def sunlight(self):
        return self._ref(sfn.get_sunlight, "_sunlight",
                         self.hour_angle_sunrise())

    def true_solar(self):
        return self._ref(sfn.get_true_solar, "_true_solar",
                         self.lon(), self.equation_of_time(), self.dt())

    def hour_angle(self):
        return self._ref(sfn.get_hour_angle, "_hour_angle",
                         self.true_solar())

    def zenith(self):
        return self._ref(sfn.get_zenith, "_zenith",
                         self.declination(), self.hour_angle(), self.lat_rad())

    def elevation(self):
        return self._ref(sfn.get_elevation, "_elevation",
                         self.zenith())

    def elevation_noatmo(self):
        return self._ref(sfn.get_elevation_noatmo, "_elevation_noatmo",
                         self.zenith())

    def azimuth(self):
        return self._ref(sfn.get_azimuth, "_azimuth",
                         self.lat_rad(), self.declination(), self.hour_angle(), self.zenith())

    def earth_distance(self):
        return self._ref(sfn.get_earth_distance, "_earth_distance",
                         self.rad_vector())

    def norm_irradiance(self):
        return self._ref(sfn.get_norm_irradiance, "_norm_irradiance",
                         self.earth_distance())

    def norm_surf_irradiance_noatmo(self):
        return self._ref(sfn.get_norm_surf_irradiance, "_norm_surf_irradiance_noatmo",
                         self.elevation_noatmo(), self.norm_irradiance())

    def norm_surf_irradiance(self):
        return self._ref(sfn.get_norm_surf_irradiance, "_norm_surf_irradiance",
                         self.elevation(), self.norm_irradiance())

    def surface_par(self):
        return self._ref(sfn.get_surface_par, "_surface_par",
                         self.norm_surf_irradiance_noatmo())

# = formatted
    def solar_noon_time(self):
        return self._ref(_vec_timedelta, '_solar_noon_time',
                         days=self.solar_noon())

    def sunlight_time(self):
        return self._ref(_vec_timedelta, '_sunlight_time',
                         days=self.sunlight())

    def sunrise_time(self):
        return self._ref(_vec_timedelta, '_sunrise_time',
                         days=self.sunrise())

    def sunset_time(self):
        return self._ref(_vec_timedelta, '_sunset_time',
                         days=self.sunset())

# = compute all
    # TODO: improve summary function
    def compute_all(self, verbose=False):
        """

        :param verbose:
        :return:
        """
        compute_list = [f for f in dir(self)
                        if callable(getattr(self, f))
                        and not f.startswith("_")
                        and not f.startswith("__")
                        and f != "compute_all"]
        for f in compute_list:
            result = getattr(self, f).__call__()

            if verbose:
                if isinstance(result, np.ndarray):
                    # for vectorized types, verbose summary prints means across all vector dimensions
                    print("{:26s}{:26s}\taveraged over matrix shape={}".format(
                        f, str(result.mean()), result.shape))
                # elif isinstance(result, list):
                #     print("{:26s}{:26s}\taveraged over matrix shape={}".format(
                #         f, str(np.median(result)), len(result)))
                else:
                    print("{:26s}{}".format(f, result))






















































