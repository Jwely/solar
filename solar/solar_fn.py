from datetime import datetime
from typing import Tuple

import numpy as np
from numpy.core.umath import radians, sin, cos, degrees, arctan2, arcsin, tan, arccos

from solar import Numeric, DateTime
from solar import CONSTANTS

"""
Every function in this module accepts scalar or vectorized inputs.
These vectorized inputs can be used to represent a range of space OR times
simply, without complex shaping and reshaping.
"""


def get_absolute_julian_day_and_century(reference_datetime: DateTime) -> Tuple[Numeric, Numeric]:
    """ computes absolute julian day and century from some UTC datetime """

    # self reference if reference_datetime is a list
    if isinstance(reference_datetime, list):
        ajd_ajc_tuples = [get_absolute_julian_day_and_century(rdt) for rdt in reference_datetime]
        _ajd, _ajc = zip(*ajd_ajc_tuples)
        ajd = np.array(list(_ajd), dtype=float)
        ajc = np.array(list(_ajc), dtype=float)
        return ajd, ajc

    else:

        rdt = reference_datetime

        # uses the reference day of january 1st 2000
        jan_1st_2000_jd = 2451545
        jan_1st_2000 = datetime(2000, 1, 1, 12, 0, 0)

        time_del = rdt - jan_1st_2000
        ajd = float(jan_1st_2000_jd) + float(time_del.total_seconds()) / 86400
        ajc = (ajd - 2451545) / 36525.0
    return ajd, ajc


def get_geomean_long(ajc: Numeric):
    """ calculates geometric mean longitude of the sun"""

    geomean_long = (280.46646 + ajc * (36000.76983 + ajc * 0.0003032)) % 360
    return geomean_long


def get_geomean_anom(ajc: Numeric):
    """calculates the geometric mean anomoly of the sun"""

    geomean_anom = (357.52911 + ajc * (35999.05029 - 0.0001537 * ajc))
    return geomean_anom


def get_earth_eccent(ajc: Numeric):
    """
    Calculates the precise eccentricity of earths orbit at ref_datetime
    """
    earth_eccent = 0.016708634 - ajc * (4.2037e-5 + 1.267e-7 * ajc)
    return earth_eccent


def get_sun_eq_of_center(ajc: Numeric,
                         geomean_anom: Numeric):
    """calculates the suns equation of center"""

    gma = radians(geomean_anom)

    sun_eq_of_center = \
        sin(gma) * (1.914602 - ajc * (0.004817 + 0.000014 * ajc)) + \
        sin(2 * gma) * (0.019993 - 0.000101 * ajc) + \
        sin(3 * gma) * 0.000289
    return sun_eq_of_center


def get_true_long(geomean_long: Numeric,
                  sun_eq_of_center: Numeric):
    """ calculates the tru longitude of the sun"""

    true_long = geomean_long + sun_eq_of_center
    return true_long


def get_true_anom(geomean_anom: Numeric,
                  sun_eq_of_center: Numeric):
    """ calculates the true anomoly of the sun"""
    true_anom = geomean_anom + sun_eq_of_center
    return true_anom


def get_rad_vector(earth_eccent: Numeric,
                   true_anom: Numeric):
    """ calculates incident radiation vector to surface at ref_datetime (AUs)"""

    ec = earth_eccent
    ta = radians(true_anom)

    rad_vector = (1.000001018 * (1 - ec ** 2)) / (1 + ec * cos(ta))
    return rad_vector


def get_app_long(true_long: Numeric,
                 ajc: Numeric):
    """ calculates apparent longitude of the sun"""
    app_long = true_long - 0.00569 - 0.00478 * sin(radians(125.04 - 1934.136 * ajc))
    return app_long


def get_oblique_mean_ellipse(ajc: Numeric):
    """ calculates the oblique mean eliptic of earth orbit """

    oblique_mean_ellipse = 23 + (26 + (21.448 - ajc * (46.815 + ajc * (0.00059 - ajc * 0.001813))) / 60) / 60
    return oblique_mean_ellipse


def get_oblique_corr(ajc: Numeric,
                     oblique_mean_ellipse: Numeric):
    """ calculates the oblique correction """

    ome = oblique_mean_ellipse

    oblique_corr = ome + 0.00256 * cos(radians(125.04 - 1934.136 * ajc))
    return oblique_corr


def get_right_ascension(app_long: Numeric,
                        oblique_corr: Numeric):
    """ calculates the suns right ascension angle """

    sal = radians(app_long)
    oc = radians(oblique_corr)

    right_ascension = degrees(arctan2(cos(oc) * sin(sal), cos(sal)))
    return right_ascension


def get_declination(app_long: Numeric,
                    oblique_corr: Numeric):
    """ solar declination angle at ref_datetime"""

    sal = radians(app_long)
    oc = radians(oblique_corr)

    declination = degrees(arcsin((sin(oc) * sin(sal))))
    return declination


def get_equation_of_time(oblique_corr: Numeric,
                         geomean_long: Numeric,
                         geomean_anom: Numeric,
                         earth_eccent: Numeric):
    """ calculates the equation of time in minutes """

    oc = radians(oblique_corr)
    gml = radians(geomean_long)
    gma = radians(geomean_anom)
    ec = earth_eccent

    vary = tan(oc / 2) ** 2

    equation_of_time = 4 * degrees(
        vary * sin(2 * gml) - 2 * ec * sin(gma) + \
        4 * ec * vary * sin(gma) * cos(2 * gml) - \
        0.5 * vary * vary * sin(4 * gml) - \
        1.25 * ec * ec * sin(2 * gma))

    return equation_of_time


def get_hour_angle_sunrise(declination: Numeric,
                           lat_r: Numeric):
    """ calculates the hour angle of sunrise """

    d = radians(declination)
    lat = lat_r     # radians

    hour_angle_sunrise = degrees(arccos((cos(radians(90.833)) / (cos(lat) * cos(d)) - tan(lat) * tan(d))))

    return hour_angle_sunrise


def get_solar_noon(lon: Numeric,
                   equation_of_time: Numeric,
                   tz=0):
    """ calculates solar noon in (local sidereal time LST)"""

    eot = equation_of_time

    solar_noon = (720 - 4 * lon - eot + tz * 60) / 1440
    return solar_noon


def get_sunrise(solar_noon: Numeric,
                hour_angle_sunrise: Numeric):
    """ calculates the time of sunrise"""

    sn = solar_noon
    ha = hour_angle_sunrise

    sunrise = (sn * 1440 - ha * 4) / 1440
    return sunrise


def get_sunset(solar_noon: Numeric,
               hour_angle_sunrise: Numeric):
    """ calculates the time of sunrise"""

    sn = solar_noon
    ha = hour_angle_sunrise

    sunset = (sn * 1440 + ha * 4) / 1440
    return sunset


def get_sunlight(hour_angle_sunrise: Numeric):
    """ calculates amount of daily sunlight in fractional days"""

    sunlight = 8 * hour_angle_sunrise / (60 * 24)
    return sunlight


def get_true_solar(lon: Numeric,
                   equation_of_time: Numeric,
                   reference_datetime: DateTime,
                   tz=0):
    """ calculates the true solar time at ref_datetime"""

    eot = equation_of_time
    rdt = reference_datetime

    # turn reference datetime into fractional days
    if isinstance(rdt, datetime):
        frac_sec = (rdt - datetime(rdt.year, rdt.month, rdt.day)).total_seconds()
    elif isinstance(rdt, list):
        frac_sec_list = [(_rdt - datetime(_rdt.year, _rdt.month, _rdt.day)).total_seconds() for _rdt in rdt]
        frac_sec = np.array(frac_sec_list)
    else:
        raise TypeError("expected reference_datetime to be datetime or "
                        "list of datetimes. got {}".format(type(reference_datetime)))

    frac_hr = frac_sec / (60 * 60) + tz
    frac_day = frac_hr / 24

    frac_day = frac_day

    # now get true solar time
    true_solar = (frac_day * 1440 + eot + 4 * lon - 60 * tz) % 1440

    return true_solar


def get_hour_angle(true_solar: Numeric):
    """ calculates the hour angle at ref_datetime"""

    ts = true_solar

    # matrix hour_angle calculations
    if isinstance(true_solar, np.ndarray):
        ha = ts
        ha[ha <= 0] = ha[ha <= 0] / 4 + 180
        ha[ha > 0] = ha[ha > 0] / 4 - 180
        hour_angle = ha

    # scalar hour_angle calculations
    else:
        if ts <= 0:
            hour_angle = ts / 4 + 180
        else:
            hour_angle = ts / 4 - 180

    return hour_angle


def get_zenith(declination: Numeric,
               hour_angle: Numeric,
               lat_r: Numeric):
    """ calculates solar zenith angle at ref_datetime"""

    d = radians(declination)
    ha = radians(hour_angle)
    lat = lat_r

    zenith = degrees(arccos(sin(lat) * sin(d) + cos(lat) * cos(d) * cos(ha)))
    return zenith


def get_elevation(zenith: Numeric):
    """ calculates solar elevation angle at ref_datetime"""

    # perform an approximate atmospheric refraction correction

    # matrix hour_angle calculations
    # these equations are hideous, but im not sure how to improve them without adding computational complexity
    if isinstance(zenith, np.ndarray) and zenith.shape:
        e = 90.0 - zenith
        ar = e * 0

        ar[e > 85] = 0
        ar[(e > 5) & (e <= 85)] = 58.1 / tan(radians(e[(e > 5) & (e <= 85)])) - \
                                  0.07 / tan(radians(e[(e > 5) & (e <= 85)])) ** 3 + \
                                  0.000086 / tan(radians(e[(e > 5) & (e <= 85)])) ** 5
        ar[(e > -0.575) & (e <= 5)] = 1735 + e[(e > -0.575) & (e <= 5)] * \
                                             (103.4 + e[(e > -0.575) & (e <= 5)] * (
                                                 -12.79 + e[(e > -0.575) & (e <= 5)] * 0.711))
        ar[e <= -0.575] = -20.772 / tan(radians(e[e <= -0.575]))

    # scalar hour_angle calculations
    else:
        e = 90.0 - zenith
        er = radians(e)

        if e > 85:
            ar = 0
        elif e > 5:
            ar = 58.1 / tan(er) - 0.07 / tan(er) ** 3 + 0.000086 / tan(er) ** 5
        elif e > -0.575:
            ar = 1735 + e * (103.4 + e * (-12.79 + e * 0.711))
        else:
            ar = -20.772 / tan(er)

    elevation_noatmo = e
    atmo_refraction = ar / 3600
    elevation = elevation_noatmo + atmo_refraction

    return elevation


def get_azimuth(
        lat_r: Numeric,
        declination: Numeric,
        hour_angle: Numeric,
        zenith: Numeric) -> Numeric:
    """
    calculates solar azimuth angle. Function requires special treatment
    of different cases of combine vectorize inputs

    Case 1: all scalar inputs
        all scalar inputs are simple, and a scalar will be output.

    Case 2: Time constant, spatial variation (declination is scalar)
        lat_r:          vector (x,y)
        declination:    scalar (z)
        hour_angle:     vector (x,y)
        zenith:         vector (x,y)
        azimuth output  vector (x,y)

    Case 3: Space constant, time variation (lat_r in radians is scalar)
        lat_r:          scalar (x)
        declination:    vector (y)
        hour_angle:     vector (y)
        zenith:         vector (y)
        azimuth output  vector (y)

    case 4: Both space and time variation
        lat_r:          vector (x,y,z)  # duplicated across z axis
        declination:    vector (z)      # will be duplicated to dims (x,y,z)
        hour_angle:     vector (x,y,z)  # duplicated across z axis
        zenith:         vector (x,y,z)  # duplicated across z axis
        azimuth output  vector (x,y,z)  # duplicated across z axis

    In either case 2 or three, all vector inputs must be identical shapes.
    """

    lat = lat_r
    d = radians(declination)    # always returns numpy arrays, even if 1x1
    ha = radians(hour_angle)    # always returns numpy arrays, even if 1x1
    z = radians(zenith)         # always returns numpy arrays, even if 1x1

    # are the inputs vectorized in space? (lat_r will be a vector)
    if isinstance(lat_r, np.ndarray):
        if lat_r.shape:
            lat_vec = True
        else:
            lat_vec = False     # ndarray with no dimensions is a scalar!
    else:
        lat_vec = False

    if isinstance(d, np.ndarray):
        if d.shape:
            time_vec = True
        else:
            time_vec = False    # ndarray with no dimensions is a scalar!
    else:
        time_vec = False

    if time_vec and lat_vec:
        # if both time and space are vectorized, we must cast declination(time)
        # to the third dimension. by duplicating for all lat/lon pairs.

        if len(lat_r.shape) == 3 and len(d.shape) == 1:
            if lat_r.shape[2] == d.shape[0]:
                d = np.repeat(d[np.newaxis, :], lat_r.shape[1], axis=0)
                d = np.repeat(d[np.newaxis, :], lat_r.shape[0], axis=0)

        az = ha * 0

        ha_p = (ha > 0)     # positive ha indices
        ha_n = (ha <= 0)    # negative ha indices

        az_ha_p = arccos(
            ((sin(lat[ha_p]) * cos(z[ha_p])) - sin(d[ha_p]))
            / (cos(lat[ha_p]) * sin(z[ha_p])))
        az[ha_p] = (degrees(az_ha_p) + 180) % 360

        az_ha_n = arccos(
            ((sin(lat[ha_n]) * cos(z[ha_n])) - sin(d[ha_n]))
            / (cos(lat[ha_n]) * sin(z[ha_n]))
        )
        az[ha_n] = (540 - degrees(az_ha_n)) % 360
        azimuth = az
        # raise NotImplementedError(
        #     " This function does not support simultaneous time and spatial vectorization."
        #     " Only lat, OR declination arguments may be vectors. If either are vectors, then"
        #     " so must be hour_angle and zenith.")

    elif lat_vec and not time_vec:
        az = ha * 0

        ha_p = (ha > 0)     # positive ha indices
        ha_n = (ha <= 0)    # negative ha indices

        az_ha_p = arccos(
            ((sin(lat[ha_p]) * cos(z[ha_p])) - sin(d))
            / (cos(lat[ha_p]) * sin(z[ha_p])))
        az[ha_p] = (degrees(az_ha_p) + 180) % 360

        az_ha_n = arccos(
            ((sin(lat[ha_n]) * cos(z[ha_n])) - sin(d))
            / (cos(lat[ha_n]) * sin(z[ha_n]))
        )
        az[ha_n] = (540 - degrees(az_ha_n)) % 360
        azimuth = az

    elif time_vec and not lat_vec:
        az = ha * 0

        ha_p = (ha > 0)     # positive ha indices
        ha_n = (ha <= 0)    # negative ha indices

        az_ha_p = arccos(
            ((sin(lat) * cos(z[ha_p])) - sin(d[ha_p]))
            / (cos(lat) * sin(z[ha_p])))
        az[ha_p] = (degrees(az_ha_p) + 180) % 360

        az_ha_n = arccos(
            ((sin(lat) * cos(z[ha_n])) - sin(d[ha_n]))
            / (cos(lat) * sin(z[ha_n]))
        )
        az[ha_n] = (540 - degrees(az_ha_n)) % 360
        azimuth = az

    else:   # all inputs are scalar

        if ha > 0:
            azimuth = (degrees(arccos(((sin(lat) * cos(z)) - sin(d)) / (cos(lat) * sin(z)))) + 180) % 360
        else:
            azimuth = (540 - degrees(arccos(((sin(lat) * cos(z)) - sin(d)) / (cos(lat) * sin(z))))) % 360

    return azimuth


def get_earth_distance(rad_vector: Numeric):
    """ distance between the earth and the sun at ref_datetime"""

    # convert rad_vector length from AU to meters
    earth_distance = rad_vector * 149597870700

    return earth_distance


def get_norm_irradiance(earth_distance: Numeric):
    """
    calculates incoming solar energy in W/m^2 to a surface normal to the sun

    inst_irradiance is calculated as
        = sun_surf_radiance *(sun_radius / earth_distance)^2

    then corrected as a function of solar incidence angle
    """

    ed = earth_distance

    sun_surf_rad = CONSTANTS.sun_surf_rad
    sun_radius = CONSTANTS.sun_radius

    # calculate irradiance to normal surface at earth distance
    norm_irradiance = sun_surf_rad * (sun_radius / ed) ** 2

    return norm_irradiance

