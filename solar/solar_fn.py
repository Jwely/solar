from datetime import datetime, timezone
from typing import Tuple

import numpy as np
from numpy.core.umath import radians, sin, cos, degrees, arctan2, arcsin, tan, arccos
from solar.types import CONSTANTS, FlexNum, FlexDate

"""
Every function in this module accepts scalar or vectorized inputs.
These vectorized inputs can be used to represent a range of space and or times
simply, without complex shaping and reshaping.
"""


def get_absolute_julian_day_and_century(reference_datetime: FlexDate) -> Tuple[FlexNum, FlexNum]:
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


def get_geomean_long(ajc: FlexNum):
    """ calculates geometric mean longitude of the sun"""

    geomean_long = (280.46646 + ajc * (36000.76983 + ajc * 0.0003032)) % 360
    return geomean_long


def get_geomean_anom(ajc: FlexNum):
    """calculates the geometric mean anomoly of the sun"""

    geomean_anom = (357.52911 + ajc * (35999.05029 - 0.0001537 * ajc))
    return geomean_anom


def get_earth_eccent(ajc: FlexNum):
    """
    Calculates the precise eccentricity of earths orbit at ref_datetime
    """
    earth_eccent = 0.016708634 - ajc * (4.2037e-5 + 1.267e-7 * ajc)
    return earth_eccent


def get_sun_eq_of_center(ajc: FlexNum,
                         geomean_anom: FlexNum):
    """calculates the suns equation of center"""

    gma = radians(geomean_anom)

    sun_eq_of_center = \
        sin(gma) * (1.914602 - ajc * (0.004817 + 0.000014 * ajc)) + \
        sin(2 * gma) * (0.019993 - 0.000101 * ajc) + \
        sin(3 * gma) * 0.000289
    return sun_eq_of_center


def get_true_long(geomean_long: FlexNum,
                  sun_eq_of_center: FlexNum):
    """ calculates the tru longitude of the sun"""

    true_long = geomean_long + sun_eq_of_center
    return true_long


def get_true_anom(geomean_anom: FlexNum,
                  sun_eq_of_center: FlexNum):
    """ calculates the true anomoly of the sun"""
    true_anom = geomean_anom + sun_eq_of_center
    return true_anom


def get_rad_vector(earth_eccent: FlexNum,
                   true_anom: FlexNum):
    """ calculates incident radiation vector to surface at ref_datetime (AUs)"""

    ec = earth_eccent
    ta = radians(true_anom)

    rad_vector = (1.000001018 * (1 - ec ** 2)) / (1 + ec * cos(ta))
    return rad_vector


def get_app_long(true_long: FlexNum,
                 ajc: FlexNum):
    """ calculates apparent longitude of the sun"""
    app_long = true_long - 0.00569 - 0.00478 * sin(radians(125.04 - 1934.136 * ajc))
    return app_long


def get_oblique_mean_ellipse(ajc: FlexNum):
    """ calculates the oblique mean eliptic of earth orbit """

    oblique_mean_ellipse = 23 + (26 + (21.448 - ajc * (46.815 + ajc * (0.00059 - ajc * 0.001813))) / 60) / 60
    return oblique_mean_ellipse


def get_oblique_corr(ajc: FlexNum,
                     oblique_mean_ellipse: FlexNum):
    """ calculates the oblique correction """

    ome = oblique_mean_ellipse

    oblique_corr = ome + 0.00256 * cos(radians(125.04 - 1934.136 * ajc))
    return oblique_corr


def get_right_ascension(app_long: FlexNum,
                        oblique_corr: FlexNum):
    """ calculates the suns right ascension angle """

    sal = radians(app_long)
    oc = radians(oblique_corr)

    right_ascension = degrees(arctan2(cos(oc) * sin(sal), cos(sal)))
    return right_ascension


def get_declination(app_long: FlexNum,
                    oblique_corr: FlexNum):
    """ solar declination angle at ref_datetime"""

    sal = radians(app_long)
    oc = radians(oblique_corr)

    declination = degrees(arcsin((sin(oc) * sin(sal))))
    return declination


def get_equation_of_time(oblique_corr: FlexNum,
                         geomean_long: FlexNum,
                         geomean_anom: FlexNum,
                         earth_eccent: FlexNum):
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


def get_hour_angle_sunrise(declination: FlexNum,
                           lat_r: FlexNum):
    """ calculates the hour angle of sunrise """

    d = radians(declination)
    lat = lat_r     # radians

    hour_angle_sunrise = degrees(arccos((cos(radians(90.833)) / (cos(lat) * cos(d)) - tan(lat) * tan(d))))

    return hour_angle_sunrise


def get_solar_noon(lon: FlexNum,
                   equation_of_time: FlexNum,
                   tz=0):
    """ calculates solar noon in (local sidereal time LST)"""
    eot = equation_of_time

    solar_noon = (720 - 4 * lon - eot + tz * 60) / 1440
    return solar_noon


def get_sunrise(solar_noon: FlexNum,
                hour_angle_sunrise: FlexNum):
    """ calculates the time of sunrise"""

    sn = solar_noon
    ha = hour_angle_sunrise

    sunrise = (sn * 1440 - ha * 4) / 1440
    return sunrise


def get_sunset(solar_noon: FlexNum,
               hour_angle_sunrise: FlexNum):
    """ calculates the time of sunrise"""

    sn = solar_noon
    ha = hour_angle_sunrise

    sunset = (sn * 1440 + ha * 4) / 1440
    return sunset


def get_sunlight(hour_angle_sunrise: FlexNum):
    """ calculates amount of daily sunlight in fractional days"""

    sunlight = 8 * hour_angle_sunrise / (60 * 24)
    return sunlight


def get_true_solar(lon: FlexNum,
                   equation_of_time: FlexNum,
                   reference_datetime: FlexDate):
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

    frac_hr = frac_sec / (60 * 60)
    frac_day = frac_hr / 24

    frac_day = frac_day

    # now get true solar time
    true_solar = (frac_day * 1440 + eot + 4 * lon) % 1440

    return true_solar


def get_hour_angle(true_solar: FlexNum):
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


def get_zenith(declination: FlexNum,
               hour_angle: FlexNum,
               lat_r: FlexNum):
    """ calculates solar zenith angle at ref_datetime"""

    d = radians(declination)
    ha = radians(hour_angle)
    lat = lat_r

    zenith = degrees(arccos(sin(lat) * sin(d) + cos(lat) * cos(d) * cos(ha)))
    return zenith


# solar power related
def get_air_mass(zenith: FlexNum, method: str = None):
    """
    air mass is a concept from the solar power world that accounts for the mass of air
    between a point on the earths surface and the sun as a function of zenith angle.

    This method implements the Kasten-Young formula as well as a simple spherical approximation
    https://www.osapublishing.org/ao/abstract.cfm?uri=ao-28-22-4735

    :param method: Either "kasten-young" or "spherical".
    """

    if method is None:
        method = "kasten-young"  # default

    if method == "spherical":
        z_r = radians(zenith)
        r = CONSTANTS.earth_radius / CONSTANTS.atm_height
        return (r * (cos(z_r)**2) + (2*r) + 1) ** 0.5 - (r * cos(z_r))

    elif method == "kasten-young":
        z_r = radians(zenith)
        return 1.0 / (cos(z_r) + 0.50572 * (96.07995 - zenith)**(-1.6364))


# solar power related
def air_mass_correction_(
        _toa_irradiance: FlexNum,
        air_mass: FlexNum,
        pollution_condition: str = None):
    """
    Uses simple air_mass relationship to scale top of atmosphere irradiance to
    bottom of atmosphere irradiance.

    Wiki reference for pollution assumptions:
    https://en.wikipedia.org/wiki/Air_mass_(solar_energy)#cite_note-interpolation-17
    """

    if pollution_condition is None:
        pollution_condition = "average"

    if pollution_condition == "clean":
        base = 0.76
        exp = 0.618
    elif pollution_condition == "average":
        base = 0.7
        exp = 0.678
    elif pollution_condition == "dirty":
        base = 0.56
        exp = 0.715
    else:
        raise ValueError(f"Invalid 'pollution_condition' {pollution_condition},"
                         f" must be from [clean, dirty, average]")

    return 1.1 * _toa_irradiance * (base ** (air_mass ** exp))


def get_elevation_noatmo(zenith: FlexNum):
    """ solar elevation angle without atmospheric correction """
    return 90.0 - zenith


def get_elevation(zenith: FlexNum):
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
        lat_r: FlexNum,
        declination: FlexNum,
        hour_angle: FlexNum,
        zenith: FlexNum) -> FlexNum:
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

    In either case 2 or 3, all vector inputs must be identical shapes.
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

    # Case 4: Both space and time variation
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

    # Case 3: Space constant, time variation (lat_r in radians is scalar)
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

    # Case 2: Time constant, spatial variation (declination is scalar)
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

    # case 1, all inputs are scalar.
    else:

        if ha > 0:
            azimuth = (degrees(arccos(((sin(lat) * cos(z)) - sin(d)) / (cos(lat) * sin(z)))) + 180) % 360
        else:
            azimuth = (540 - degrees(arccos(((sin(lat) * cos(z)) - sin(d)) / (cos(lat) * sin(z))))) % 360

    return azimuth


def get_earth_distance(rad_vector: FlexNum):
    """ distance between the earth and the sun at ref_datetime"""

    # convert rad_vector length from AU to meters
    earth_distance = rad_vector * 149597870700

    return earth_distance


def get_sun_norm_toa_irradiance(earth_distance: FlexNum):
    """
    Calculates incoming solar energy in W/m^2 to a surface which is normal to the sun.

    instantaneous norm_irradiance is calculated as
        = sun_surf_radiance *(sun_radius / earth_distance)^2

    then corrected as a function of solar incidence angle
    """

    ed = earth_distance

    sun_surf_rad = CONSTANTS.sun_surf_rad
    sun_radius = CONSTANTS.sun_radius

    # calculate irradiance to normal surface at earth distance
    norm_irradiance = sun_surf_rad * (sun_radius / ed) ** 2

    return norm_irradiance


def _polar_to_cartesian(phi: FlexNum, theta: FlexNum):
    """
    Converts polar to cartesian coordinate unit vector (vector length of one).

    :param phi: phi is measured as an azimuth or aspect angle in the x,y plane. (degrees)
    :param theta: theta is measured as an angle from vertical, as zenith or slope angle. (degrees)

    note that phi may be 'Nan' where theta is zero (perfect z vector)
    """
    theta_r = radians(theta)
    phi_r = radians(phi)

    y = np.sin(theta_r) * np.cos(phi_r)
    x = np.sin(theta_r) * np.sin(phi_r)
    z = np.cos(theta_r)
    return np.stack([x, y, z], axis=-1)


def get_sun_vec_cartesian(azimuth: FlexNum, elevation: FlexNum):
    """
    Returns a vector pointing directly at the sun, from the surface on the earth, in cartesian coordinates.
        x axis = east positive, west negative
        y axis = north positive, south negative
        z axis = normal and up from earth spheroid positive.
    """
    # NOTE: elevation is used instead of zenith here to allow atmospheric corrections to impact sun vector.
    # polar coordinate to cartesian vector equations
    return _polar_to_cartesian(azimuth, 90.0 - elevation)


def get_surf_vec_cartesian(aspect: FlexNum, slope: FlexNum):
    """
    computes a surface normal vector in cartesian. Only used when slope and aspect are given.
        x axis = east positive, west negative
        y axis = north positive, south negative
        z axis = normal and up from earth spheroid positive.
    """
    if slope is None and aspect is None:
        return None

    else:
        # the angle off vertical of the surface normal vector is = to the slope
        return _polar_to_cartesian(aspect, slope)


def get_sun_surf_angle(sun_vec: FlexNum, surface_vec: FlexNum):
    """
    Gets the angle between the solar vector and the surface normal vector
    """

    # regardless of the simulation dimensionality, the inputs highest dimension should have
    # 3 elements for the X,Y,Z plane. So we can unstack all other dimensions.
    dims = sun_vec.shape
    out_dims = dims[:-1]
    flat_len = np.prod(out_dims)

    # if there are fewer than two dims, no need for einsum manipulations
    if len(sun_vec.shape) < 2:
        cosines = np.dot(sun_vec, surface_vec)
        angles = np.arccos(np.clip(cosines, -1.0, 1.0))
        return np.degrees(angles)

    else:
        if len(sun_vec.shape) == 2:
            v1 = sun_vec
            v2 = surface_vec

        else:   # len(sun_vec.shape) > 2
            v1 = sun_vec.reshape((flat_len, 3))
            v2 = surface_vec.reshape((flat_len, 3))

        cosines = np.einsum("ij,ij->i", v1.reshape(-1, 3), v2.reshape(-1, 3))
        angles = np.arccos(np.clip(cosines, -1.0, 1.0))
        return np.degrees(angles).reshape(out_dims)


def get_surf_norm_toa_irradiance(sun_surf_angle: FlexNum, sun_norm_irradiance: FlexNum):
    """
    gets the irradiance intensity normal to a surface.
    If slope and aspect inputs were used this accounts for the angle of the earths surface.
    """

    ssa_r = radians(sun_surf_angle)
    srf = cos(ssa_r) * sun_norm_irradiance  # a theta of zero results in full norm irradiance

    # set negative values to zero.
    if isinstance(srf, np.ndarray) and srf.shape:
        srf[srf <= 0] = 0
        return srf
    else:
        return max([srf, 0.0])


def get_surf_norm_toa_par(norm_surf_irradiance: FlexNum):
    """
    photosynthetically active radiation flux
    NOTE: using atmospheric surface irradiance creates discontinuities.
    """
    return norm_surf_irradiance * CONSTANTS.etta_par




