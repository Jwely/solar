"""
This module is for doing quick solar parameter calculations as they may impact
plant activity on the ground: it is a strictly computational approach based
on the geometry of the earth-sun system.

Other approaches include NASA data products regarding photosynthetically active radiation
such as MCD18A2 [1].


[1] https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mcd18a2_v006

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
    rad_vector              radiation vector (distance in AU)           (scalar)
    earth_distance          earths distance to sun in meters            (scalar)
    sun_norm_irradiance     incident solar energy at earth distance     (scalar)
    surf_norm_irradiance    incident energy to surface

Units unless otherwise labeled.
    angle    = degrees
    distance = meters
    energy   = watts or joules
    time     = mostly in datetime objects or absolute julian time. ALL TIME is in UTC time zone!

"""

from jwelysolar.calculator import SolarCalculator

