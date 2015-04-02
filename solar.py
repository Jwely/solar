"""
solar object class, built in my spare time for a personal project, but
will eventually be merged into the dnppy module
"""


__author__ = "Jeffry Ely, Jeff.ely.08@gmail.com"


from datetime import datetime, timedelta
from numpy import *

class solar:
    """
    Object class for handling solar calculations. Many equations are taken from the
    excel sheet at this url : [http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html]

    It requires a physical location on the earth and a datetime object

    Inputs:
        lat             lattitude (+ to N)
        lon             longitude (+ to E)
        time_zone       integer such as (-5) for eastern time (+ to E)
        date_time_obj   local time datetime object, may be datestamp with matching fmt
        fmt             format for interpreting "date_time_obj" if it is of string type

    Units used by this class unless otherwise labled:
        angle    = degrees
        distance = meters
        energy   = watts or joules
        time     = mostly in datetime objects.

    PLANNED IMPROVEMENT:
    1) inputs of numpy arrays for lat and lon needs to be allowed.
    2) inputs of a numpy array DEM for slope/aspect effects on incident solar energy
    intensity needto be allowed as well.
    
    """

    def __init__(self, lat, lon, time_zone, date_time_obj, fmt = False, metadata = None):
        """
        initializes critical spatial and temporal information

        lat             decimal degrees latitude
        lon             decimal degrees longitude
        time_zone       float of time shift from GMT (such as "-5" for EST)
        date_time_obj   either a timestamp string following fmt or a datetime obj
        fmt             if date_time_obj is a string, fmt is required to interpret it
        """

        # attribute list
        self.rdt            = []                # ref_datetime. Precise GMT time of solar obj
        self.tz             = 0.0               # timezone, offset from GMT - float
        self.ajd            = 0.0               # absolute julian day
        self.ajc            = 0.0               # absolute julian century
        self.earth_eccent   = 0.0               # eccentricity of earths orbit


        # Constants as attributes
        self.sun_surf_rad   = 63156942.6        # radiation at suns surface (W/m^2)
        self.sun_radius     = 695800000.        # radius of the sun in meters
        self.orbital_period = 365.2563630       # num of days it takes earth to revolve
        self.altitude       = -0.01448623       # altitude of center of solar disk


        # sets up the object with some subfunctions
        self.set_datetime(date_time_obj, fmt)
        self.set_location(lat, lon)
        self.set_timezone(time_zone)


        # compute solar attributes
        if isinstance(lat, ndarray) and isinstance(lon, ndarray):
            self.is_numpy   = True
            self.metadata   = metadata
            self.compute_all_numpy()
            
        else:
            self.is_numpy   = False
            self.compute_all_scalar()
            


    def set_datetime(self, date_time_obj, fmt = False):
        """
        sets the critical time information

        accepts datetime objects or a datetime string with format
        """

        # if input is datetime_obj set it
        if isinstance(date_time_obj, datetime):
            self.rdt = date_time_obj

        elif isinstance(date_time_obj, str) and isinstance(fmt, str):
            self.rdt = datetime.strptime(date_time_obj,fmt)
        else:
            raise Exception("bad datetime!")

        self.abs_julian()
        
        return


    def abs_julian(self):
        """
        calculates the absolute decimal julian day, century of ref_datetime

        This is invoked whenever set_datetime or set_timezone is invoked
        """

        # uses the reference day of january 1st 2000
        jan_1st_2000_jd   = 2451545
        jan_1st_2000      = datetime(2000,1,1,12,0,0)

        time_del = self.rdt - jan_1st_2000
        self.ajd = float(jan_1st_2000_jd) + float(time_del.total_seconds())/86400
        self.ajc = (self.ajd - 2451545)/36525.0

        return

    
    def set_location(self, lat, lon):
        """ sets spatial information"""

        self.lat    = lat           # lattitude (E positive)- float
        self.lon    = lon           # longitude (N positive)- float
        self.lat_r  = radians(lat)  # lattitude in radians
        self.lon_r  = radians(lon)  # longitude in radians
        return


    def set_timezone(self, GMT_hour_offset):
        """ amends the datetime object with time zone information"""

        # sets time zone in a way to prevent accidental repeated offsetting
        if self.tz != GMT_hour_offset:
            time_del = timedelta(hours = (self.tz - GMT_hour_offset))
            self.rdt = self.rdt + time_del
            self.tz  = GMT_hour_offset

        self.abs_julian()
        return


    def get_geomean_long(self):
        """ calculates geometric mean longitude of the sun"""

        self.geomean_long = (280.46646 + self.ajc * (36000.76983 + self.ajc*0.0003032)) % 360
        return self.geomean_long


    def get_geomean_anom(self):
        """calculates the geometric mean anomoly of the sun"""

        self.geomean_anom = (357.52911 + self.ajc * (35999.05029 - 0.0001537 * self.ajc))
        return self.geomean_anom

    
    def get_earth_eccent(self):
        """
        Calculates the precise eccentricity of earths orbit at ref_datetime
        """
        
        self.earth_eccent = 0.016708634 - self.ajc * (4.2037e-5 + 1.267e-7 * self.ajc)
        return self.earth_eccent


    def get_sun_eq_of_center(self):
        """calculates the suns equation of center"""

        ajc = self.ajc
        gma = radians(self.geomean_anom)

        self.sun_eq_of_center = sin(gma) * (1.914602 - ajc*(0.004817 + 0.000014 * ajc)) + \
                                sin(2*gma) * (0.019993 - 0.000101 * ajc) + \
                                sin(3*gma) * 0.000289

        return self.sun_eq_of_center


    def get_true_long(self):
        """ calculates the tru longitude of the sun"""

        self.true_long = self.geomean_long + self.sun_eq_of_center
        return self.true_long


    def get_true_anom(self):
        """ calculates the true anomoly of the sun"""

        self.true_anom = self.geomean_anom + self.sun_eq_of_center
        return self.true_anom

        
    def get_rad_vector(self):
        """ calculates incident radiation vector to surface at ref_datetime (AUs)"""

        ec = self.earth_eccent
        ta = radians(self.true_anom)
        
        
        self.rad_vector  = (1.000001018*(1 - ec**2)) / (1 + ec *cos(ta))
        return self.rad_vector


    def get_app_long(self):
        """ calculates apparent longitude of the sun"""

        stl = self.true_long
        ajc = self.ajc

        self.app_long = stl - 0.00569 - 0.00478 * sin(radians(125.04 - 1934.136 * ajc))
        return self.app_long
        

    def get_oblique_mean_elip(self):
        """ calculates the oblique mean eliptic of earth orbit """

        ajc = self.ajc

        self.oblique_mean_elip = 23 + (26 + ((21.448 - ajc * (46.815 + ajc * (0.00059 - ajc * 0.001813))))/60)/60
        return self.oblique_mean_elip


    def get_oblique_corr(self):
        """ calculates the oblique correction """

        ome = self.oblique_mean_elip
        ajc = self.ajc
        
        self.oblique_corr = ome + 0.00256 * cos(radians(125.04 - 1934.136 * ajc))
        return self.oblique_corr


    def get_right_ascension(self):
        """ calculates the suns right ascension angle """

        sal = radians(self.app_long)
        oc  = radians(self.oblique_corr)

        self.right_ascension = degrees(arctan2(cos(oc) * sin(sal), cos(sal)))
        return self.right_ascension
        
    
    def get_declination(self):
        """ solar declination angle at ref_datetime"""

        sal = radians(self.app_long)
        oc  = radians(self.oblique_corr)
        
        self.declination = degrees(arcsin((sin(oc) * sin(sal))))
        return self.declination

        pass


    def get_equation_of_time(self):
        """ calculates the equation of time in minutes """

        oc  = radians(self.oblique_corr)
        gml = radians(self.geomean_long)
        gma = radians(self.geomean_anom)
        ec  = self.earth_eccent

        vary = tan(oc/2)**2

        self.equation_of_time = 4 * degrees(vary * sin(2*gml) - 2 * ec * sin(gma) + \
                                4 * ec * vary * sin(gma) * cos(2 * gml) -           \
                                0.5 * vary * vary * sin(4 * gml) -                  \
                                1.25 * ec * ec * sin(2 * gma))
        
        return self.equation_of_time

    def get_hour_angle_sunrise(self):
        """ calculates the hour hangle of sunrise """

        d   = radians(self.declination)
        lat = radians(self.lat)

        self.hour_angle_sunrise = degrees(arccos((cos(radians(90.833)) /      \
                                  (cos(lat) * cos(d)) - tan(lat) * tan(d))))

        return self.hour_angle_sunrise


    def get_solar_noon(self):
        """ calculats solar noon in (local sidereal time LST)"""

        lat = self.lat
        lon = self.lon
        eot = self.equation_of_time
        tz  = self.tz

        self.solar_noon = (720 - 4 * lon - eot + tz * 60)/1440

        # format this as a time for display purposes (Hours:Minutes:Seconds)
        if self.is_numpy:
            self.solar_noon_time = timedelta(days = self.solar_noon.mean())
        else:
            self.solar_noon_time = timedelta(days = self.solar_noon)
        
        return self.solar_noon_time
    

    def get_sunrise(self):
        """ calculates the time of sunrise"""

        sn = self.solar_noon
        ha = self.hour_angle_sunrise

        self.sunrise = (sn * 1440 - ha * 4)/1440

        # format this as a time for display purposes (Hours:Minutes:Seconds)
        if self.is_numpy:
            self.sunrise_time = timedelta(days = self.sunrise.mean())
        else:
            self.sunrise_time = timedelta(days = self.sunrise)
        
        return self.sunrise_time


    def get_sunset(self):
        """ calculates the time of sunrise"""

        sn = self.solar_noon
        ha = self.hour_angle_sunrise

        
        self.sunset = (sn * 1440 + ha * 4)/1440

        # format this as a time for display purposes (Hours:Minutes:Seconds)
        if self.is_numpy:
            self.sunset_time = timedelta(days = self.sunset.mean())
        else:
            self.sunset_time = timedelta(days = self.sunset)

        return self.sunset_time


    def get_sunlight(self):
        """ calculates amount of daily sunlight in fractional days"""

        self.sunlight = 8 * self.hour_angle_sunrise / (60 * 24)

        # format this as a time for display purposes (Hours:Minutes:Seconds)
        if self.is_numpy:
            self.sunlight_time = timedelta(days = self.sunlight.mean())
        else:
            self.sunlight_time = timedelta(days = self.sunlight)
        
        return self.sunlight_time

    
    def get_true_solar(self):
        """ calculates the true solar time at ref_datetime"""

        lat = self.lat
        lon = self.lon
        eot = self.equation_of_time

        # turn reference datetime into fractional days
        frac_sec = (self.rdt - datetime(self.rdt.year, self.rdt.month, self.rdt.day)).total_seconds() 
        frac_hr  = frac_sec / (60 * 60) + self.tz
        frac_day = frac_hr / 24

        self.frac_day = frac_day

        # now get true solar time
        self.true_solar = (frac_day * 1440 + eot + 4 * lon - 60 * self.tz) % 1440

        # format this as a time for display purposes (Hours:Minutes:Seconds)
        if self.is_numpy:
            self.true_solar_time = timedelta(days = self.true_solar.mean() / (60*24))
        else:
            self.true_solar_time = timedelta(days = self.true_solar / (60*24))

        return self.true_solar_time
        

    def get_hour_angle(self):
        """ calculates the hour angle at ref_datetime"""

        ts = self.true_solar

        # matrix hour_angle calculations
        if self.is_numpy:
            ha = ts
            ha[ha <= 0] = ha[ha <= 0]/4 + 180
            ha[ha >  0] = ha[ha >  0]/4 - 180
            self.hour_angle = ha

        # scalar hour_angle calculations
        else:
            if ts <= 0:
                self.hour_angle = ts/4 + 180
            else:
                self.hour_angle = ts/4 - 180

        return self.hour_angle
        
        
    def get_zenith(self):
        """ calculates solar zenith angle at ref_datetime"""

        d   = radians(self.declination)
        ha  = radians(self.hour_angle)
        lat = radians(self.lat)

        self.zenith = degrees(arccos(sin(lat) * sin(d) + \
                                    cos(lat) * cos(d) * cos(ha)))

        return self.zenith
    

    def get_elevation(self):
        """ calculates solar elevation angle at ref_datetime"""
        
        # perform an approximate atmospheric refraction correction
        
        # matrix hour_angle calculations
        if self.is_numpy:
            e = ((self.zenith * 0) + 90) - self.zenith
            ar = e
            
            ar[ar >  85]     = 0
            ar[ar >   5]     = 58.1 / tan(radians(ar[ar >5])) - 0.07 / tan(radians(ar[ar >5]))**3 + 0.000086 / tan(radians(ar[ar >5]))**5
            ar[ar >  -0.575] = 1735 + ar[ar > -0.575] * (103.4 + ar[ar > -0.575] * (-12.79 + ar[ar > -0.575] * 0.711))
            ar[ar <= -0.575] = -20.772 / tan(radians(ar[ar <= -0.575]))

        # scalar hour_angle calculations
        else:
            e  = 90 - self.zenith
            er = radians(e)
            
            if   e > 85:        ar = 0
            elif e > 5:         ar = 58.1 / tan(er) - 0.07 / tan(er)**3 + 0.000086 / tan(er)**5  
            elif e > -0.575:    ar = 1735 + e * (103.4 + e * ( -12.79 + e * 0.711)) 
            else:               ar = -20.772 / tan(er)
            
        self.elevation_noatmo = e
        self.atmo_refraction  = ar / 3600
        self.elevation        = self.elevation_noatmo + self.atmo_refraction
        
        return self.elevation

    
    def get_azimuth(self):
        """ calculates solar azimuth angle at ref_datetime"""

        lat = radians(self.lat)        
        d   = radians(self.declination)
        ha  = self.hour_angle                
        z   = radians(self.zenith)     

        # matrix azimuth calculations
        if self.is_numpy: pass
            
        if ha > 0:
            self.azimuth = (degrees(arccos(((sin(lat) * cos(z)) - sin(d)) / cos(lat) * sin(z))) + 180) % 360

        else:
            self.azimuth = (540 - degrees(arccos(((sin(lat) * cos(z)) -sin(d))/ (cos(lat) * sin(z))))) % 360

        return self.azimuth


    def get_earth_distance(self):
        """ distance between the earth and the sun at ref_datetime"""

        # temporarily set to simply return 1 AU (in meters)
        self.earth_dist = float(149567871000)

        # this should eventually be a function of datetime for better precision
        
        return self.earth_dist


    def get_inst_irradiance(self):
        """
        calculates instantaneous incoming solar energy in W/m^2

        inst_irradiance is calculated as 
            = sun_surf_radiance *(sun_radius / earth_distance)^2

        then corrected as a function of solar incidence angle
        """

        ed = self.get_earth_distance()

        # calculate irradiance to normal surface at earth distance
        self.inst_irradiance = self.sun_surf_rad * (self.sun_radius / ed)**2

        # correct for surface orientation based on sun angles
        
        return self.inst_irradiance


    def compute_all_scalar(self):
        """ computes and prints all the attributes of this solar object"""

        print("="*50)
        print("latitude, longitude \t{0}, {1}".format(self.lat, self.lon))
        
        print("datetime \t\t{0} (GMT)".format(self.rdt))
        print("abs julian day \t\t{0}\t (day)".format(self.ajd))
        print("abs julian century \t{0}\t (cen)".format(self.ajc))
        print("suns goemean long \t{0}\t (deg)".format(self.get_geomean_long()))
        print("suns goemean anom \t{0}\t (deg)".format(self.get_geomean_anom()))
        print("earth eccentricity \t{0}".format(self.get_earth_eccent()))
        print("suns eq of center \t{0}".format(self.get_sun_eq_of_center()))
        print("suns true long \t\t{0}\t (deg)".format(self.get_true_long()))
        print("suns true anom \t\t{0}\t (deg)".format(self.get_true_anom()))
        print("radiation vector \t{0}\t (AU)".format(self.get_rad_vector()))
        print("suns apparent long \t{0}\t (deg)".format(self.get_app_long()))
        print("earth obliq mean elip \t{0}\t (deg)".format(self.get_oblique_mean_elip()))
        print("earth obliq correction\t{0}\t (deg)".format(self.get_oblique_corr()))
        print("sun right ascension \t{0}\t (deg)".format(self.get_right_ascension()))
        print("solar declination angle {0}\t (deg)".format(self.get_declination()))
        print("equation of time \t{0}\t (min)".format(self.get_equation_of_time()))
        print("hour angle sunrise\t{0}\t (deg)".format(self.get_hour_angle_sunrise()))
        print("")
        print("solar noon \t\t{0}\t (HMS - LST)".format(self.get_solar_noon()))
        print("sunrise \t\t{0}\t (HMS - LST)".format(self.get_sunrise()))
        print("sunset  \t\t{0}\t (HMS - LST)".format(self.get_sunset()))
        print("sunlight durration \t{0}\t (HMS)".format(self.get_sunlight()))
        print("")
        print("true solar time \t{0}\t (HMS - LST)".format(self.get_true_solar()))
        print("hour angle \t\t{0}\t (deg)".format(self.get_hour_angle()))
        print("solar zenith angle \t{0}\t (deg)".format(self.get_zenith()))
        print("solar elevation angle \t{0}\t (deg)".format(self.get_elevation()))
        print("solar azimuth angle \t{0}\t (deg)".format(self.get_azimuth()))

        print("")
        print("irradiance \t\t{0}\t (W/m*m)".format(self.get_inst_irradiance()))
        print("="*50)


    def compute_all_numpy(self):
        """ computes and prints all the attributes of this solar object"""

        print("="*50)
        print("latitude, longitude \t{0}, {1}".format(self.lat, self.lon))
        
        print("datetime \t\t{0} (GMT)".format(self.rdt))
        print("abs julian day \t\t{0}\t (day)".format(self.ajd))
        print("abs julian century \t{0}\t (cen)".format(self.ajc))
        print("suns goemean long \t{0}\t (deg)".format(self.get_geomean_long()))
        print("suns goemean anom \t{0}\t (deg)".format(self.get_geomean_anom()))
        print("earth eccentricity \t{0}".format(self.get_earth_eccent()))
        print("suns eq of center \t{0}".format(self.get_sun_eq_of_center()))
        print("suns true long \t\t{0}\t (deg)".format(self.get_true_long()))
        print("suns true anom \t\t{0}\t (deg)".format(self.get_true_anom()))
        print("radiation vector \t{0}\t (AU)".format(self.get_rad_vector()))
        print("suns apparent long \t{0}\t (deg)".format(self.get_app_long()))
        print("earth obliq mean elip \t{0}\t (deg)".format(self.get_oblique_mean_elip()))
        print("earth obliq correction\t{0}\t (deg)".format(self.get_oblique_corr()))
        print("sun right ascension \t{0}\t (deg)".format(self.get_right_ascension()))
        print("solar declination angle {0}\t (deg)".format(self.get_declination()))
        print("equation of time \t{0}\t (min)".format(self.get_equation_of_time()))
        print("hour angle sunrise\t{0}\t (deg)".format(self.get_hour_angle_sunrise()))
        print("")
        print("solar noon \t\t{0}\t (HMS - LST)".format(self.get_solar_noon()))
        print("sunrise \t\t{0}\t (HMS - LST)".format(self.get_sunrise()))
        print("sunset  \t\t{0}\t (HMS - LST)".format(self.get_sunset()))
        print("sunlight durration \t{0}\t (HMS)".format(self.get_sunlight()))
        print("")
        print("true solar time \t{0}\t (HMS - LST)".format(self.get_true_solar()))
        print("hour angle \t\t{0}\t (deg)".format(self.get_hour_angle()))
        print("solar zenith angle \t{0}\t (deg)".format(self.get_zenith()))
        print("solar elevation angle \t{0}\t (deg)".format(self.get_elevation()))
        print("solar azimuth angle \t{0}\t (deg)".format(self.get_azimuth()))

        print("")
        print("irradiance \t\t{0}\t (W/m*m)".format(self.get_inst_irradiance()))
        print("="*50)

if __name__ == "__main__":

    
    # scalar test
    lat         = 37
    lon         = -76.4
    tz          = -4
    datestamp   = datetime.now()
    
    sc = solar(lat, lon, tz, datestamp)
    
    # numpy array test
    lat         = zeros((1, 2)) + 37
    lon         = zeros((1, 2)) -76.4
    tz          = -4
    datestamp   = datetime.now()
    
    sm = solar(lat, lon, tz, datestamp)







    
        
