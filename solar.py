
__author__ = "Jeffry Ely, Jeff.ely.08@gmail.com"

from datetime import datetime, timedelta
import math

class solar:
    """
    Object class for handling solar calculations

    It requires a physical location on the earth and a datetime object

    Units used by this class:
        angle    = radians
        distance = meters
        energy   = watts or joules
        time     = mostly in datetime objects. typically specified
    """

    def __init__(self, lat, lon, time_zone, date_time_obj, fmt = False):
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
        
        pass


    def set_datetime(self, date_time_obj, fmt = False):
        """
        sets the critical time information

        accepts datetime objects or a datetime string with format
        """

        # if input is datetime_obj set it
        if isinstance(date_time_obj, datetime):
            self.rdt = date_time_obj

        if isinstance(date_time_obj, str) and isinstance(fmt, str):
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

        return self.ajd

    
    def set_location(self, lat, lon):
        """ sets spatial information"""

        self.lat    = float(lat)        # lattitude (E positive)- float
        self.lon    = float(lon)        # longitude (N positive)- float
        self.lat_r  = math.radians(lat) # lattitude in radians
        self.lon_r  = math.radians(lon) # longitude in radians
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


    def earth_eccentricity(self):
        """
        Calculates the precise eccentricity of earths orbit at ref_datetime

        requires use of absolute julian century
        """
        
        self.earth_eccent = 0.016708634 - self.ajc*(4.2037e-5 + 1.267e-7 * self.ajc)
        return self.earth_eccent
    

    def earth_distance(self):
        """ distance between the earth and the sun at ref_datetime"""

        # temporarily set to simply return 1 AU (in meters)
        self.earth_dist = float(149567871000)

        # this should eventually be a function of datetime for better precision
        
        return self.earth_dist


    def inst_irradiance(self):
        """
        calculates instantaneous incoming solar energy in W/m^2

        inst_irradiance is calculated as 
            = sun_surf_radiance *(sun_radius / earth_distance)^2

        then corrected as a function of solar incidence angle
        """

        ed = self.earth_distance()

        # calculate irradiance to normal surface at earth distance
        self.inst_irradiance = self.sun_surf_rad * (self.sun_radius / ed)**2

        # correct for surface orientation based on sun angles
        

        return self.inst_irradiance


    def mean_long(self):
        """ calculates mean longitude of the sun"""

        # L  = mean anom + lon of ascending node + arg of periapsis

        
    
    def rad_vector(self):
        """ incident radiation vector to surface at ref_datetime"""

        pass


    def declination(self):
        """ solar declination angle at ref_datetime"""

        pass


    def zenith(self):
        """ calculates solar zenith angle at ref_datetime"""

        pass


    def elevation(self):
        """ calculates solar elevation angle at ref_datetime"""

        pass


    def azimuth(self):
        """ calculates solar azimuth angle at ref_datetime"""

        pass


    def atmo_refraction(self):
        """ approximates atmospheric refraction at ref_datetime"""

        pass
    
    
    def all_daily(self):
        """
        calculates all of the daily stats for ref_datetime

        executes methods:
            noon
            sunrise
            sunset
            daily_irradiance
            daily_light
        """

        pass


    def noon(self):
        """solar noon datetime obj for  day of ref_datetime"""

        pass


    def sunrise(self):
        """calculates datetime obj for sunrise on day of ref_datetime"""
        

        pass


    def sunset(self):
        """calculates datetime obj for sunset on day of ref_datetime"""

        pass


    def daily_irradiance(self):
        """ approximates the days quantity of incident solar energy in joules/m^2"""

        pass
    

    def daily_light(self):
        """calculates durration of sunlight on day of ref_datetime"""

        pass


if __name__ == "__main__":

    lat         = 40
    lon         = -105
    tz          = -7
    datestamp   = "20030101120000"
    fmt         = "%Y%m%d%H%M%S"

    s = solar(lat, lon, tz, datestamp, fmt)
    
    print("abs julian day:",s.ajd)
    print("irradiance:",s.inst_irradiance())
    print("eccentricity:",s.earth_eccentricity())









    
        
