

from datetime import datetime, timedelta

class solar:
    """
    Object class for handling solar calculations

    It requires a physical location on the earth and a datetime object
    """

    def __init__(self, lat, lon):
        """initializes critical  spatial information"""

        # attributes
        self.lat    = float(lat)        # lattitude (E positive)- float
        self.lon    = float(lon)        # longitude (N positive)- float
        self.rdt    = []                # ref_datetime. Precise GMT time of solar obj
        self.tz     = 0.0               # timezone, offset from GMT - float
        self.adj    = 0.0               # absolute julian day
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

        self.abs_julian_day()
        return


    def set_timezone(self, GMT_hour_offset):
        """ amends the datetime object with time zone information"""

        # sets time zone in a way to prevent accidental repeated offsetting
        if self.tz != GMT_hour_offset:
            time_del = timedelta(hours = (GMT_hour_offset - self.tz))
            self.rdt = self.rdt + time_del
            self.tz  = GMT_hour_offset

        self.abs_julian_day()
        return 

        
    def abs_julian_day(self):
        """
        calculates the ABSOLUTE decimal julian day of ref_datetime

        This is invoked whenever set_datetime or set_timezone is invoked
        """

        # uses the reference day of january 1st 2000
        jan_1st_2000_jd   = 2451545
        jan_1st_2000      = datetime(2000,1,1,12,0,0)

        time_del = self.rdt - jan_1st_2000
        self.ajd = float(jan_1st_2000_jd) + float(time_del.total_seconds())/86400
        pass
    

    def eccentricity(self):
        """ calculates earths orbital eccentricity at ref_datetime"""

        pass


    def distance(self):
        """ distance between the earth and the sun at ref_datetime"""
        

    def rad_vector(self):
        """ incident radiation vector to surface at ref_datetime"""

        pass


    def declination(self):
        """ solar declination angle at ref_datetime"""

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


    def light_duration(self):
        """ calculates durration of sunlight on day of ref_datetime"""

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


if __name__ == "__main__":

    s = solar(37.0349, -76.3601)
    s.set_datetime("20150322120000","%Y%m%d%H%M%S")
    s.set_timezone(-5)
    print("abs julian day:",s.ajd)









    
        
