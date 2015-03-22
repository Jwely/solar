

from datetime import datetime

class solar:
    """
    Object class for handling solar calculations

    It requires a physical location on the earth and a datetime object
    """

    def __init__(self, lat, lon):
        """initializes critical  spatial information"""

        # attributes
        self.lat    = float(lat)        # lattitude - float
        self.lon    = float(lon)        # longitude - float
        self.rdt    = []                # exact datetime of solar object
        self.tz     = 0.0               # time zone (from GTM) - float

        pass


    def set_datetime(self, date_time_obj, fmt = False):
        """
        sets the critical time information

        accepts datetime objects or a datetime string with format
        """

        # if input is datetime_obj set it
        if isinstance(date_time_obj, datetime):
            self.rdt = date_time_obj
            return

        if isinstance(date_time_obj, str) and isinstance(fmt, str):
            self.rdt = datetime.strptime(date_time_obj,fmt)
            return
        else:
            raise Exception("bad datetime!")
    
        
    def abs_julian_day(self):
        """ calculates the ABSOLUTE julian day"""

        pass
    

    def timezone(self, time_shift):
        """ amends datetime object with timezone information"""

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

        
