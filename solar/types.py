from datetime import datetime
from typing import Union, List

import numpy as np


# condense verbose type hints for vectorized or scalar inputs.
class Constants(object):
    """ Constants (on any normal human time scale) """

    def __init__(self):
        self.sun_surf_rad = 63156942.6  # radiation at suns surface (W/m^2)
        self.sun_radius = 695800000.  # radius of the sun in meters
        self.orbital_period = 365.2563630  # num of days it takes earth to revolve
        self.altitude = -0.01448623  # altitude of center of solar disk

        # solar energy constants
        self.earth_radius = 6371000  # radius of simple spherical earth (m)
        self.atm_height = 9000       # effective height of atmosphere

        # photosynthetic constants
        self.etta_photon = 4.56  # micro mol / J
        self.etta_par = 0.368  # fraction of 5800K black body solar radiation that is photosynthetically active (W/W)


CONSTANTS = Constants()
FlexNum = Union[float, int, np.ndarray]
FlexDate = Union[datetime, List[datetime]]