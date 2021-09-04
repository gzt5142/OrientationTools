# -*- coding: utf-8 -*-
#
# License: GNU GPL v3 -- See 'LICENCE' file.
#
#
import numpy as np
from math import sin, cos, pi


def lightVector(Az, El):
    """
    :param Az: Compass direction (from north, increasing clockwise).  In degrees [0-359]
    :param El: Elevation above idealized horizon. In degrees [0-90].  0=horizon, 90=overhead
    :return: normalized vector (as 1-D numpy array) representing light direction
    """
    Az_rad = float(Az) * (pi / 180)
    El_rad = float(El) * (pi / 180)

    L = np.array([
        sin(Az_rad) * cos(El_rad),  # X
        cos(Az_rad) * cos(El_rad),  # Y
        sin(El_rad)                 # Z
    ])
    ## In theory, this vector is already unit length one.  But let's normalize anyway...
    L /= np.sqrt((L**2).sum())
    return L

