# -*- coding: utf-8 -*-
#
# License: GNU GPL v3 -- See 'LICENCE' file.
#
#
import numpy as np
import csv
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


def lightList(f):
    """
    :param f: File name/path to the sky config file.
    :return: list of light directions and their weights.
    """
    lts = []
    with open(f, "r") as skyConfigFile:
        csvReadFile = csv.reader(skyConfigFile)
        for line in csvReadFile:
            # brute force -- if this line does not have three comma-separated values, we skip.
            if (len(line) < 3):
                continue
            if line[0].startswith("Format"): #special case... a 'comment' line in the header happens to have two commas in it.
                continue
            lts.append(line)
    return lts


def normalize(bumpMap):
    # "...Casting {} to range -1 to 1".format(bumpMap.dtype)
    if bumpMap.dtype == 'uint8':
        return ((bumpMap / (2 ** 8 - 1)) * 2) - 1
    if bumpMap.dtype == 'uint16':
        return ((bumpMap / (2 ** 16 - 1)) * 2) - 1
    if bumpMap.dtype == 'float32':
        return (bumpMap - 0.5) * 2
    if bumpMap.dtype == 'float64':
        return (bumpMap - 0.5) * 2
    return bumpMap