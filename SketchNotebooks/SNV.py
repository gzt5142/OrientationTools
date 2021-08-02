from math import pi as PI
from math import sin, cos
import numpy as np
import csv


def lightVector(Az, El):
    Az_rad = float(Az) * (PI / 180)
    El_rad = float(El) * (PI / 180)

    L = np.array([
        sin(Az_rad) * cos(El_rad),  # X
        cos(Az_rad) * cos(El_rad),  # Y
        sin(El_rad)  # Z
    ])
    ## In theory, this vector is already unit length one.  But let's normalize anyway...
    L /= np.sqrt((L**2).sum())
    return L


def castShadows(d, az, el, cellwidth):
    A = np.ones(d.shape)
    S = lightVector(az, el)

    lineMax=A.shape[0]-1
    colMax=A.shape[1]-1

    del_col= S[0]
    del_row = -S[1]
    del_H = S[2]

    shadowLine = np.pad(d, 1, 'edge')

    if (az <= 45.0 or az > 315.0):
        #print("Light from North")
        dH_drow =  abs(del_H / del_row) * cellwidth
        p = abs(del_col / del_row)

        for line in range(1, lineMax):
            dem = d[line-1,:]
            n=shadowLine[line-1,1:-1]
            sl=shadowLine[line,1:-1]
            a = A[line-1]
            if del_col <= 0:
                nw=shadowLine[line-1,:-2]
            else:
                nw=shadowLine[line-1,2:]
            tstHeights = (n * (1-p) ) + (nw * ( p)) - dH_drow

            sl[dem < tstHeights] = tstHeights[dem < tstHeights]
            a[dem < tstHeights] = 0

    if (225 < az <= 315):
        #print("Light from West")
        dH_dcol =  abs(del_H / del_col) * cellwidth
        p = abs(del_row / del_col)

        for col in range(1, colMax):
            dem = d[:,col-1]
            w=shadowLine[1:-1, col-1]
            sl=shadowLine[1:-1,col]
            a = A[:,col-1]
            if del_row <= 0:
                nw=shadowLine[:-2, col-1]
            else:
                nw=shadowLine[2:,col-1]
            tstHeights = (w * (1-p) ) + (nw * ( p)) - dH_dcol

            sl[dem < tstHeights] = tstHeights[dem < tstHeights]
            a[dem < tstHeights] = 0

    if (135 < az <= 225):
        #print("Light from South")
        dH_drow =  abs(del_H / del_row) * cellwidth
        p = abs(del_col / del_row)

        for line in range(lineMax, 0, -1):
            dem = d[line,:]
            s=shadowLine[line+2,1:-1]
            sl=shadowLine[line+1,1:-1]
            a = A[line]
            if del_col <= 0:
                sw=shadowLine[line+2,:-2]
            else:
                sw=shadowLine[line+2,2:]

            tstHeights = (s * (1-p) ) + (sw * ( p)) - dH_drow

            sl[dem < tstHeights] = tstHeights[dem < tstHeights]
            a[dem < tstHeights] = 0

    if (45 < az <= 135):
        #print("Light from East")
        dH_dcol =  abs(del_H / del_col) * cellwidth
        p = abs(del_row / del_col)

        for col in range(colMax, 0, -1):
            dem = d[:,col]
            e=shadowLine[1:-1, col+2]
            sl=shadowLine[1:-1,col+1]
            a = A[:,col-1]
            if del_row <= 0:
                ne=shadowLine[:-2, col+2]
            else:
                ne=shadowLine[2:,col+2]
            tstHeights = (e * (1-p) ) + (ne * ( p)) - dH_dcol

            sl[dem < tstHeights] = tstHeights[dem < tstHeights]
            a[dem < tstHeights] = 0
            
    return A


def surfaceNormals(elevArray, cellSize):
    A = np.pad(elevArray, 1, 'edge')  # adds a 1-px border; covers edge cases.  Will be trimmed in the slicing manipulations.

    a = A[0:-2, 0:-2]
    b = A[0:-2, 1:-1]
    c = A[0:-2, 2:]

    d = A[1:-1, 0:-2]
    f = A[1:-1, 2:]

    g = A[2:, 0:-2]
    h = A[2:, 1:-1]
    i = A[2:, 2:]

    Z = np.ones(elevArray.shape, dtype=A.dtype)

    dx = ((a + 2*d + g) - (c + 2*f + i)) / ( 8 * cellSize)
    dy = ((g + 2*h + i) - (a + 2*b + c)) / ( 8 * cellSize)

    mag = np.sqrt(dx ** 2 + dy ** 2 + Z ** 2)
    sn = np.stack([dx / mag, dy / mag, Z / mag], 2)
    return sn


def lambert(S, L):
    hillshadeArray = (S[:,:,0] * L[0]) + (S[:,:,1] * L[1]) + (S[:,:,2] * L[2])
    return hillshadeArray



def lightList(f):
    lts = []
    with open(f, "r") as skyConfigFile:
        csvReadFile = csv.reader(skyConfigFile)
        for line in csvReadFile:
            # brute force -- if this line does not have three comma-separated values, we skip it.
            if (len(line) == 0) or (len(line) < 3):
                continue
            if line[0].startswith("Format"):  # special case... a comment line in the header has 2 commas in it.
                continue
            lts.append(line)
    return lts


