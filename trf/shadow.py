
import numpy as np
from . import utils


def shadowLine(d, az, el, cellwidth):
    returnArray = np.zeros(d.shape)
    lightVector = utils.lightVector(az, el)

    # X,Y as geographic coordinates are not the same as X,Y in numpy arrays.
    # switch to row/col names for numpy manipulations to save our sanity.
    lineMax = returnArray.shape[0] - 1
    colMax = returnArray.shape[1] - 1

    delta_col = lightVector[0]
    delta_row = -lightVector[1]
    deltha_ht = lightVector[2]

    shadowLine = np.pad(d, 1, 'edge')

    # The counters increment or decrement based on the dominant light direction.
    # I am not fiddling with the math to get those iterables .... Using if statements to run an unrolled
    # bit of code depending on dominant lighting direction.  This is ugly, but works.
    # TODO: generalize the formula to use iterables with the right sign.
    if (az <= 45.0 or az > 315.0):
        dH_drow = abs(deltha_ht / delta_row) * cellwidth
        p = abs(delta_col / delta_row)

        for line in range(1, lineMax):
            dem = d[line - 1, :]
            n = shadowLine[line - 1, 1:-1]
            sl = shadowLine[line, 1:-1]
            a = returnArray[line - 1]
            if delta_col <= 0:
                nw = shadowLine[line - 1, :-2]
            else:
                nw = shadowLine[line - 1, 2:]
            tstHeights = (n * (1 - p)) + (nw * (p)) - dH_drow

            sl[dem < tstHeights] = tstHeights[dem < tstHeights]
            a[dem < tstHeights] = tstHeights[dem < tstHeights] - dem[dem < tstHeights]

    if (225 < az <= 315):
        dH_dcol = abs(deltha_ht / delta_col) * cellwidth
        p = abs(delta_row / delta_col)

        for col in range(1, colMax):
            dem = d[:, col - 1]
            w = shadowLine[1:-1, col - 1]
            sl = shadowLine[1:-1, col]
            a = returnArray[:, col - 1]
            if delta_row <= 0:
                nw = shadowLine[:-2, col - 1]
            else:
                nw = shadowLine[2:, col - 1]
            tstHeights = (w * (1 - p)) + (nw * (p)) - dH_dcol

            sl[dem < tstHeights] = tstHeights[dem < tstHeights]
            a[dem < tstHeights] = tstHeights[dem < tstHeights] - dem[dem < tstHeights]

    if (135 < az <= 225):
        dH_drow = abs(deltha_ht / delta_row) * cellwidth
        p = abs(delta_col / delta_row)

        for line in range(lineMax, 0, -1):
            dem = d[line, :]
            s = shadowLine[line + 2, 1:-1]
            sl = shadowLine[line + 1, 1:-1]
            a = returnArray[line]
            if delta_col <= 0:
                sw = shadowLine[line + 2, :-2]
            else:
                sw = shadowLine[line + 2, 2:]

            tstHeights = (s * (1 - p)) + (sw * (p)) - dH_drow

            sl[dem < tstHeights] = tstHeights[dem < tstHeights]
            a[dem < tstHeights] = tstHeights[dem < tstHeights] - dem[dem < tstHeights]

    if (45 < az <= 135):
        dH_dcol = abs(deltha_ht / delta_col) * cellwidth
        p = abs(delta_row / delta_col)

        for col in range(colMax, 0, -1):
            dem = d[:, col]
            e = shadowLine[1:-1, col + 2]
            sl = shadowLine[1:-1, col + 1]
            a = returnArray[:, col - 1]
            if delta_row <= 0:
                ne = shadowLine[:-2, col + 2]
            else:
                ne = shadowLine[2:, col + 2]
            tstHeights = (e * (1 - p)) + (ne * (p)) - dH_dcol

            sl[dem < tstHeights] = tstHeights[dem < tstHeights]
            #a[dem < tstHeights] = 1
            a[dem < tstHeights] = tstHeights[dem < tstHeights] - dem[dem < tstHeights]

    return returnArray

