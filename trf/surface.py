
import numpy as np
from math import sqrt


def gradient(A, cellwidth, method="N4"):
    """
    :param A: The array of heights. For terrains, this is the DEM.
    :param cellwidth: cell height and width.  This should be in same units as the heights encoded in A
    :param m: Which method to estimate gradient? A list of valid method identifiers is below.
    :return: A two-element array holding rasters for nabla X and nabla Y for the given height field.
    """
    #
    A = np.pad(A, 1, 'edge')  # adds a 1-px border; covers edge cases.  Will be trimmed in the slicing manipulations.

    # Create slices of the height array, shifted in each direction.
    a = A[0:-2, 0:-2]
    b = A[0:-2, 1:-1]
    c = A[0:-2, 2:]

    d = A[1:-1, 0:-2]
    e = A[1:-1, 1:-1]
    f = A[1:-1, 2:]

    g = A[2:, 0:-2]
    h = A[2:, 1:-1]
    i = A[2:, 2:]

    #
    # +-------+-------+-------+
    # |       |       |       |
    # |   a   |   b   |   c   |
    # |       |       |       |
    # +-------+-------+-------+
    # |       |       |       |
    # |   d   |   e   |   f   |
    # |       |       |       |
    # +-------+-------+-------+
    # |       |       |       |
    # |   g   |   h   |   i   |
    # |       |       |       |
    # +-------+-------+-------+
    #
    # For a given cell at 'e', we can now compute differences among its neighbors by referring to the arrays
    # at a, b, c, etc.  Subtracting d from f, for example, yields an array where the cell values are the difference
    # between the west and east orthogonal cells for all cells in the array.

    # I'm using these different views into the array rather than convolutions, because it is marginally faster
    # than importing ndimage and using ndimage.convolve with the relevant kernels.  Slicing in this way also
    # avoids making copies of the DEM array.

    m = method.upper()
    # just in case we get handed a bogus method code:
    if m not in ['FFD', 'SimpleD', 'N4', '2FD', 'N82', '3FDWRSD', 'N8R', '3FDWRD', 'N8E', '3FD']:
        m = 'N4'

    # Various methods to compute the partial derivatives for change in z with respect to x and y
    # are available.  dz/dx is dz_dx and dz/dy is dz_dy

    # N4 - Four Neighbors.  This is a second order, finite difference (i.e. the Rook's Case).
    #      TAPES-G calls this "2FD".
    #      Fleming and Hoffer, 1979; Ritter, 1987; Zevenbergen and Thorne (1987); O'Neill & Mark (1987)
    #
    # This is an easy, intuitive case, and is our default if no method specified (or if the method given is
    # not in our list of implemented methods).
    if (m == "N4" or m == "2FD"):
        dz_dx = (f - d) / (2 * cellwidth)
        dz_dy = (b - h) / (2 * cellwidth)

    # N82 - 8 neighbors, weighted. Third order finite difference weighted by reciprocal of squared distance
    #      differential weights (Horn, 1981); Sobel operator (Richards, 1986).
    # TAPES-G name is 3FDWRSD
    # This is the method Esri uses for computing gradients (and from gradients, slope and aspect).
    if (m == "N82" or m == "3FDWRSD"):
        dz_dx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * cellwidth)
        dz_dy = ((a + 2 * b + c) - (g + 2 * h + i)) / (8 * cellwidth)

    # SimpleD - Simple Difference
    if (m == "SimpleD"):
        dz_dx = (e - h) / cellwidth
        dz_dy = (e - d) / cellwidth

    #
    # FFD - Finite Frame Difference
    if (m == "FFD"):
        dz_dx = (a - g + c - i) / (4 * cellwidth)
        dz_dy = (i - g + c - a) / (4 * cellwidth)

    # N8R - 8 neighbors, weighted. Third order finite difference weighted by reciprocal  distance.
    #      Unwin, 1981 -  3FDWRD
    #
    if (m == "N8R" or m == "3FDWRD"):
        r = sqrt(2)
        dz_dx = ((c + r * f + i) - (a + r * d + g)) / ((4 + 2 * r) * cellwidth)
        dz_dy = ((a + r * b + c) - (g + r * h + i)) / ((4 + 2 * r) * cellwidth)

    # N8E - 8 neighbors, even weighting. "Queen's case": Equivalent to fitting a second order trend surface.
    #      Identical weights. TAPES-G calls this "3FD" -- "Third order finite difference"
    # 	    Horn, 1981; Heerdegen and Beran, 1982;  Wood, 1996
    if (m == "N8E" or m == "3FD"):
        dz_dx = ((c + f + i) - (a + d + g)) / (6 * cellwidth)
        dz_dy = ((a + b + c) - (g + h + i)) / (6 * cellwidth)

    return [dz_dx, dz_dy]


def normals_by_method(DEM, cellwidth, m):
    # Z = np.ones(DEM.shape, dtype=np.double)
    [dx, dy] = gradient(DEM, cellwidth, method=m)

    # save a few cycles by using the constant '1' instead of Z**2, since Z is always ones before normalizing
    mag = np.sqrt(dx ** 2 + dy ** 2 + 1)
    return np.stack([-dx / mag, -dy / mag, 1 / mag], 0)

