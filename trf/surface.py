
import numpy as np
from math import sqrt

from . import utils

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



def bump_extent_by_mask(maskArray, bumpMapTile):
    """
    :param maskArray:  boolean array of the cells where the bumpmap is to be applied.
    :param bumpMapTile: the bumpmap to tile over the extent of maskarray.
    :return:
    """
    (x, y) = maskArray.shape
    returnBumpMap = np.zeros((3, x, y))
    returnBumpMap[2,:,:] = 1.0 # all vectors in the returned bump-map are now [0,0,1]


    # IMPORTANT NOTE:  An RGB image is read in to an array of shape (row, column, band).
    # In manipulating multi-band rasters in arcpy, the band is the first index.  Needs to be
    # re-shaped to (band, row, col).  We'll do that after tiling the image...

    x_tiles = (x // bumpMapTile.shape[0]) + 2
    y_tiles = (y // bumpMapTile.shape[1]) + 2
    img = np.tile(bumpMapTile, (x_tiles, y_tiles, 1))
    center_x = img.shape[0] // 2
    center_y = img.shape[1] // 2
    Bm = img[center_x - (x // 2):center_x + (x // 2) + 1, center_y - (y // 2):center_y + (y // 2) + 1, :]
    Bm = Bm[0:x, 0:y, 0:3] # discard alpha channel, if present
    B = utils.normalize(Bm)
    B = np.moveaxis(B, 2, 0) ##<<<<< change from (row,col,band) to (band,row,col)
    n = np.broadcast_to(maskArray, (3, x, y))
    returnBumpMap[n] = B[n]
    return  returnBumpMap


def applyBumpMap(S, B):
    """
    :param S: Surface normal vector array.  (band, row, column)
    :param B: Bump map surface normal vectors (band, row, column)
    :return: Bumpified surface normal vectors for surface
    """
    Uz = S[0,:, :] / S[2, :, :]
    mag = np.sqrt(1 + Uz ** 2)
    U = np.stack([1 / mag, np.zeros((S.shape[1], S.shape[2])), -Uz / mag], 0)

    Vz = S[1, :, :] / S[2, :, :]
    mag = np.sqrt(1 + Vz ** 2)
    V = np.stack([np.zeros((S.shape[1], S.shape[2])), 1 / mag, -Vz / mag], 0)

    # N is a straightforward change-of-basis calculation.  We are going from tangent space
    # to world space.  Tangent space is defined by U, V, and S vectors:
    #    - U is i-hat (vector representing one unit in "x" direction;
    #    - V is j-hat (vector representing one unit in "y" direction;
    #    - S is the original surface normal, a.k.a. k-hat -- one unit in the "z" direction.
    # Dot each of these vectors with the bump-map vector (which is in tangent space)
    # to get the i,j,k components of that same vector in world space:
    Nx = (U[0, :, :] * B[0, :, :]) + (V[0,:, :] * B[1,:, :]) + (S[0, :, :] * B[2, :, :])
    Ny = (U[1, :, :] * B[0, :, :]) + (V[1,:, :] * B[1,:, :]) + (S[1, :, :] * B[2, :, :])
    Nz = (U[2, :, :] * B[0, :, :]) + (V[2,:, :] * B[1,:, :]) + (S[2, :, :] * B[2, :, :])

    mag = np.sqrt(Nx ** 2 + Ny ** 2 + Nz ** 2)  # In theory, N is already normalized.  This guarantees it.
    # N is the surface normal vector S, with bump map B applied to it.
    return np.stack([Nx / mag, Ny / mag, Nz / mag], 0)