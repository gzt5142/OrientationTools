

import numpy as np
from scipy import ndimage


def lambert(S, L):
    """
    :param S: 3-band array of X, Y, and Z components of surface normal vectors. Band is the first index.
    :param L: 3D vector of light direction
    :return: A single band, 2D array representing the shade values according to Lambert Cosine Emission Law
    """
    BV = (S[0] * L[0]) + (S[1] * L[1]) + (S[2] * L[2])
    return BV


def oren_nayer(normal, dirToLight, sigma):
    """
    :param normal: 3-band array of x, y, and z components of surface normal vectors.  Band is first index.
    :param lightDir:3D vector of light direction
    :param sigma: variation of the roughness distribution
    :return: a single-band, 2D array of brightness values.
    """
    # This code adapted from the implementation at
    # https://github.com/ranjak/opengl-tutorial/blob/master/shaders/illumination/diramb_orennayar_pcn.vert

    # we assume for GIS/cartography that the view direction is overhead.
    dirToCamera = np.array([0.0, 0.0, 1.0])

    # float termA = 1.0f - 0.5f * sigma2 / (sigma2 + 0.57f);
    A = 1 - 0.5 * sigma**2 / (sigma**2 + 0.57)

    # float termB = 0.45f * sigma2 / (sigma2 + 0.09f);
    B = 0.45 * sigma**2 / (sigma**2 + 0.09)


    # Angle from surface to light
    # float cosZi = dot(normCamSpace, dirToLight);
    cosZi = lambert(normal, dirToLight)
    ## if <0, light is from 'behind'.  Reflectance is zero.
    cosZi[cosZi<0] = 0

    # Angle from surface to camera
    # float cosZr = dot(normCamSpace, viewDir);
    cosZr = lambert(normal, dirToCamera)

    # float cosAzimuthSinaTanb = (dot(dirToLight, viewDir) - cosZr * cosZi) / max(cosZr, cosZi);
    cosAz_sinA_tanB = (np.dot(dirToLight, dirToCamera) - (cosZr * cosZi)) / np.maximum(cosZr, cosZi)

    # vec4 orenNayarColor = diffuseColor * cosZi * (termA + termB * max(0.0f, cosAzimuthSinaTanb)) *lightIntensity;
    BV = cosZi * (A + B * np.maximum(0, cosAz_sinA_tanB))

    return BV

def lommel_seeliger(normal, dirToLight):
    dirToCamera = np.array([0, 0, 1.0])
    cos_e = lambert(normal, dirToCamera)
    cos_i = lambert(normal, dirToLight)
    BV = np.zeros(cos_i.shape)
    BV[cos_i >= 0] = 1 / (1 + (cos_e[cos_i >= 0] / cos_i[cos_i >= 0]))
    return BV

def mixedLambert(normal, dirToLight):
    w1 = 15
    w2 = 85
    gkernel = np.array([
        [1, 4, 6, 4, 1],
        [4, 16, 24, 16, 4],
        [6, 24, 36, 24, 6],
        [4, 16, 24, 16, 4],
        [1, 4, 6, 4, 1]
    ])
    k = gkernel / np.sum(gkernel.ravel())
    blurred = np.zeros(normal.shape)
    blurred[0] = ndimage.convolve(normal[0], k)
    blurred[1] = ndimage.convolve(normal[1], k)
    blurred[2] = 0# ndimage.convolve(normal[2], k)
    mag = np.sqrt(blurred[0] ** 2 + blurred[1] ** 2 + blurred[2] ** 2)
    mag[mag==0] = 1
    blurred = blurred / mag

    bv1= normal[0] * blurred[0] + normal[1]*blurred[1]
    bv2 = lambert(normal, dirToLight)
    BV = bv1*(w1/(w1+w2)) + bv2*(w2/(w1+w2))
    BV[BV<0] = 0
    return BV
