import numpy as np
cimport numpy as np
cimport cython


cdef extern from "math.h":
    double sqrt(double x)


@cython.boundscheck(False)
@cython.wraparound(False)
def createRay(int iDim, np.ndarray[np.float_t, ndim=1] rRay, np.ndarray[np.float_t, ndim=1] p):

    cdef float dr
    cdef float distance
    xTravel = []
    distance = sqrt(rRay[0]*rRay[0] + rRay[1]*rRay[1] + rRay[2]*rRay[2])

    while ((distance < 1e8) and (distance > 2000)):
        xTravel.append((rRay[0], rRay[1], rRay[2]))

        if distance < 10000:
            dr = distance / 40.0
            if distance < 4000:
                dr = distance / 100.0
                if distance < 2500:
                    dr = distance / 250.0
                    if distance < 2030:
                        dr = 1
        else:
            dr = distance / 10.0

        rRay[0] = rRay[0] + p[0] * dr
        rRay[1] = rRay[1] + p[1] * dr
        rRay[2] = rRay[2] + p[2] * dr
        distance = sqrt(rRay[0]*rRay[0] + rRay[1]*rRay[1] + rRay[2]*rRay[2])

    return xTravel
