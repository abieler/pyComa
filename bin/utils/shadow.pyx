import numpy as np
cimport numpy as np
cimport cython

cdef extern from "math.h":
    double sqrt(double x)

@cython.boundscheck(False)
@cython.wraparound(False)
def calculate_shadow(np.ndarray[np.float_t, ndim=3] triangles, np.ndarray[np.float_t, ndim=2] n_hat, np.ndarray[np.float_t, ndim=1] r, np.ndarray[np.float_t, ndim=1] p0, int nTriangles, int currIndex):
    # return 0 if in shadow and 1 if not
    cdef int i = 0
    cdef int j = 0
    cdef int k = 0
    cdef float a = 0.0
    cdef float b = 0.0
    cdef float rI
    cdef float dot_uv, dot_uu, dot_vv, dot_wv, dot_wu
    cdef float divisor, sI, tI
    cdef np.ndarray[np.float_t, ndim=1] pI = np.zeros(3)
    cdef np.ndarray[np.float_t, ndim=1] u = np.zeros(3)
    cdef np.ndarray[np.float_t, ndim=1] v = np.zeros(3)
    cdef np.ndarray[np.float_t, ndim=1] w = np.zeros(3)
    
    for i in range(nTriangles):
        a = 0.0
        b = 0.0
        for j in range(3):
            a = a + n_hat[i,j] * (triangles[i,0,j] - p0[j])
            b = b + n_hat[i,j] * r[j]
        println("a:", a)
        println("b:", b)
        if (b != 0.0) and (a != 0.0):
            rI = a / b
            if rI >= 0:
                dot_uv = 0.0
                dot_uu = 0.0
                dot_vv = 0.0
                dot_wu = 0.0
                dot_wv = 0.0
                for k in range(3):
                    pI[k] = p0[k] + rI * r[k]
                    u[k] = triangles[i,1,k] - triangles[i,0,k]
                    v[k] = triangles[i,2,k] - triangles[i,0,k]
                    w[k] = pI[k] - triangles[i,0,k]
                    
                    dot_uv = dot_uv + u[k]*v[k]
                    dot_uu = dot_uu + u[k]*u[k]
                    dot_vv = dot_vv + v[k]*v[k]
                    dot_wu = dot_wu + w[k]*u[k]
                    dot_wv = dot_wv + w[k]*v[k]
                
                divisor = dot_uv*dot_uv - dot_uu * dot_vv
                sI = (dot_uv*dot_wv - dot_vv*dot_wu) / divisor
                tI = (dot_uv*dot_wu - dot_uu*dot_wv) / divisor
                
                if ((sI >= 0.00) and (tI >= 0.00) and (sI + tI < 1.0)):
                    if i != currIndex:
                        return 0
    return 1
