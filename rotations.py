#!/Users/abieler/Library/Enthought/Canopy_64bit/User/bin/python
#!/Users/abieler/anaconda/bin/python2.7
#!/usr/bin/env python
from __future__ import division
import sys
#sys.path.append("/Users/abieler/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages")
#for p in sys.path:
#    print p
import spice
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plotVector(x0, x_hat, ax):
    x_hat = x_hat / np.sqrt(np.sum(x_hat**2))
    x = x0 + x_hat
    xplt = [x0[0], x[0]]
    yplt = [x0[1], x[1]]
    zplt = [x0[2], x[2]]
    ax.plot(xplt, yplt, zplt)


def rotateVector(x0, phaseAngle=0, latitude=0):
    theta = phaseAngle / 180 * np.pi
    latitude = latitude / 180 * np.pi

    xr = np.array([1, 0, 0])
    yr = np.array([0, 1, 0])
    zr = np.array([0, 0, 1])

    xNew = spice.vrotv(x0, zr, theta)
    yNew = spice.vrotv(yr, zr, theta)

    xNew = spice.vrotv(xNew, yNew, latitude)

    return np.array(xNew)


def rotateCoordinateSystem(phaseAngle=0, latitude=5, el=0, az=0):

    theta = phaseAngle / 180 * np.pi
    latitude = latitude / 180 * np.pi
    el = el / 180 * np.pi
    az = az / 180 * np.pi

    xr = np.array([1, 0, 0])
    yr = np.array([0, 1, 0])
    zr = np.array([0, 0, 1])

    ########################################
    # rotate around z axis for phase angle
    ########################################
    xr = spice.vrotv(xr, zr, theta)
    yr = spice.vrotv(yr, zr, theta)
    zr = spice.vrotv(zr, zr, theta)

    rotAxis1 = spice.vcrss(xr, -np.array(zr))

    xr = spice.vrotv(xr, rotAxis1, latitude)
    yr = spice.vrotv(yr, rotAxis1, latitude)
    zr = spice.vrotv(zr, rotAxis1, latitude)

    xr = np.array(xr)
    yr = np.array(yr)
    zr = np.array(zr)

    xr = spice.vrotv(xr, -zr, az)
    yr = spice.vrotv(yr, -zr, az)

    xr = np.array(xr)
    yr = np.array(yr)

    xr = spice.vrotv(xr, -yr, el)
    zr = spice.vrotv(zr, -yr, el)

    xr = np.array(xr)
    zr = np.array(zr)

    return xr, yr, zr


def rotateCoordinateSystem2(phaseAngle=0, latitude=5, alpha=0, beta=0, gamma=0):

    theta = phaseAngle / 180 * np.pi
    latitude = latitude / 180 * np.pi
    alpha = alpha / 180 * np.pi
    beta = beta / 180 * np.pi
    gamma = gamma / 180 * np.pi

    xr = np.array([1, 0, 0])
    yr = np.array([0, 1, 0])
    zr = np.array([0, 0, 1])

    ########################################
    # rotate around z axis for phase angle
    ########################################
    xr = spice.vrotv(xr, zr, theta)
    yr = spice.vrotv(yr, zr, theta)
    zr = spice.vrotv(zr, zr, theta)

    rotAxis1 = spice.vcrss(xr, -np.array(zr))

    xr = spice.vrotv(xr, rotAxis1, latitude)
    yr = spice.vrotv(yr, rotAxis1, latitude)
    zr = spice.vrotv(zr, rotAxis1, latitude)

    xr = np.array(xr)
    yr = np.array(yr)
    zr = np.array(zr)

    # rotate around x-axis by alpha
    yr = spice.vrotv(yr, -xr, alpha)
    zr = spice.vrotv(zr, -xr, alpha)
    yr = np.array(yr)
    zr = np.array(zr)

    # rotate around y-axis by beta
    xr = spice.vrotv(xr, yr, beta)
    zr = spice.vrotv(zr, yr, beta)
    xr = np.array(xr)
    zr = np.array(zr)

    # rotate around z-axis by gamma
    xr = spice.vrotv(xr, zr, gamma)
    yr = spice.vrotv(yr, zr, gamma)
    xr = np.array(xr)
    yr = np.array(yr)

    return xr, yr, zr


def createRotationMatrix(i,j,k):
    C = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]
    Cr = [i, j, k]
    R = np.zeros((3, 3))
    for l in range(3):
        for m in range(3):
            R[l][m] = C[l].dot(Cr[m])
    return R


'''
phaseAngle = 180
latitude = 0
az = 0
el = 20

alpha = 20
beta = 0
gamma = 0


i,j,k = rotateCoordinateSystem2(phaseAngle,latitude,alpha,beta,gamma)


C = [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])]
Cr = [i,j,k]

R = np.zeros((3,3))
for l in range(3):
    for m in range(3):
        R[l][m] = C[l].dot(Cr[m])
        

print R

print R.dot(np.array([1,0,0]))


x0 = np.array([-2,0,0])
x0 = rotateVector(x0,phaseAngle,latitude)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plotVector(x0, i, ax)
plotVector(x0, j, ax)
plotVector(x0, k, ax)

ax.plot([0],[0],[0], 'ok',markersize=20)
ax.plot([-3],[0],[0], 'oy',markersize=20)
ax.set_xlim((-3,3))
ax.set_ylim((-3,3))
ax.set_zlim((-3,3))


plt.show()
'''