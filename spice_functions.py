import numpy as np
import spice


def get_coordinates(UtcStartTime, KernelMetaFile, Target, RefFrame, Abcorr, Observer, UtcStopTime=None, nDeltaT=100):

    spice.furnsh(KernelMetaFile)

    date = datetime.datetime.strptime(UtcStartTime, '%Y-%m-%dT%H:%M:%S')
    x, y, z, r = [], [], [], []

    if UtcStopTime:
        EndDate = datetime.datetime.strptime(UtcStopTime, '%Y-%m-%dT%H:%M:%S')
    else:
        EndDate = datetime.datetime.strptime(UtcStartTime, '%Y-%m-%dT%H:%M:%S')

    dt = datetime.timedelta(seconds=nDeltaT)

    while date < EndDate:
        UtcString = date.strftime("%Y-%m-%dT%H:%M:%S")
        et = spice.str2et(UtcString)
        #rObserver, LightTime = spice.spkpos("ROSETTA", et, "J2000", "NONE", "CHURYUMOV-GERASIMENKO")
        rObserver, LightTime = spice.spkpos(Target, et, Frame, Abcorr, Observer)

        x.append(rTarget[0])
        y.append(rTarget[1])
        z.append(rTarget[2])
        r.append(np.sqrt(np.sum(np.array(rTarget)**2)))
        date += dt

    return x, y, z, r
