import sys
import numpy as np
import spice
import datetime


def get_coordinates(UtcStartTime, KernelMetaFile, Target, RefFrame, Abcorr, Observer, UtcStopTime=None, nDeltaT=100):

    try:
        spice.furnsh(KernelMetaFile)

        date = datetime.datetime.strptime(UtcStartTime, '%Y-%m-%dT%H:%M:%S')
        x, y, z, r, dd = [], [], [], [], []

        if UtcStopTime:
            EndDate = datetime.datetime.strptime(UtcStopTime, '%Y-%m-%dT%H:%M:%S')
        else:
            EndDate = datetime.datetime.strptime(UtcStartTime, '%Y-%m-%dT%H:%M:%S')

        dt = datetime.timedelta(seconds=nDeltaT)
        while date < EndDate:
            UtcString = date.strftime("%Y-%m-%dT%H:%M:%S")
            et = spice.str2et(UtcString)
            #rObserver, LightTime = spice.spkpos("ROSETTA", et, "J2000", "NONE", "CHURYUMOV-GERASIMENKO")
            rTarget, LightTime = spice.spkpos(Target, et, RefFrame, Abcorr, Observer)

            x.append(rTarget[0])
            y.append(rTarget[1])
            z.append(rTarget[2])
            r.append(np.sqrt(np.sum(np.array(rTarget)**2)))
            dd.append(datetime.datetime.strptime(UtcString, '%Y-%m-%dT%H:%M:%S'))
            date += dt

        # return values in meters, hence * 1000
        return np.array(x) * 1000, np.array(y) * 1000, np.array(z) * 1000, np.array(r) * 1000, dd
    except Exception, e:
        print "#####################################################################"
        print "SPICE ERROR: Loaded kernel does not contain information for the full"
        print "range of dates that were specified by the user. Please make sure you "
        print "selected the proper spice kernel or shorten your selected data range."
        print "To get newest version of the operationalKernel you can  go to"
        print "ices-dev.engin.umich.edu"
        if len(x) >=1:
            print ""
            print " --> Continue run with limited data that is available in spice kernels:"
            print "tStart: %s" %(UtcStartTime)
            print "tStop : %s" %(UtcString)
            return np.array(x) * 1000, np.array(y) * 1000, np.array(z) * 1000, np.array(r) * 1000, dd
            print "#####################################################################"

        else:
            print " - Run unsuccessful, exiting now."
            print "#####################################################################"
            sys.exit()
