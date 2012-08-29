#!/usr/bin/env python
""" Code to read input in format

#update number 5
t1,    data
t2,    data
...
#update number 10
t1,    data
t2,    data

"""

import numpy as np
import configtimeobj
import correlator


def read_config_time_data_real_imag(filename, configs=None, times=None):
    """Takes a file ofr the form in the header of this source file and
    returns a configtimeobject

    """
    f = open(filename)

    print "reading data from %s" % filename

    rawdata = np.genfromtxt(f, delimiter=",", comments="#", autostrip=True, dtype='int,float,float')
    f.close()
    if  configs and times:
        print "using configs: %d and times: %d" % (configs, times)
    else:
        print "dimensions not set will guess from data"
        times, configs = guess_dimensions(rawdata)

    data = rawdata.reshape(configs, times)
    #data = data[:10]
    #print data
    return configtimeobj.Cfgtimeobj.fromListTuple(data)


def read_config_time_data_real(filename, configs=None, times=None):
    """Takes a file ofr the form in the header of this source file and
    returns a configtimeobject

    """
    f = open(filename)

    print "reading data from %s" % filename

    rawdata = np.genfromtxt(f, delimiter=",", comments="#",
                            autostrip=True, dtype='int,float', usecols=(0, 1))
    f.close()

    if  configs and times:
        print "using configs: %d and times: %d" % (configs, times)
    else:
        print "dimensions not set will guess from data"
        times, configs = guess_dimensions(rawdata)

    data = rawdata.reshape(configs, times)
    #data = data[:10]
    #print data
    return configtimeobj.Cfgtimeobj.fromListTuple(data)


def read_correlator(filename, configs=None, times=None):
    """Takes a file ofr the form in the header of this source file and
    returns a configtimeobject

    """
    f = open(filename)

    print "reading data from %s" % filename

    rawdata = np.genfromtxt(f, delimiter=",", comments="#",
                            autostrip=True, dtype='int,float', usecols=(0, 1))
    f.close()

    if  configs and times:
        print "using configs: %d and times: %d" % (configs, times)
    else:
        print "dimensions not set will guess from data"
        times, configs = guess_dimensions(rawdata)

    data = rawdata.reshape(configs, times)
    #data = data[:10]
    #print data
    return correlator.Correlator.fromListTuple(data)


def read_config_vev(filename, configs=None):
    f = open(filename)

    print "reading vev from %s" % filename

    rawdata = np.genfromtxt(f, delimiter=",", comments="#",
                            autostrip=True, dtype='float', usecols=0)
    f.close()

    #print rawdata

    if  configs:
        print "using configs: %d" % (configs)
    else:
        print "vev dimensions not set will guess from data"
        configs = len(rawdata)

    # rawdata = rawdata[:10]
    # configs = 10
    vevdict = {}
    for config in range(configs):
        vevdict[config] = rawdata[config]

    return vevdict


def guess_dimensions(rawdata):
    print rawdata[0][0]
    if(rawdata[0][0] != 0):
        raise Exception("does not start on time 0 can not guess")

    time = 0
    for datum in rawdata:
        if(datum[0] == time):
            time += 1
        elif(datum[0] == 0):
            configs = ((rawdata.shape[0]) / time)
            print "guessing time = %d" % time
            print "guessing configs = %d" % configs
            return (time, configs)
        else:
            raise Exception("number of times could not be guessed")
