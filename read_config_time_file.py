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
import logging


def read_config_time_data_real_imag(filename, configs=None, times=None):
    """Takes a file ofr the form in the header of this source file and
    returns a configtimeobject

    """
    f = open(filename)

    logging.info("reading data from %s", filename)

    rawdata = np.genfromtxt(f, delimiter=",", comments="#", autostrip=True, dtype='int,float,float')
    f.close()
    if  configs and times:
        logging.info("using configs: %d and times: %d", configs, times)
    else:
        logging.debug("dimensions not set will guess from data")
        times, configs = guess_dimensions(rawdata)

    data = rawdata.reshape(configs, times)

    return configtimeobj.Cfgtimeobj.fromListTuple(data)


def read_config_time_data_real(filename, configs=None, times=None):
    """Takes a file ofr the form in the header of this source file and
    returns a configtimeobject

    """
    return configtimeobj.Cfgtimeobj.fromListTuple(
        read_config_time_data_real_dict(filename, configs, times))


def read_config_time_data_real_dict(filename, configs=None, times=None):
    """Takes a file ofr the form in the header of this source file and
    returns a dict

    """
    f = open(filename)

    logging.info("reading data from %s", filename)

    rawdata = np.genfromtxt(f, delimiter=",", comments="#",
                            autostrip=True, dtype='int,float', usecols=(0, 1))
    f.close()

    if  configs and times:
        logging.info("using configs: %d and times: %d", configs, times)
    else:
        logging.debug("dimensions not set will guess from data")
        times, configs = guess_dimensions(rawdata)

    data = rawdata.reshape(configs, times)
    #data = data[:10]
    return data


def read_correlator(filename, configs=None, times=None):
    """Takes a file ofr the form in the header of this source file and
    returns a configtimeobject

    """
    f = open(filename)

    logging.info("reading data from %s", filename)

    rawdata = np.genfromtxt(f, delimiter=",", comments="#",
                            autostrip=True, dtype='int,float', usecols=(0, 1))
    f.close()

    if  configs and times:
        logging.info("using configs: %d and times: %d", configs, times)
    else:
        logging.debug("dimensions not set will guess from data")
        times, configs = guess_dimensions(rawdata)

    data = rawdata.reshape(configs, times)
    #data = data[:10]
    return correlator.Correlator.fromListTuple(data)


def read_config_vev(filename, configs=None):
    f = open(filename)

    logging.info("reading vev from %s", filename)

    rawdata = np.genfromtxt(f, delimiter=",", comments="#",
                            autostrip=True, dtype='float', usecols=0)
    f.close()

    if  configs:
        logging.info("using configs: %d", configs)
    else:
        logging.debug("vev dimensions not set will guess from data")
        configs = len(rawdata)

    # rawdata = rawdata[:10]
    # configs = 10
    vevdict = {}
    for config in range(configs):
        vevdict[config] = rawdata[config]

    return vevdict


def guess_dimensions(rawdata):

    if(rawdata[0][0] != 0):
        raise Exception("does not start on time 0 can not guess")

    time = 0
    for datum in rawdata:
        if(datum[0] == time):
            time += 1
        elif(datum[0] == 0):
            configs = ((rawdata.shape[0]) / time)
            logging.debug("guessing time = %d \tguessing configs = %d", time, configs)
            return (time, configs)
        else:
            raise Exception("number of times could not be guessed")
