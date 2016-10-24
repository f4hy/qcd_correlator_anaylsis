""" Read input of eigenvalues from file where the values are set up as
time step \t lambda0 \t lambda1 \t lambda2 ...
"""
import os
import numpy as np
import logging
import math
import cPickle as pickle

#DEBUG = True
DEBUG = False


def readfile_neigenvalues(basedir, N):
    """reads datafile, returns data object
    of the form
    data[configuration file][time step][ith eigen value]
    """
    data = {}

    assert N > 0, "must give at least one eigenvalue"

    files = list(os.listdir(basedir))
    files.sort()

    if DEBUG:
        files = files[:10]

    for filename in [x for x in files if x[-3:] == "dat"]:
        f = open(basedir + filename)

        rawdata = np.loadtxt(f, comments="#", usecols=(tuple(range(1, N + 1))))
        data[filename] = dict(enumerate(rawdata))

        logging.debug("read %s", filename)

    logging.info("done reading files")
    return data


def readfullfile_config_time(filename):
    """Read-in chimera formated corrilator files
    These are formated so the first line is
    #configs #numbertimes ?unknown ?unknown ?unknown
    then
    timestep  corrilator  imaginarypart
    """

    logging.info("reading chimera formated file %s\n", filename)

    f = open(filename)

    names = ["time", "real", "imag"]

    header = f.readline()       # e.g. 100 40 1 0 1
    numberofconfigs, numberoftimes, _, _, _ = map(int, header.split())

    f.seek(0)                   # Go back to the start

    rawdata = np.genfromtxt(f, comments="#", skiprows=1, dtype=None, names=names)

    data = {}
    for confnumber, config in enumerate(rawdata.reshape(numberofconfigs, numberoftimes)):
        data[confnumber] = dict({time: complex(real, imag) for time, real, imag in config})

    return data


def read_configdir_time_evalues(dirname, N=-1):
    """Read in eigenvalues from a directory where the directory
    represents a config.

    Each file represents a time via filename.time.  Each eigenvalue is
    on a new line.

    E.G
    dirname/1000/smeared_quark_filed_1000_time.0
    """
    assert N > 0, "must give at least one eigenvalue"

    data = {}
    configdata = {}

    logging.debug("Reading from %s", (os.path.relpath(dirname)))

    # Hopfully no files that are digits!
    configs = [x for x in os.listdir(dirname) if x.isdigit()]
    logging.debug("found configs %s", str(configs))

    for config in configs:

        files = list(os.listdir(dirname + config))

        filebasename = files[0].split('.')[0]

        # if DEBUG:
        #     files = files[:128]

        # This is a stupid hack to sort the file names properly, also
        # ensures all of the files are actually there in order

        times = range(len(files))

        configdata = {}
        for time in times:

            filepath = dirname + config + "/" + filebasename + '.' + str(time)
            f = open(filepath)

            rawdata = np.loadtxt(f, comments="#")
            configdata[time] = rawdata[0:N]
            logging.info("Done reading file %s", filepath)

        data[config] = configdata
    return data


def read_configdir_timeorlevel_evalues(dirname, N=-1, store=True, recall=True):
    """Read in eigenvalues from a directory where the directory
    represents a config.

    EITHER
    Each file represents a time via filename_time.#  Each eigenvalue is
    on a new line.

    E.G
    dirname/1000/smeared_quark_filed_1000_time.0

    OR
    Each file represents a eigenvalue over all times via filename_level.#  Each time is
    on a new line.

    E.G
    dirname/1000/smeared_quark_filed_1000_level.0
    """
    assert N > 0, "must give at least one eigenvalue"

    if(recall):
        try:
            logging.info("opening pickle")
            data = pickle.load(open(dirname + "save.p", "rb"))
            return data
        except IOError:
            logging.error("no picked file, building from scratch")

    data = {}
    configdata = {}

    logging.debug("Reading from %s", (os.path.relpath(dirname)))
    # Hopfully no files that are digits!
    configs = [x for x in os.listdir(dirname) if x.isdigit()]

    logging.info(configs)
    # configs = ['1000', '1220']

    for config in configs:
        files = list(os.listdir(dirname + config))

        filebasename = files[0].split('.')[0]

        # if DEBUG:
        #     files = files[:128]

        # This is a stupid hack to sort the file names properly, also
        # ensures all of the files are actually there in order

        if("time" in filebasename):
            logging.debug("Config %s is in time format", config)
            times = range(len(files))
            configdata = {}
            for time in times:
                filepath = dirname + config + "/" + filebasename + '.' + str(time)
                f = open(filepath)

                rawdata = np.loadtxt(f, comments="#")
                configdata[time] = rawdata[0:N]

        elif("level" in filebasename):
            logging.debug("Config %s is in level format", config)
            levels = range(len(files))

            rawdata = None
            for level in levels[0:N]:
                filepath = dirname + config + "/" + filebasename + '.' + str(level)
                logging.debug("reading file %s", filepath)
                f = open(filepath)

                readdata = np.loadtxt(f, comments="#")
                if rawdata is not None:
                    rawdata = np.column_stack((rawdata, readdata))
                else:
                    rawdata = readdata

            configdata = dict(enumerate(rawdata))

        else:
            raise TypeError("unsure which format, exiting")

        data[config] = configdata
    if(store):
        pickle.dump(data, open(dirname + "save.p", "wb") , protocol=pickle.HIGHEST_PROTOCOL)
        logging.info("stored data as %s save.p", dirname)
    return data


def reduce_to_first_eig_val(data):
    newdata = {}
    for config, values in data.iteritems():
        newconfvalue = {}
        for time, array in values.iteritems():
            newconfvalue[time] = array[0]
        newdata[config] = newconfvalue

    return newdata


def reduce_to_trace(data):
    newdata = {}
    for config, values in data.iteritems():
        newconfvalue = {}
        for time, array in values.iteritems():
            newconfvalue[time] = sum(array)
            newdata[config] = newconfvalue

    return newdata


def reduce_to_weighted_trace(data):
    newdata = {}
    for config, values in data.iteritems():
        newconfvalue = {}
        for time, array in values.iteritems():
            array = [lamb*math.exp(-64*lamb*lamb) for lamb in array]
            newconfvalue[time] = sum(array)
            newdata[config] = newconfvalue

    return newdata


def reduce_to_trace_n_n(data, n, m):
    newdata = {}
    for config, values in data.iteritems():
        newconfvalue = {}
        for time, array in values.iteritems():
            newconfvalue[time] = np.sum(array[n:m])
            newdata[config] = newconfvalue

    return newdata


def reduce_to_exponentiate_trace(data):
    newdata = {}
    for config, values in data.iteritems():
        newconfvalue = {}
        for time, array in values.iteritems():
            #newconfvalue[time] = np.sum(np.exp(-1.0*array[5:20]))
            #newconfvalue[time] = np.average(array, weights=np.exp(-2*array))
            newconfvalue[time] = np.average(array, weights=np.arange(len(array)))
            newdata[config] = newconfvalue

    return newdata


def reduce_to_nth_eig_val(data, n):
    newdata = {}
    for config, values in data.iteritems():
        newconfvalue = {}
        for time, array in values.iteritems():
            newconfvalue[time] = array[n]
        newdata[config] = newconfvalue

    return newdata


def diagonalize(data):
    newdata = {}
    for config, values in data.iteritems():
        newconfvalue = {}
        for time, array in values.iteritems():
            newconfvalue[time] = np.outer(array, array)
        newdata[config] = newconfvalue

    return newdata

if __name__ == "__main__":
    logging.debug("no longer runable")
#     testdir = "./alleigens860/"
#     data = readfile_neigenvalues(testdir, 1)
#     import tools
#     x, trash = tools.computecors(data, 1)
#     #print(x)
#     raise SystemExit
#     testarray = np.ones(1)
#     time = 71
#     t_o = 7
#     for config, values in data.iteritems():
#         testarray = {}
#         testarray[t_o] = values[t_o]
#         testarray[time] = values[time]
# #         for time, array in values.iteritems():
# #             testarray = array
# # #            break
#         break
# #    testarray = testarray[:2]
#     #print(testarray[t_o])
#     #print(testarray[time])
#     M= np.mat(np.outer(testarray[t_o], testarray[time]))
#     #print("M")
#     lamb, evec = np.linalg.eig(M)
#     # print "evec", len(evec)
#     # print evec

#     #print "eigenvalues", lamb
#     # print evec[0]
#     # print np.vdot(evec[0], evec[0])
#     # e1= evec[0]/math.sqrt(np.vdot(evec[0], evec[0]))
#     # print e1
#     # print np.vdot(e1, e1)
#     # print np.linalg.inv(evec) * (M * evec)
