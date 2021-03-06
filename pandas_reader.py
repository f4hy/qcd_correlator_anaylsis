#!/usr/bin/env python
import logging
import pandas as pd
import re
import numpy as np
from cStringIO import StringIO

import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)


pair = re.compile(r'\(([^,\)]+),([^,\)]+)\)')


def parse_pair(s):
    return complex(*map(float, pair.match(s).groups()))


def complex_converter(txt):
    txt = txt.strip("()") + "j"
    txt = txt.replace(",-", "-")
    txt = txt.replace(",", "+")
    # print txt
    return complex(txt)


def lines_without_comments(filename, comment="#"):
    s = StringIO()
    with open(filename) as f:
        for line in f:
            if not line.startswith(comment):
                s.write(line)
    s.seek(0)
    return s


def read_full_paraenformat_mean(filename):
    f = lines_without_comments(filename)
    df = pd.read_csv(f, delimiter=' ', names=["time", "correlator"],
                     converters={1: parse_pair})
    return df.groupby('time').mean()


def read_normal_paraenformat(filename):
    f = lines_without_comments(filename)
    df = pd.read_csv(f, delimiter=' ', names=["time", "correlator"],
                     converters={1: parse_pair})
    return df

def read_normal(filename):
    try:
        f = lines_without_comments(filename)
        df = pd.read_csv(f, delimiter=' ', names=["time", "correlator"],
                         converters={1: parse_pair})
    except AttributeError:
        logging.info("Failed to read with pandas, reading normal")
        f = lines_without_comments(filename)
        df = pd.read_csv(f, delimiter=',', names=["time", "correlator"])
    return df

def read_configcols_paraenformat(filename):
    f = lines_without_comments(filename)
    df = pd.read_csv(f, delimiter=' ', names=["time", "correlator"],
                     converters={1: parse_pair})
    grouped = df.groupby('time')
    joined = df.groupby('time').apply(
        lambda g: pd.Series(g['correlator'].values)
    ).rename(columns=lambda x: 'correlator%s' % x)
    return joined

def read_configcols_normal(filename):
    f = lines_without_comments(filename)
    df = pd.read_csv(f, delimiter=',', names=["time", "correlator", "correlator_imag"])
    grouped = df.groupby('time')
    joined = df.groupby('time').apply(
        lambda g: pd.Series(g['correlator'].values)
    ).rename(columns=lambda x: 'correlator%s' % x)
    return joined

def read_full_paraenformat_multiindex(filename):
    f = lines_without_comments(filename)
    df = pd.read_csv(f, delimiter=' ', names=["time", "correlator"],
                     converters={1: parse_pair})
    times = set(df.time)
    df["config"] = pd.Series(map(lambda x: x/len(times), df.index))
    multiindexed = df.set_index(["config","time"])
    grouped = multiindexed.groupby(level=1)
    return grouped


def read_single_time_paraenformat(filename, t):
    f = lines_without_comments(filename)
    logging.info("reading time {} from file {}".format(t, filename))
    df = pd.read_csv(f, delimiter=' ', names=["time", "data"],
                     converters={1: parse_pair})
    single_time = df.ix[df["time"] == t, "data"]
    return single_time.values

def read_single_time_commaformat(filename, t):
    f = lines_without_comments(filename)
    logging.info("reading time {} from file {}".format(t, filename))
    df = pd.read_csv(f, delimiter=',', names=["time", "data"] )
    single_time = df.ix[df["time"] == t, "data"]
    return single_time.values


def read_datadict_paraenformat_real(filename, real=True):
    f = lines_without_comments(filename)
    df = pd.read_csv(f, delimiter=' ', names=["time", "correlator"],
                     converters={1: parse_pair})

    if(real):
        df = df.apply(np.real, reduce=False)  # Take the real part

    vc = df.time.value_counts()
    times = vc.index
    time_counts = vc.values
    df["config"] = pd.Series(map(lambda x: x/len(times), df.index))
    vcc = df.config.value_counts()
    cfgs = vcc.index
    cfgs_counts = vcc.values

    df = df.set_index("config")
    df = df.set_index("time", append=True)

    logging.info("Read file, got {} configs and {} times".format(time_counts[0], cfgs_counts[0]))

    logging.debug("checking consistancy")
    assert all(time_counts[0] == count for count in time_counts), "Inconsistant time counts!"
    assert all(cfgs_counts[0] == count for count in cfgs_counts), "Inconsistant cfgs counts!"

    data_dict = {}
    for c in cfgs:
        # data_dict[c] = {}
        tmp_dict = {}
        for t in times:
            tmp_dict[t] = df.ix[(c, t)][0]

        data_dict[c] = tmp_dict

    return data_dict
    # for row in df.iterrows():
    #     print row


def read_datadict_commacomplex(filename, real=True):

    f = lines_without_comments(filename)
    df = pd.read_csv(f, delimiter=',', comment="#", names=["time", "correlator", "correlator_imag"])

    if(real):
        del df["correlator_imag"]
    else:
        df.correlator = df.correlator+(df.correlator_imag*complex(0,1))

    vc = df.time.value_counts()
    times = vc.index
    time_counts = vc.values
    df["config"] = pd.Series(map(lambda x: x/len(times), df.index))
    vcc = df.config.value_counts()
    cfgs = vcc.index
    cfgs_counts = vcc.values

    df = df.set_index("config")
    df = df.set_index("time", append=True)

    logging.info("Read file, got {} configs and {} times".format(time_counts[0], cfgs_counts[0]))

    logging.debug("checking consistancy")
    assert all(time_counts[0] == count for count in time_counts), "Inconsistant time counts!"
    assert all(cfgs_counts[0] == count for count in cfgs_counts), "Inconsistant cfgs counts!"

    data_dict = {}
    for c in cfgs:
        # data_dict[c] = {}
        tmp_dict = {}
        for t in times:
            tmp_dict[t] = df.ix[(c, t)][0]

        data_dict[c] = tmp_dict

    return data_dict

def read_baryon(filename, real=True):

    f = lines_without_comments(filename)
    df = pd.read_csv(f, delimiter=',', comment="#", names=["time", "correlator", "correlator_imag"])


    vc = df.time.value_counts()
    times = vc.index
    time_counts = vc.values
    df["config"] = pd.Series(map(lambda x: x/len(times), df.index))
    vcc = df.config.value_counts()
    cfgs = vcc.index
    cfgs_counts = vcc.values


    df = df.set_index("config")
    df = df.set_index("time", append=True)

    newtimes =  [0] + list(reversed(range(1,len(times))))
    newindex = [(i,j)  for i in sorted(cfgs) for j in newtimes]
    reordered = df["correlator_imag"].reindex(newindex)
    reordered.index = newindex
    timelist = list(range(len(times))) * len(cfgs)
    result = (df["correlator"] + reordered.values)/2.0
    final = pd.DataFrame(data={"time": timelist, "correlator":result.values})
    return final



def read_vev_parenformat(filename, real=True):
    f = lines_without_comments(filename)
    df = pd.read_csv(f, delimiter=' ', names=["garbage", "vev"],
                     converters={1: parse_pair})

    if(real):
        df = df.apply(np.real)  # Take the real part

    cfgs = df.index
    datadict = {}
    for c in cfgs:
        datadict[c] = df.vev[c]
    return datadict


def read_datadict_ambiguouscomplex(filename, real=True):
    if determine_format(filename) == "complex pairs":
        return read_datadict_paraenformat_real(filename, real)
    else:
        return read_datadict_commacomplex(filename, real)


def determine_format(filename):
    f = open(filename)
    for linenum, l in enumerate(f):
        if any([pair.match(x) for x in l.split()]):
            logging.debug("filetype complex pairs")
            return "complex pairs"
        if linenum > 5:
            break
    else:
        logging.debug("filetype comma columns")
        return "comma columns"


if __name__ == "__main__":
    filename = "corsnk-etap000DDL7Egp1_src-etap000DDL7Egp1.dat"
    determine_format(filename)
    # print read_datadict_paraenformat_real(filename)
    filename = "/home/bfahy/r2/pruning/special/atrestpions/correlators_myformat/corsnk-pionp000SD0A1um1_src-pionp000SD0A1um1.dat"
    #determine_delimiter(filename)
    #print read_datadict_commacomplex(filename)
