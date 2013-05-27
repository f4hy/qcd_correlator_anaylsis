#!/usr/bin/env python
import logging
import pandas as pd
import re
import numpy as np

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


def read_datadict_paraenformat_real(filename, real=True):
    f = open("corsnk-etap000DDL7Egp1_src-etap000DDL7Egp1.dat")
    df = pd.read_csv(f, delimiter=' ', names=["time", "correlator"],
                     converters={1: parse_pair})

    if(real):
        df = df.apply(np.real)  # Take the real part

    vc = df.time.value_counts()
    times = vc.index
    time_counts = vc.values
    df["config"] = pd.Series(map(lambda x: x/len(times), df.index))
    vcc = df.config.value_counts()
    cfgs = vcc.index
    cfgs_counts = vcc.values

    df = df.set_index("config")
    df = df.set_index("time", append=True)

    logging.debug("Read file, got {} configs and {} times".format(cfgs_counts[0], time_counts[0]))

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


if __name__ == "__main__":
    filename = "corsnk-etap000DDL7Egp1_src-etap000DDL7Egp1.dat"
    print read_datadict_paraenformat_real(filename)
