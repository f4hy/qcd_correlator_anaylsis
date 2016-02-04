import math
import numpy as np

def jackknife_errors(data, datajk):
    """Compute the error bars from the data and the data with specific
    configs removed"""
    errors = {}
    for shift in data:

        total = 0.0
        M = 0
        for _, val in datajk[shift].iteritems():
            total += ((data[shift]) - (val)) ** 2
            M += 1
        fm = float(M)
        errors[shift] = math.sqrt(((fm - 1) / (fm)) * total)
    return errors


def errorbars(data, datajk):
    total = math.fsum((np.array(datajk.values()) - data)**2)
    fm = float(len(datajk))
    return math.sqrt(((fm - 1.0) / (fm)) * total)
