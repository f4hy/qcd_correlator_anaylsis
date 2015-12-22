#!/usr/bin/env python
import math
import logging

def ncols(N):
    logging.info("Determing columns for {} plots".format(N))
    if N < 6:
        return N
    else:
        if N % 4 == 0:
            return 4
        if N % 3 == 0:
            return 3
    return int(math.sqrt(N))+1

def auto_fit_range(minval, maxval):
    spread = maxval - minval
    return minval-spread*0.1, maxval+spread*0.1

def print_paren_error(value,error):

    try:
        digits = math.floor(math.log10(error))
    except ValueError:
        return value
    if digits >= 1.0:
        formated_string = "{m:d}({e:d})".format(m=int(value), e=int(error))
        return formated_string
    if digits == 0.0:
        formated_string = "{m:0.1f}({e:.1f})".format(m=value, e=error)
        return formated_string

    digits= -1.0*digits
    formated_error = int(round(error * (10**(digits + 1))))
    formated_value = "{m:.{d}f}".format(d=int(digits) + 1, m=value)

    formated_string = "{m}({e})".format(m=formated_value, e=formated_error)
    return formated_string
