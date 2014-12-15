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
