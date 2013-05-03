#!/usr/bin/env python

import numpy as np
import math
import logging
import correlator
import matplotlib.pyplot as plt
import pylab
from scipy.optimize import leastsq

def fit(cor, tmin, tmax):
    fun = lambda v, x, y: ((v[1] * np.exp((-1.0)*v[0]*x)) -y)
    fn = lambda v, x: ((v[1] * np.exp((-1.0)*v[0]*x)) )
    initial_guess = [2.0, 2.0]
    x = np.arange(32)[tmin:tmax]
    y = np.exp((-1.0)*np.arange(32))
    y = [cor[0][t] for t in range(tmin,tmax)]
    # y = cor[0][tmin:tmax]
    print y
    v, success = leastsq(fun,initial_guess, args=(x,y) )
    if not success:
        raise ValueError()
    print 'Estimater parameters: ', v
    # print 'Real parameters: ', v_real
    X = np.linspace(tmin, tmax ,200*5)
    print cor.times
    # pylab.plot(cor.times,cor[0].values(),'ro', X, fn(v,X))
    # pylab.show()
    emass =  cor.effective_mass(1)
    pylab.plot(emass.keys(),emass.values(),'ro', X, v[0]*np.ones_like(X))
    pylab.show()

def make_fake_cor():
    cfgs = list(range(1))
    times = list(range(32))
    data = {}
    vev = {}
    for c in cfgs:
        vev[c] = 0.0
        tmp = {}
        for t in times:
            tmp[t] = 5.0*np.exp((-2.5)*t)+1.0*np.exp((-5.5)*t)
        data[c] = tmp
    return correlator.Correlator.fromDataDicts(data, vev, vev)

if __name__ == "__main__":
    print(make_fake_cor())
    # pylab.plot(np.arange(32), make_fake_cor()[0])
    # pylab.show()
    cor = make_fake_cor()
    fit(cor, 2, 20)
