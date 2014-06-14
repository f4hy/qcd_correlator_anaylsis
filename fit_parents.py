import numpy as np
from iminuit import Minuit
import logging

mass_bounds = (0.005, 5.0)
amp_bounds = (0.0, 1.0e5)
const_bounds = (-5.0, 1.0e8)

class InvalidFit(RuntimeError):
    pass

def massamp_guess(cor, tmax, *args):
    dt = 3
    maxt = tmax - dt
    ave = cor.average_sub_vev()
    emass = cor.effective_mass(dt)
    if not emass[maxt] > 0:
        for t in range(maxt,0,-1):
            if emass[t] > 0:
                mass_guess = emass[t]
                break

    else:
        mass_guess = emass[maxt]
    amp_guess = ave[maxt]*np.exp(mass_guess*(maxt))
    return [mass_guess, amp_guess]


def const_guess(cor, tmax, *args):
    dt = 3
    maxt = tmax - dt
    ave = cor.average_sub_vev()
    emass = cor.effective_mass(dt)
    if not emass[maxt] > 0:
        for t in range(maxt,0,-1):
            if emass[t] > 0:
                mass_guess = emass[t]
                break

    else:
        mass_guess = emass[maxt]
    amp_guess = ave[maxt]*np.exp(mass_guess*(maxt))
    return [mass_guess, amp_guess, 0.01]


def twoexp_sqr_guess(cor, tmax, tmin):
    dt = 3
    maxt = tmax - dt
    ave = cor.average_sub_vev()
    emass = cor.effective_mass(dt)
    mass_guess = np.median(emass.values())
    amp_guess = ave[maxt]*np.exp(mass_guess*(maxt))
    mass2_guess = np.sqrt(emass[tmin])
    amp2_guess = (abs(ave[tmin] - amp_guess*np.exp(-mass_guess*tmin)) /
                  (amp_guess*np.exp(-(mass_guess+mass2_guess**2)*tmin)))/2
    return [mass_guess, amp_guess, mass2_guess, amp2_guess]

def twoexp_sqr_const_guess(cor, tmax, tmin):
    return twoexp_sqr_guess(cor, tmax, tmin) + [0.001]


class periodic(object):
    """ Parent class for functions which are periodic and need to know the time extent"""
    def setNt(self, Nt):
        self.Nt = Nt
        while not self.Nt:
            try:
                self.Nt = int(raw_input('Time period not specified, please enter Nt:'))
            except ValueError:
                print "Not a valid number"


class mass_amp(object):
    """Parent class for functions which take a mass and an amplitude"""
    def __init__(self):
        self.starting_guess = massamp_guess
        self.bounds = [mass_bounds, amp_bounds]
        self.parameter_names = ["mass", "amp"]
        self.subtract = False

    def my_cov_fun(self, mass, amp):
        vect = self.aoc - self.formula((mass, amp), self.times)
        return vect.dot(self.inv_cov).dot(vect)

    def valid(self, *kargs):
        return True

    def custom_minuit(self, data, invmatrix, times, guess):
        self.aoc = data
        self.inv_cov = invmatrix
        self.times = times
        m = Minuit(self.my_cov_fun, mass=guess[0], amp=guess[1],
                   print_level=0, pedantic=False)
        return m


class mass_amp_subtracted(object):
    """Parent class for functions which take a mass and an amplitude"""
    def __init__(self):
        self.starting_guess = subtracted_guess
        self.bounds = [mass_bounds, amp_bounds]
        self.parameter_names = ["mass", "amp"]
        self.subtract = True

    def my_cov_fun(self, mass, amp):
        vect = self.aoc - self.formula((mass, amp), self.times)
        return vect.dot(self.inv_cov).dot(vect)

    def valid(self, *kargs):
        return True

    def custom_minuit(self, data, invmatrix, times, guess):
        self.aoc = data
        self.inv_cov = invmatrix
        self.times = times
        m = Minuit(self.my_cov_fun, mass=guess[0], amp=guess[1],
                   print_level=0, pedantic=False)
        return m


class mass_amp_const(object):
    """Parent class for functions which take a mass and an amplitude"""
    def __init__(self):
        self.starting_guess = const_guess
        self.bounds = [mass_bounds, amp_bounds, const_bounds]
        self.parameter_names = ["mass", "amp", "const"]
        self.subtract = False

    def my_cov_fun(self, mass, amp, const):
        vect = self.aoc - self.formula((mass, amp, const), self.times)
        return vect.dot(self.inv_cov).dot(vect)

    def valid(self, *kargs):
        return True

    def custom_minuit(self, data, invmatrix, times, guess):
        self.aoc = data
        self.inv_cov = invmatrix
        self.times = times
        m = Minuit(self.my_cov_fun, mass=guess[0], amp=guess[1], const=guess[2],
                   print_level=0, pedantic=False)
        return m


class twice_mass_amp(object):
    """Parent class for functions which take a mass and an amplitude"""
    def __init__(self):
        self.starting_guess = twoexp_sqr_guess
        self.bounds = [mass_bounds, amp_bounds, mass_bounds, amp_bounds]
        self.parameter_names = ["mass", "amp", "mass2", "amp2"]
        self.subtract = False

    def my_cov_fun(self, mass, amp, mass2, amp2):
        vect = self.aoc - self.formula((mass, amp, mass2, amp2), self.times)
        return vect.dot(self.inv_cov).dot(vect)

    def valid(self, params):
        mass, amp, mass2, amp2 = params
        return True
        if amp2 > 10*amp:
            logging.error("Invalid fit with paramters {}".format(repr(params)))
            # raise InvalidFit
            return False
        else:
            return True

    def custom_minuit(self, data, invmatrix, times, guess):
        self.aoc = data
        self.inv_cov = invmatrix
        self.times = times
        m = Minuit(self.my_cov_fun, mass=guess[0], amp=guess[1], mass2=guess[2], amp2=guess[3],
                   print_level=0, pedantic=False, limit_amp2=amp_bounds, limit_mass2=mass_bounds,
                   limit_mass=mass_bounds, limit_amp=amp_bounds)
        return m


class twice_mass_amp_const(object):
    """Parent class for functions which take a mass and an amplitude"""
    def __init__(self):
        self.starting_guess = twoexp_sqr_const_guess
        self.bounds = [mass_bounds, amp_bounds, mass_bounds, amp_bounds, const_bounds]
        self.parameter_names = ["mass", "amp", "mass2", "amp2", "const"]
        self.subtract = False

    def my_cov_fun(self, mass, amp, mass2, amp2, const):
        vect = self.aoc - self.formula((mass, amp, mass2, amp2, const), self.times)
        return vect.dot(self.inv_cov).dot(vect)

    def valid(self, params):
        mass, amp, mass2, amp2, const = params
        if amp2 > 10*amp:
            logging.error("Invalid fit with paramters {}".format(repr(params)))
            raise InvalidFit
            return False
        else:
            return True

    def custom_minuit(self, data, invmatrix, times, guess):
        self.aoc = data
        self.inv_cov = invmatrix
        self.times = times
        m = Minuit(self.my_cov_fun, mass=guess[0], amp=guess[1], mass2=guess[2], amp2=guess[3], const=guess[4],
                   print_level=0, pedantic=False, limit_amp2=amp_bounds, limit_mass2=mass_bounds,
                   limit_mass=mass_bounds, limit_amp=amp_bounds)
        return m
