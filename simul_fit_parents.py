import numpy as np
from iminuit import Minuit

from fit_parents import mass_bounds, amp_bounds, const_bounds

class sharedmass_amp(object):
    """Parent class for functions which take a mass and an amplitude"""
    def __init__(self):
        self.starting_guess = self.thisguess
        self.bounds = [mass_bounds, amp_bounds, amp_bounds]
        self.parameter_names = ["mass", "amp1", "amp2"]
        self.subtract = False
        self.stride = 1

    def thisguess(self, cor, period, *args):
        dt = 1
        ave = cor.average_sub_vev()
        emass = cor.periodic_effective_mass(dt, fast=True, period=period)
        mass_guess = np.mean([emass[i[1]-dt-1] for i in self.indexes])

        mid1 = (self.indexes[0][0]+self.indexes[0][1])/2
        mid2 = (self.indexes[1][0]+self.indexes[1][1])/2
        rmid1 = (self.ranges[0][0]+self.ranges[0][1])/2
        rmid2 = (self.ranges[1][0]+self.ranges[1][1])/2
        amp_guess1 = ave[mid1]*np.exp(mass_guess*(rmid1))
        amp_guess2 = ave[mid2]*np.exp(mass_guess*(rmid2))
        return [mass_guess, amp_guess1, amp_guess2]


    def my_cov_fun(self, mass, amp1, amp2):
        vect = self.aoc - self.formula((mass, amp1, amp2), self.times)
        return vect.dot(self.inv_cov).dot(vect)

    def valid(self, *kargs):
        return True

    def custom_minuit(self, data, invmatrix, times, guess):
        self.aoc = data
        self.inv_cov = invmatrix
        self.times = times
        dof = len(guess)+len(data)
        m = Minuit(self.my_cov_fun, mass=guess[0], error_mass=guess[0]*0.1,
                   amp1=guess[1], error_amp1=guess[1]*0.1,
                   amp2=guess[2], error_amp2=guess[2]*0.1,
                   errordef=dof,
                   print_level=0, pedantic=True)
        return m

class shared_twice_mass_amp(object):
    """Parent class for functions which take a mass and an amplitude"""
    def __init__(self):
        self.starting_guess = self.thisguess
        self.bounds = [mass_bounds, amp_bounds, amp_bounds, mass_bounds, amp_bounds, amp_bounds]
        self.parameter_names = ["massa", "amp1a", "amp2a", "massb", "amp1b", "amp2b"]
        self.subtract = False
        self.stride = 1

    def thisguess(self, cor, period, *args):
        dt = 1
        ave = cor.average_sub_vev()
        emass = cor.periodic_effective_mass(dt, fast=False, period=period)
        massa_guess = np.mean([emass[i[1]-dt-1] for i in self.indexes])
        massb_guess = massa_guess

        amp_guess1a = ave[self.indexes[0][0]]*np.exp(massa_guess*(self.ranges[0][0]))
        amp_guess2a = ave[self.indexes[1][0]]*np.exp(massa_guess*(self.ranges[1][0]))
        amp_guess1b = ave[self.indexes[0][0]]*np.exp(massb_guess*(self.ranges[0][0]))
        amp_guess2b = ave[self.indexes[1][0]]*np.exp(massb_guess*(self.ranges[1][0]))
        return [massa_guess, amp_guess1a, amp_guess2a, massb_guess, amp_guess1b, amp_guess2b]


    def my_cov_fun(self, massa, amp1a, amp2a, massb, amp1b, amp2b):
        vect = self.aoc - self.formula((massa, amp1a, amp2a, massb, amp1b, amp2b), self.times)
        return vect.dot(self.inv_cov).dot(vect)

    def valid(self, *kargs):
        return True

    def custom_minuit(self, data, invmatrix, times, guess):
        self.aoc = data
        self.inv_cov = invmatrix
        self.times = times
        dof = len(guess)+len(data)
        m = Minuit(self.my_cov_fun, massa=guess[0], error_massa=guess[0]*0.1,
                   amp1a=guess[1], error_amp1a=guess[1]*0.1,
                   amp2a=guess[2], error_amp2a=guess[2]*0.1,
                   massb=guess[0], error_massb=guess[0]*0.1,
                   amp1b=guess[1], error_amp1b=guess[1]*0.1,
                   amp2b=guess[2], error_amp2b=guess[2]*0.1,
                   errordef=dof,
                   print_level=0, pedantic=True)
        return m


class multirange(object):
    """ Parent class for functions which are periodic and need to know the time extent"""
    def setranges(self, ranges):
        self.ranges = ranges
        indexes = []
        prev = 0
        for i in ranges:
            length = (i[1]-i[0])
            indexes.append((prev, prev+length))
            prev = prev+length+1
        self.indexes = indexes
