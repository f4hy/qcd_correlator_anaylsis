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

    def thisguess(self, cor, *args):
        dt = 1
        ave = cor.average_sub_vev()
        emass = cor.cosh_effective_mass(dt)
        mass_guess = np.mean([emass[i[1]-dt-1] for i in self.indexes])

        amp_guess1 = ave[self.indexes[0][1]]*np.exp(mass_guess*(self.ranges[0][1]))
        amp_guess2 = ave[self.indexes[1][1]]*np.exp(mass_guess*(self.ranges[1][1]))
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
        m = Minuit(self.my_cov_fun, mass=guess[0], amp1=guess[1], amp2=guess[2],
                   print_level=0, pedantic=False)
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
