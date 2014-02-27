import numpy as np
import iminuit
from fit_parents import periodic, mass_amp, mass_amp_const, twice_mass_amp

mass_bounds = (0.005, 2.0)
amp_bounds = (0.0, 1.0e8)
const_bounds = (-5.0, 1.0e8)




class cosh(mass_amp, periodic):
    def __init__(self, Nt=None):
        super(cosh, self).__init__()
        self.setNt(Nt)
        self.description = "cosh"
        self.template = "{1: f}Cosh(-{0: f}*(t-%d/2))" % self.Nt

    def formula(self, v, x):
        #return ((2*v[1])/np.exp(v[0]*Nt/2.0) * np.cosh((-1.0)* v[0]*((x-(Nt/2.0)))))
        return (v[1] * np.cosh((-1.0)*v[0]*((x-(self.Nt/2.0)))))


class single_exp(mass_amp):
    def __init__(self, **kargs):
        super(single_exp, self).__init__()
        self.description = "exp"
        self.template = "{1: f}exp(-{0: f}*t)"

    def formula(self, v, x):
        return (v[1] * np.exp((-1.0) * v[0] * x))

class periodic_exp(mass_amp, periodic):
    def __init__(self, Nt=None):
        super(periodic_exp, self).__init__()
        self.setNt(Nt)
        self.description = "fwd-back-exp"
        self.template = "{1: f}(exp(-{0: f}*t)+exp(-{0: f}*(t-%d))" % self.Nt

    def formula(self, v, x):
        return (v[1] * (np.exp((-1.0) * v[0] * x) + np.exp(v[0] * (x-(self.Nt)))))

class periodic_exp_const(mass_amp_const, periodic):
    def __init__(self, Nt=None):
        super(periodic_exp_const, self).__init__()
        self.setNt(Nt)
        self.description = "fwd-back-exp_const"
        self.template = "{1: f}(exp(-{0: f}*t)+exp(-{0: f}*(t-%d))+{2: f}" % self.Nt

    def formula(self, v, x):
        return (v[1] * (np.exp((-1.0) * v[0] * x) + np.exp(v[0] * (x-(self.Nt)))))+v[2]

class cosh_const(mass_amp_const, periodic):
    def __init__(self, Nt=None):
        super(cosh_const, self).__init__()
        self.setNt(Nt)
        self.description = "cosh+const"
        self.template = "{1: f}Cosh(-{0: f}*(t-%d/2))+{2: f}" % self.Nt

    def formula(self, v, x):
        return (v[1] * np.cosh((-1.0)*v[0]*((x-(self.Nt/2.0)))))+v[2]


class two_exp(twice_mass_amp):
    def __init__(self, **kargs):
        super(two_exp, self).__init__()
        self.description = "two_exp"
        self.template = "{1: f}exp(-{0: f}*t)(1+{3: f}exp(-{2: f}^2*t)"

    def formula(self, v, x):
        return (v[1] * np.exp((-1.0) * v[0] * x)*(1.0 + v[3]*np.exp((-1.0)*(v[2]**2)*x)))


class periodic_two_exp(twice_mass_amp):
    def __init__(self, Nt=None):
        super(periodic_two_exp, self).__init__()
        self.setNt(Nt)
        self.description = "periodic_two_exp"
        self.template = "{1: f}exp(-{0: f}*t)(1+{3: f}exp(-{2: f}^2*t)"
    def formula(self, v, x):
                return ((v[1]*np.exp((-1.0)*v[0]*x)*(1.0 + v[3]*np.exp((-1.0)*(v[2]**2)*x))) +
                        (v[1]*np.exp(v[0]*(x-(self.Nt)))*(1.0 + v[3]*np.exp((v[2]**2)*(x-(self.Nt))))))


class jlab:
    def __init__(self, **kargs):
        self.starting_guess = twoexp_sqr_guess
        self.parameter_names = ["mass", "amp", "mass2"]
        self.description = "jlab"
        self.template = "{1: f}exp(-{0: f}*t)+{1: f}exp(-{2: f}*t)"

    def formula(self, v, x):
        return ((1 - v[1]) * np.exp((-1.0) * v[0] * x))+(v[1] * np.exp((-1.0) * v[2] * x))



def pade_guess(**kargs):
    return [0.05, 100, 1.0, 10]

class pade:
    """ exp( - E * t )  *    A /  ( 1 +   a1* t + a2 * t^2 ...  ) """
    def __init__(self, Nt=None):
        self.starting_guess = pade_guess
        self.parameter_names = ["mass", "amp", "B", "C"]
        self.description = "Pade"
        self.template = "{1: f}exp(-{0: f}*t)/(1+{3: f}t +{2: f}*t^2)"
        self.Nt = Nt
        # while not self.Nt:
        #     try:
        #         self.Nt = int(raw_input('Time period not specified, please enter Nt:'))
        #     except ValueError:
        #         print "Not a valid number"

    def formula(self, v, x):
        return (v[1]*np.exp((-1.0)*v[0]*x)) / (1.0 + v[2]*x+v[3]*x**2)
