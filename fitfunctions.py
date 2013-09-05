import numpy as np


class cosh:
    def __init__(self, Nt=None):
        self.starting_guess = [0.5, 50.0]
        self.parameter_names = ["mass", "amp"]
        self.description = "cosh"
        self.Nt = Nt
        while not self.Nt:
            try:
                self.Nt = int(raw_input('Time period not specified, please enter Nt:'))
            except ValueError:
                print "Not a valid number"
        self.template = "{1: f}Cosh(-{0: f}*(t-%d/2))" % self.Nt

    def formula(self, v, x):
        #return ((2*v[1])/np.exp(v[0]*Nt/2.0) * np.cosh((-1.0)* v[0]*((x-(Nt/2.0)))))
        return (v[1] * np.cosh((-1.0)*v[0]*((x-(self.Nt/2.0)))))


class single_exp:
    def __init__(self, **kargs):
        self.starting_guess = [0.5, 10.0]
        self.parameter_names = ["mass", "amp"]
        self.description = "exp"
        self.template = "{1: f}exp(-{0: f}*t)"

    def formula(self, v, x):
        return (v[1] * np.exp((-1.0) * v[0] * x))


class periodic_exp:
    def __init__(self, Nt=None):
        self.starting_guess = [1.0, 1000.0]
        self.parameter_names = ["mass", "amp"]
        self.description = "fwd-back-exp"
        self.Nt = Nt
        while not self.Nt:
            try:
                self.Nt = int(raw_input('Time period not specified, please enter Nt:'))
            except ValueError:
                print "Not a valid number"
        self.template = "{1: f}(exp(-{0: f}*t)+exp(-{0: f}*(t-%d))" % self.Nt

    def formula(self, v, x):
        return (v[1] * (np.exp((-1.0) * v[0] * x) + np.exp(v[0] * (x-(self.Nt)))))


class two_exp:
    def __init__(self, **kargs):
        self.starting_guess = [0.05, 100, 1.0, 10]
        self.parameter_names = ["mass", "amp", "mass2", "amp2"]
        self.description = "two_exp"
        self.template = "{1: f}exp(-{0: f}*t)(1+{3: f}exp(-{2: f}^2*t)"

    def formula(self, v, x):
        return (v[1] * np.exp((-1.0) * v[0] * x)*(1.0 + v[3]*np.exp((-1.0)*(v[2]**2)*x)))


class periodic_two_exp:
    def __init__(self, Nt=None):
        self.starting_guess = [0.05, 100, 1.0, 10]
        self.parameter_names = ["mass", "amp", "mass2", "amp2"]
        self.description = "periodic_two_exp"
        self.template = "{1: f}exp(-{0: f}*t)(1+{3: f}exp(-{2: f}^2*t)"
        self.Nt = Nt
        while not self.Nt:
            try:
                self.Nt = int(raw_input('Time period not specified, please enter Nt:'))
            except ValueError:
                print "Not a valid number"

    def formula(self, v, x):
                return ((v[1]*np.exp((-1.0)*v[0]*x)*(1.0 + v[3]*np.exp((-1.0)*(v[2]**2)*x))) +
                        (v[1]*np.exp(v[0]*(x-(self.Nt)))*(1.0 + v[3]*np.exp((v[2]**2)*(x-(self.Nt))))))


class jlab:
    def __init__(self, **kargs):
        self.starting_guess = [0.5, 100, 1.0, 10]
        self.parameter_names = ["mass", "amp", "mass2"]
        self.description = "jlab"
        self.template = "{1: f}exp(-{0: f}*t)+{1: f}exp(-{2: f}*t)"

    def formula(self, v, x):
        return ((1 - v[1]) * np.exp((-1.0) * v[0] * x))+(v[1] * np.exp((-1.0) * v[2] * x))


class cosh_const:
    def __init__(self, Nt=None):
        self.starting_guess = [0.1, 10, 0.0]
        self.parameter_names = ["mass", "amp", "const"]
        self.description = "cosh+const"
        self.Nt = Nt
        while not self.Nt:
            try:
                self.Nt = int(raw_input('Time period not specified, please enter Nt:'))
            except ValueError:
                print "Not a valid number"
        self.template = "{1: f}Cosh(-{0: f}*(t-%d/2))+{2: f}" % self.Nt

    def formula(self, v, x):
        return (v[1] * np.cosh((-1.0)*v[0]*((x-(self.Nt/2.0)))))+v[2]
