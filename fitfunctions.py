import numpy as np

mass_bounds = (0.005, 2.0)
amp_bounds = (0.0, 1.0e8)
const_bounds = (-5.0, 1.0e8)

def massamp_guess(cor, tmax, *args):
    dt = 3
    maxt = tmax - dt
    ave = cor.average_sub_vev()
    emass = cor.effective_mass(dt)
    mass_guess = emass[maxt]
    amp_guess = ave[maxt]*np.exp(mass_guess*(maxt))
    return [mass_guess, amp_guess]

def const_guess(cor, tmax, *args):
    dt = 3
    maxt = tmax - dt
    ave = cor.average_sub_vev()
    emass = cor.effective_mass(dt)
    mass_guess = emass[maxt]
    amp_guess = ave[maxt]*np.exp(mass_guess*(maxt))
    return [mass_guess, amp_guess, 0.01]

def twoexp_sqr_guess(cor, tmax, tmin):
    dt = 3
    maxt = tmax - dt
    ave = cor.average_sub_vev()
    emass = cor.effective_mass(dt)
    mass_guess = emass[maxt]
    amp_guess = ave[maxt]*np.exp(mass_guess*(maxt))
    mass2_guess = np.sqrt(emass[tmin])
    amp2_guess = (ave[tmin] - amp_guess*np.exp(-mass_guess*tmin))/(amp_guess*np.exp(-(mass_guess+mass2_guess**2)*tmin))
    return [mass_guess, amp_guess, mass2_guess, amp2_guess]


class cosh:
    def __init__(self, Nt=None):
        self.starting_guess = massamp_guess
        self.bounds = [mass_bounds, amp_bounds]
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
        self.starting_guess = massamp_guess
        self.bounds = [mass_bounds, amp_bounds]
        self.parameter_names = ["mass", "amp"]
        self.description = "exp"
        self.template = "{1: f}exp(-{0: f}*t)"

    def formula(self, v, x):
        return (v[1] * np.exp((-1.0) * v[0] * x))


class periodic_exp:
    def __init__(self, Nt=None):
        self.starting_guess = massamp_guess
        self.bounds = [mass_bounds, amp_bounds]
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

class periodic_exp_const:
    def __init__(self, Nt=None):
        self.starting_guess = const_guess
        self.bounds = [mass_bounds, amp_bounds, const_bounds]
        self.parameter_names = ["mass", "amp", "const"]
        self.description = "fwd-back-exp_const"
        self.Nt = Nt
        while not self.Nt:
            try:
                self.Nt = int(raw_input('Time period not specified, please enter Nt:'))
            except ValueError:
                print "Not a valid number"
        self.template = "{1: f}(exp(-{0: f}*t)+exp(-{0: f}*(t-%d))+{2: f}" % self.Nt

    def formula(self, v, x):
        return (v[1] * (np.exp((-1.0) * v[0] * x) + np.exp(v[0] * (x-(self.Nt)))))+v[2]


class two_exp:
    def __init__(self, **kargs):
        self.starting_guess = twoexp_sqr_guess
        self.bounds = [mass_bounds, amp_bounds, mass_bounds, amp_bounds]
        self.parameter_names = ["mass", "amp", "mass2", "amp2"]
        self.description = "two_exp"
        self.template = "{1: f}exp(-{0: f}*t)(1+{3: f}exp(-{2: f}^2*t)"

    def formula(self, v, x):
        return (v[1] * np.exp((-1.0) * v[0] * x)*(1.0 + v[3]*np.exp((-1.0)*(v[2]**2)*x)))


class periodic_two_exp:
    def __init__(self, Nt=None):
        self.starting_guess = twoexp_sqr_guess
        self.bounds = [mass_bounds, amp_bounds, mass_bounds, amp_bounds]
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
        self.starting_guess = twoexp_sqr_guess
        self.parameter_names = ["mass", "amp", "mass2"]
        self.description = "jlab"
        self.template = "{1: f}exp(-{0: f}*t)+{1: f}exp(-{2: f}*t)"

    def formula(self, v, x):
        return ((1 - v[1]) * np.exp((-1.0) * v[0] * x))+(v[1] * np.exp((-1.0) * v[2] * x))


class cosh_const:
    def __init__(self, Nt=None):
        self.starting_guess = const_guess
        self.bounds = [mass_bounds, amp_bounds, const_bounds]
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
