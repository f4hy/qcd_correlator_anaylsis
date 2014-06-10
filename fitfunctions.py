import numpy as np
from fit_parents import periodic, mass_amp, mass_amp_const, twice_mass_amp, twice_mass_amp_const, InvalidFit


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


class periodic_exp_subtracted(mass_amp, periodic):
    def __init__(self, Nt=None):
        super(periodic_exp_subtracted, self).__init__()
        self.setNt(Nt)
        self.description = "fwd-back-exp subtracted"
        self.subtract = 1
        self.template = "{1: f}(exp(-{0: f}*t)+exp(-{0: f}*(t-%d)) - {1: f}(exp(-{0: f}*%d)+exp(-{0: f}*(%d-%d)) - " % (self.Nt, self.subtract, self.subtract, self.Nt)
        self.fallback = None

    def formula(self, v, x):
        return (v[1] * (np.exp((-1.0) * v[0] * x) + np.exp(v[0] * (x-(self.Nt))))) - (v[1] * (np.exp((-1.0) * v[0] * self.subtract) + np.exp(v[0] * (self.subtract-(self.Nt)))))


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
        self.fallback = "single_exp"

    def formula(self, v, x):
        return (v[1] * np.exp((-1.0) * v[0] * x)*(1.0 + v[3]*np.exp((-1.0)*(v[2]**2)*x)))


class periodic_two_exp(twice_mass_amp, periodic):
    def __init__(self, Nt=None):
        super(periodic_two_exp, self).__init__()
        self.setNt(Nt)
        self.description = "periodic_two_exp"
        self.template = "{1: f}exp(-{0: f}*t)(1+{3: f}exp(-{2: f}^2*t)"
        self.fallback = "periodic_exp"


    def formula(self, v, x):
                return ((v[1]*np.exp((-1.0)*v[0]*x)*(1.0 + v[3]*np.exp((-1.0)*(v[2]**2)*x))) +
                        (v[1]*np.exp(v[0]*(x-(self.Nt)))*(1.0 + v[3]*np.exp((v[2]**2)*(x-(self.Nt))))))  # noqa


class periodic_two_exp_subtracted(twice_mass_amp, periodic):
    def __init__(self, Nt=None):
        super(periodic_two_exp_subtracted, self).__init__()
        self.setNt(Nt)
        self.description = "periodic_two_exp_subtracted"
        self.subtract = 1
        self.template = "{1: f}exp(-{0: f}*t)(1+{3: f}exp(-{2: f}^2*t)"
        self.fallback = None


    def formula(self, v, x):
                return (((v[1]*np.exp((-1.0)*v[0]*x)*(1.0 + v[3]*np.exp((-1.0)*(v[2]**2)*x))) +
                        (v[1]*np.exp(v[0]*(x-(self.Nt)))*(1.0 + v[3]*np.exp((v[2]**2)*(x-(self.Nt)))))) -
                ((v[1]*np.exp((-1.0)*v[0]*self.subtract)*(1.0 + v[3]*np.exp((-1.0)*(v[2]**2)*self.subtract))) +
                 (v[1]*np.exp(v[0]*(self.subtract-(self.Nt)))*(1.0 + v[3]*np.exp((v[2]**2)*(self.subtract-(self.Nt)))))))  # noqa


class periodic_two_exp_const(twice_mass_amp_const, periodic):
    def __init__(self, Nt=None):
        super(periodic_two_exp_const, self).__init__()
        self.setNt(Nt)
        self.description = "periodic_two_exp_const"
        self.template = "{1: f}exp(-{0: f}*t)(1+{3: f}exp(-{2: f}^2*t)+{4: f}"
        self.fallback = "periodic_exp_const"


    def formula(self, v, x):
                return ((v[1]*np.exp((-1.0)*v[0]*x)*(1.0 + v[3]*np.exp((-1.0)*(v[2]**2)*x))) +
                        (v[1]*np.exp(v[0]*(x-(self.Nt)))*(1.0 + v[3]*np.exp((v[2]**2)*(x-(self.Nt))))))+v[4]  # noqa


# def pade_guess(*args, **kargs):
#     first_two = fit_parents.massamp_guess(args[0], args[1])
#     return [first_two[0], first_two[1], 0.0]

# mass_bounds = (0.005, 2.0)
# amp_bounds = (0.0, 1000.0)


# class pade:
#     """ exp( - E * t )  *    A /  ( 1 +   a1* t + a2 * t^2 ...  ) """
#     def __init__(self, Nt=None):
#         self.starting_guess = pade_guess
#         self.bounds = [mass_bounds, amp_bounds, (-100.0,100.0)]
#         self.parameter_names = ["mass", "amp", "B"]
#         self.description = "Pade"
#         self.template = "{1: f}exp(-{0: f}*t)/(1+{3: f}t)"

#     def formula(self, v, x):
#         return (v[1]*np.exp((-1.0)*v[0]*x)) / (1.0 + v[2]*x)

#     def my_cov_fun(self, mass, amp, B):
#         vect = self.aoc - self.formula((mass, amp, B), self.times)
#         return vect.dot(self.inv_cov).dot(vect)

#     def custom_minuit(self, data, invmatrix, times, guess):
#         self.aoc = data
#         self.inv_cov = invmatrix
#         self.times = times
#         m = Minuit(self.my_cov_fun, mass=guess[0], amp=guess[1], B=guess[2],
#                    print_level=0, pedantic=False, limit_mass=mass_bounds)
#         return m

# class jlab:
#     def __init__(self, **kargs):
#         self.starting_guess = twoexp_sqr_guess
#         self.parameter_names = ["mass", "amp", "mass2"]
#         self.description = "jlab"
#         self.template = "{1: f}exp(-{0: f}*t)+{1: f}exp(-{2: f}*t)"

#     def formula(self, v, x):
#         return ((1 - v[1]) * np.exp((-1.0) * v[0] * x))+(v[1] * np.exp((-1.0) * v[2] * x))
