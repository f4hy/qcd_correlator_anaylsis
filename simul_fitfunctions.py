import numpy as np
import simul_fit_parents as sfp
import fit_parents as fp
import fitfunctions as singlefitfunctions

class twocor_periodic_twoexp(sfp.shared_twice_mass_amp, fp.periodic, sfp.multirange):
    individual = singlefitfunctions.periodic_two_exp
    def __init__(self, Nt=None, ranges=None):
        super(twocor_periodic_twoexp, self).__init__()
        self.setNt(Nt)
        self.setranges(ranges)
        self.description = "two_cor-fwd-back-exp"
        self.template = "m{0: f}, A1{1: f} A2{2: f}, m_b^2{3: f}, A1b{4: f} A2b{5: f}"
        self.multi = True

    def formula(self, v, x):
        ys = []
        for i in range(len(self.ranges)):
            r = self.ranges[i]
            tx = np.arange(r[0], r[1]+1, 1)
            ys.append( (v[1+i] * (np.exp((-1.0) * v[0] * tx)*(1.0 + v[4+i]*np.exp((-1.0)*(v[3]**2)*tx))))
                    + (v[1+i] *(np.exp(v[0] * (tx-(self.Nt)))*(1.0 + v[4+i]*np.exp((v[3]**2) * (tx-(self.Nt))))))
            )
            # (v[1] * (np.exp((-1.0) * v[0] * x) + np.exp(v[0] * (x-(self.Nt)))))
        return np.concatenate(ys)


class twocor_periodic_exp(sfp.sharedmass_amp, fp.periodic, sfp.multirange):
    individual = singlefitfunctions.periodic_exp
    def __init__(self, Nt=None, ranges=None):
        super(twocor_periodic_exp, self).__init__()
        self.setNt(Nt)
        self.setranges(ranges)
        self.description = "two_cor-fwd-back-exp"
        self.template = "m{0: f}, A1{1: f} A2{2: f}"
        self.multi = True

    def formula(self, v, x):
        ys = []
        for i in range(len(self.ranges)):
            r = self.ranges[i]
            tx = np.arange(r[0], r[1]+1, 1)
            ys.append(v[1+i] * (np.exp((-1.0) * v[0] * tx) + np.exp(v[0] * (tx-(self.Nt)))))
            # (v[1] * (np.exp((-1.0) * v[0] * x) + np.exp(v[0] * (x-(self.Nt)))))
        return np.concatenate(ys)

class twocor_antiperiodic_exp(sfp.sharedmass_amp, fp.periodic, sfp.multirange):
    individual = singlefitfunctions.periodic_exp
    def __init__(self, Nt=None, ranges=None):
        super(twocor_antiperiodic_exp, self).__init__()
        self.setNt(Nt)
        self.setranges(ranges)
        self.description = "two_cor-fwd-back-exp"
        self.template = "m{0: f}, A1{1: f} A2{2: f}"
        self.multi = True

    def formula(self, v, x):
        ys = []
        for i in range(len(self.ranges)):
            r = self.ranges[i]
            tx = np.arange(r[0], r[1]+1, 1)
            ys.append(v[1+i] * (np.exp((-1.0) * v[0] * tx) - np.exp(v[0] * (tx-(self.Nt)))))
            # (v[1] * (np.exp((-1.0) * v[0] * x) + np.exp(v[0] * (x-(self.Nt)))))
        return np.concatenate(ys)

class twocor_antiperiodic_periodic_exp(sfp.sharedmass_amp, fp.periodic, sfp.multirange):
    individual = singlefitfunctions.periodic_exp
    def __init__(self, Nt=None, ranges=None):
        super(twocor_antiperiodic_periodic_exp, self).__init__()
        self.setNt(Nt)
        self.setranges(ranges)
        self.description = "two_cor-fwd-back-exp"
        self.template = "m{0: f}, A1{1: f} A2{2: f}"
        self.multi = True

    def formula(self, v, x):
        ys = []

        r = self.ranges[0]
        tx = np.arange(r[0], r[1]+1, 1)
        ys.append(v[1] * (np.exp((-1.0) * v[0] * tx) - np.exp(v[0] * (tx-(self.Nt)))))
        r = self.ranges[1]
        tx = np.arange(r[0], r[1]+1, 1)
        ys.append(v[2] * (np.exp((-1.0) * v[0] * tx) + np.exp(v[0] * (tx-(self.Nt)))))

        return np.concatenate(ys)
