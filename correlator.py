import configtimeobj
import jackknife
import math
import vev
import logging
import newton

class Correlator(configtimeobj.Cfgtimeobj):

    made_symmetric = False
    symmetry = None
    vevdata = None
    period = None
    # op1 = None
    # op2 = None

    asv = None
    jkasv = None

    emass_skip_times = []

    @classmethod
    def fromOpvalCTO(cls, opval1, opval2, dts=None):
        if not opval1.compatible(opval2):
            raise Exception("Incompatible opbjects")
        if dts is None:
            logging.warning("No dts given, using all")
            dts = opval1.times
        times = opval1.times
        configs = opval1.configs
        numtimes = opval1.numtimes
        get1 = opval1.get
        get2 = opval2.get

        data = {}
        for cfg in configs:
            opc1 = get1(config=cfg)
            opc2 = get2(config=cfg)
            inerdata = {}
            for dt in dts:
                #acc = 0.0
                acc = math.fsum(opc1[(t + dt) % numtimes] * opc2[t] for t in times)
                # for t in times:
                #     shifted_t = (t + dt) % numtimes
                #     acc += opc1[shifted_t] * opc2[t]
                inerdata[dt] = acc / float(numtimes)
            data[cfg] = inerdata

        vev1 = opval1.average_over_times()
        vev2 = opval2.average_over_times()

        # data = {cfg:
        #         {dt:
        #          math.fsum(get1(config=cfg, time=((t+dt)%numtimes))*get2(config=cfg, time=t)
        #                    for t in times)/float(numtimes)
        #          for dt in dts }
        #          for cfg in configs}

        return cls(data, vev1, vev2)

    def __init__(self, datadict, vev1, vev2):
        if vev1 is not None and vev2 is not None:
            self.vev1 = vev.Vev(vev1)
            self.vev2 = vev.Vev(vev2)
        else:
            self.vev1 = None
            self.vev2 = None
        super(Correlator, self).__init__(datadict)

    @classmethod
    def fromDataDicts(cls, corr, vev1, vev2):
        """ Create a correlator from a dictionaries for the correlator, and vevs
        """
        return cls(corr, vev1, vev2)

    def verify(self):
        logging.debug("verifying correlator")

        if self.vev1 is not None:
            assert self.configs == self.vev1.configs
            assert self.configs == self.vev2.configs

        super(Correlator, self).verify()

    def writefullfile(self, filename, comp=False):
        if self.vev1 is not None:
            logging.debug("writting vevs to %s", filename + ".vev1,2")
            self.vev1.writefullfile(filename + ".vev1")
            self.vev2.writefullfile(filename + ".vev2")
        logging.debug("writting correlator to %s", filename + ".cor")
        super(Correlator, self).writefullfile(filename + ".cor", comp=comp)

    def average_sub_vev(self):
        if not self.asv:
            if self.vev1 is None:
                self.asv = self.average_over_configs()
            else:
                vev1 = self.vev1.average()
                vev2 = self.vev2.average()
                self.asv = {t: corr - vev1 * vev2
                            for t, corr in self.average_over_configs().iteritems()}
        return self.asv

    def jackknife_average_sub_vev(self):
        if not self.jkasv:
            if self.vev1 is None:
                self.jkasv = self.jackknifed_averages()
            else:

                jkvev1 = self.vev1.jackknife()
                jkvev2 = self.vev2.jackknife()
                #corrjk = self.jackknifed_averages()
                jk = configtimeobj.Cfgtimeobj.fromDataDict(self.jackknifed_averages())
                self.jkasv = {c: {t: jk.get(config=c, time=t) - jkvev1[c] * jkvev2[c]
                                  for t in self.times}
                              for c in self.configs}

        return self.jkasv

    def jackknifed_errors(self):
        jk = configtimeobj.Cfgtimeobj.fromDataDict(self.jackknife_average_sub_vev())
        asv = self.average_sub_vev()
        return {t: jackknife.errorbars(asv[t], jk.get(time=t)) for t in self.times}

    def prune_invalid(self, sigma=1, delete=False):
        logging.info("original times {}-{}".format( min(self.times), max(self.times) ) )
        asv = self.average_sub_vev()
        errors = self.jackknifed_errors()

        # new_times = [t for t in self.times if asv[t] - 2.0 * errors[t] > 0.0]
        new_times = []
        for t in self.times:
            if (abs(asv[t]) - errors[t]*sigma) > 0.0:
                new_times.append(t)
            else:
                break
        if len(new_times) < 2:
            logging.error("this correlator has less than 2 valid times!! Skipping pruing")
        else:
            logging.info("pruned times down to {}-{}".format( min(new_times), max(new_times) ) )
            if delete:
                for removed_times in [t for t in self.times if t > max(new_times)]:
                    logging.info("removing data for time {}".format(removed_times))
                    for cfg in self.configs:
                        del self.data[cfg][removed_times]
            self.times = new_times
            self.asv = None
            self.jkasv = None

    def effective_mass(self, dt):
        asv = self.average_sub_vev()
        emass = {}
        for t in self.times[:-dt]:
            try:
                emass[t] = (1.0 / float(dt)) * math.log(asv[t] / asv[t + dt])
            except ValueError:
                logging.debug("invalid argument to log at t={}, setting to NaN".format(t))
                emass[t] = float('NaN')
            except KeyError:
                logging.error("index out of range")
            except ZeroDivisionError:
                logging.error("Div by zero either dt:{} or average value sub vev {}".format(dt,asv[t + dt]))
                emass[t] = float('NaN')
        return emass

    def effective_amp(self, dt):
        asv = self.average_sub_vev()

        emass = {}
        eamp = {}
        for t in self.times[:-dt]:
            try:
                emass[t] = (1.0 / float(dt)) * math.log(asv[t] / asv[t + dt])
                eamp[t] = asv[t] * math.exp(emass[t]*t)
            except ValueError:
                logging.debug("invalid argument to log at t={}, setting to NaN".format(t))
                emass[t] = float('NaN')
                eamp[t] = float('NaN')
            except KeyError:
                logging.error("index out of range")
            except ZeroDivisionError:
                logging.error("Div by zero either dt:{} or average value sub vev {}".format(dt,asv[t + dt]))
                emass[t] = float('NaN')
                eamp[t] = float('NaN')
        #exit()
        return eamp


    def periodic_effective_mass(self, dt, fast=True, period=None):
        if self.symmetry is None:
            logging.warning("Called periodic effective mass without symmetry determined")
            self.determine_symmetry()
            if self.symmetry is None:
                logging.error("Called periodic effective mass and symmetry can not be found")
                raise RuntimeError("Could not determine symmetry")
        if self.symmetry == "symmetric":
            logging.info("Calling cosh emass")
            return self.cosh_effective_mass(dt, fast=fast, period=period)
        if self.symmetry == "anti-symmetric":
            logging.info("Calling sinh emass")
            return self.sinh_effective_mass(dt, fast=fast, period=period)

        logging.error("Symmetry is not 'symmetric' nor 'anti-symmetric'")
        raise RuntimeError("Could not determine symmetry")

    def periodic_effective_mass_errors(self, dt, fast=True, period=None):

        if self.symmetry is None:
            logging.warning("Called periodic effective mass without symmetry determined")
            self.determine_symmetry()
            if self.symmetry is None:
                logging.error("Called periodic effective mass and symmetry can not be found")
                raise RuntimeError("Could not determine symmetry")
        if self.symmetry == "symmetric":
            return self.cosh_effective_mass_errors(dt, fast=fast, period=period)
        if self.symmetry == "anti-symmetric":
            return self.sinh_effective_mass_errors(dt, fast=fast, period=period)

        logging.error("Symmetry is not 'symmetric' nor 'anti-symmetric'")
        raise RuntimeError("Could not determine symmetry")



    def cosh_effective_mass(self, dt, fast=True, period=None):
        if fast: logging.info("cosh emass computed fast method")

        T = self.period_check(period)
        asv = self.average_sub_vev()
        emass = {}
        for t in self.times[dt:-dt]:
            if t in self.emass_skip_times:
                emass[t] = 0.0
                continue
            try:
                guess = (1.0 / float(dt))*math.acosh((asv[t+dt] + asv[t-dt])/(2.0*asv[t]))
                if fast:
                    emass[t] = guess
                else:
                    emass[t] = newton.newton_cosh_for_m(t,t+dt,asv, guess, T)
            except ValueError:
                logging.debug("invalid argument to acosh, setting to zero")
                emass[t] = float('NaN')
            except KeyError:
                logging.error("index out of range")
            except ZeroDivisionError:
                logging.error("Div by zero either dt:{} or average value sub vev {}".format(dt,asv[t]))
                emass[t] = float('NaN')
        return emass

    def sinh_effective_mass(self, dt, fast=True, period=None):
        if fast: logging.info("sinh emass computed fast method")
        T = self.period_check(period)
        asv = self.average_sub_vev()
        emass = {}
        for t in self.times[dt:-dt]:
            if t in self.emass_skip_times:
                emass[t] = 0.0
                continue
            try:
                guess = (1.0 / float(dt))*math.asinh((asv[t+dt] + asv[t-dt])/(2.0*asv[t]))
                if fast:
                    emass[t] = guess
                else:
                    emass[t] = newton.newton_sinh_for_m(t,t+dt,asv, guess, T)
            except ValueError:
                logging.debug("invalid argument to asinh, setting to zero")
                emass[t] = float('NaN')
            except KeyError:
                logging.error("index out of range")
            except ZeroDivisionError:
                logging.error("Div by zero either dt:{} or average value sub vev {}".format(dt,asv[t]))
                emass[t] = float('NaN')
        return emass


    def cosh_effective_amp(self, dt, period, mass):
        asv = self.average_sub_vev()
        emass = {}
        eamp = {}
        for t in self.times[dt:-dt]:
            if t in self.emass_skip_times:
                emass[t] = 0.0
                continue
            try:
                # emass[t] = (1.0 / float(dt))*math.acosh((asv[t+dt] + asv[t-dt])/(2.0*asv[t]))
                # eamp[t] = asv[t] / (math.exp(emass[t]*t)+math.exp(emass[t]*(period-t)))
                eamp[t] = asv[t] / (math.exp(mass*t)+math.exp(mass*(period-t)))
            except ValueError:
                logging.debug("invalid argument to acosh, setting to zero")
                emass[t] = float('NaN')
                eamp[t] = asv[t] * math.exp(emass[t]*t)
            except KeyError:
                logging.error("index out of range")
            except ZeroDivisionError:
                logging.error("Div by zero either dt:{} or average value sub vev {}".format(dt,asv[t]))
                emass[t] = float('NaN')
                eamp[t] = asv[t] * math.exp(emass[t]*t)
        return eamp


    def cosh_const_effective_mass(self, dt):
        asv = self.average_sub_vev()
        asvt = {t: asv[t+dt]-asv[t] for t in self.times[:-dt] }
        emass = {}
        for t in self.times[dt:-(dt+dt)]:
            if t in self.emass_skip_times:
                emass[t] = 0.0
                continue
            try:
                emass[t] = (1.0 / float(dt))*math.acosh((asvt[t+dt] + asvt[t-dt])/(2.0*asvt[t]))
            except ValueError:
                logging.debug("invalid argument to acosh, setting to nan")
                emass[t] = float('NaN')
            except KeyError:
                logging.error("index out of range")
            except ZeroDivisionError:
                logging.error("Div by zero either dt:{} or average value sub vev {}".format(dt,asv[t]))
                emass[t] = 0.0
        return emass


    def effective_mass_errors(self, dt):

        jkasv = self.jackknife_average_sub_vev()
        jkemass = {}
        for cfg in self.configs:
            asvc = jkasv[cfg]
            emass = {}
            for t in self.times[:-dt]:
                if t in self.emass_skip_times:
                    emass[t] = 0.0
                    continue

                try:
                    emass[t] = (1.0 / float(dt)) * math.log(asvc[t] / asvc[t + dt])
                except ValueError:
                    #logging.debug("invalid argument to log, setting to zero")
                    emass[t] = 0.0
                except ZeroDivisionError:
                    logging.debug("div by zero, setting to zero")
                    emass[t] = 0.0
                except KeyError:
                    logging.error("index out of range")
            jkemass[cfg] = emass
        jkemassobj = configtimeobj.Cfgtimeobj.fromDataDict(jkemass)
        effmass_dt = self.effective_mass(dt)
        return {t: jackknife.errorbars(effmass_dt[t], jkemassobj.get(time=t))
                for t in self.times[:-dt]}

    def cosh_effective_mass_errors(self, dt, fast=True, period=None):
        if fast: logging.info("cosh emass computed fast method")
        period = self.period

        jkasv = self.jackknife_average_sub_vev()
        jkemass = {}
        for cfg in self.configs:
            asvc = jkasv[cfg]
            emass = {}
            for t in self.times[dt:-dt]:
                if t in self.emass_skip_times:
                    emass[t] = 0.0
                    continue

                try:
                    guess = (1.0 / float(dt))*math.acosh((asvc[t+dt] + asvc[t-dt])/(2.0*asvc[t]))
                    if fast:
                        emass[t] = guess
                    else:
                        emass[t] = newton.newton_cosh_for_m(t,t+dt,asvc, guess,period)
                except ValueError:
                    #logging.debug("invalid argument to log, setting to zero")
                    emass[t] = 0.0
                except ZeroDivisionError:
                    logging.debug("div by zero, setting to zero")
                    emass[t] = 0.0
                except KeyError:
                    logging.error("index out of range")
            jkemass[cfg] = emass
        jkemassobj = configtimeobj.Cfgtimeobj.fromDataDict(jkemass)
        effmass_dt = self.cosh_effective_mass(dt, fast=fast, period=period)
        return {t: jackknife.errorbars(effmass_dt[t], jkemassobj.get(time=t))
                for t in self.times[dt:-dt]}



    def sinh_effective_mass_errors(self, dt, fast=True, period=None):
        if fast: logging.info("sinh emass computed fast method")

        T = self.period_check(period)

        jkasv = self.jackknife_average_sub_vev()
        jkemass = {}
        for cfg in self.configs:
            asvc = jkasv[cfg]
            emass = {}
            for t in self.times[dt:-dt]:
                if t in self.emass_skip_times:
                    emass[t] = 0.0
                    continue
                try:
                    guess = (1.0 / float(dt))*math.asinh((asvc[t+dt] + asvc[t-dt])/(2.0*asvc[t]))
                    if fast:
                        emass[t] = guess
                    else:
                        emass[t] = newton.newton_sinh_for_m(t,t+dt,asvc, guess, T)
                except ValueError:
                    #logging.debug("invalid argument to log, setting to zero")
                    emass[t] = 0.0
                except ZeroDivisionError:
                    logging.debug("div by zero, setting to zero")
                    emass[t] = 0.0
                except KeyError:
                    logging.error("index out of range")
            jkemass[cfg] = emass
        jkemassobj = configtimeobj.Cfgtimeobj.fromDataDict(jkemass)
        effmass_dt = self.sinh_effective_mass(dt, fast=fast, period=T)
        return {t: jackknife.errorbars(effmass_dt[t], jkemassobj.get(time=t))
                for t in self.times[dt:-dt]}


    def cosh_const_effective_mass_errors(self, dt):

        jkasv = self.jackknife_average_sub_vev()
        jkemass = {}
        for cfg in self.configs:
            asvc = jkasv[cfg]
            emass = {}
            asvct = {t: asvc[t+dt] - asvc[t] for t in self.times[:-dt]}
            for t in self.times[dt:-(dt+1)]:
                if t in self.emass_skip_times:
                    emass[t] = 0.0
                    continue
                try:
                    emass[t] = (1.0 / float(dt))*math.acosh((asvc[t+dt] + asvc[t-dt])/(2.0*asvc[t]))
                except ValueError:
                    #logging.debug("invalid argument to log, setting to zero")
                    emass[t] = 0.0
                except ZeroDivisionError:
                    logging.debug("div by zero, setting to zero")
                    emass[t] = 0.0
                except KeyError:
                    logging.error("index out of range")
            jkemass[cfg] = emass
        jkemassobj = configtimeobj.Cfgtimeobj.fromDataDict(jkemass)
        effmass_dt = self.cosh_const_effective_mass(dt)
        return {t: jackknife.errorbars(effmass_dt[t], jkemassobj.get(time=t))
                for t in self.times[dt:-(dt+dt)]}


    def reduce_to_bins(self, n):
        reduced = {}
        binedvev1 = {}
        binedvev2 = {}
        for i, b in enumerate(self.bins(n)):
            # logging.debug("bin:")
            # logging.debug(b)
            size = float(len(b))
            reduced[i] = {t: math.fsum((self.get(config=c, time=t) for c in b)) / size
                          for t in self.times}

            if self.vev1 is None:
                binedvev1 = binnedvev2 = None
            else:
                binedvev1[i] = math.fsum((self.vev1[c] for c in b)) / size
                binedvev2[i] = math.fsum((self.vev2[c] for c in b)) / size

        logging.info("Binned from %d, reduced to %d bins", self.numconfigs, len(reduced.keys()))
        # Make a new correlator for the bined data
        return Correlator.fromDataDicts(reduced, binedvev1, binedvev2)

    def bins(self, n):
        """ Yield successive n-sized chunks from configs.
        """
        if self.numconfigs % n is not 0:
            logging.warning("Bin size %d not factor of num configs %d !!!", n, self.numconfigs)
        for i in xrange(0, self.numconfigs, n):
            yield self.configs[i:i + n]

    def writeasv(self, filename, header=None):
        logging.info("Writing asv_cor to {}".format(filename))
        asv = self.average_sub_vev()
        error = self.jackknifed_errors()
        with open(filename, 'w') as outfile:
            if header:
                outfile.write(header)
                outfile.write("\n")
            for t, a in asv.iteritems():
                outfile.write("{!r},   {!r}, {!r}\n".format(t, a, error[t]))

    def writeemass(self, filename, dt=3, header=None, periodic=True):
        logging.info("Writing emass{} to {}".format(dt,filename))
        if periodic:
            emass = self.periodic_effective_mass(dt, fast=False)
            error = self.periodic_effective_mass_errors(dt, fast=False)
        with open(filename, 'w') as outfile:
            if header:
                outfile.write(header)
                outfile.write("\n")
            for t, e in emass.iteritems():
                outfile.write("{!r},   {!r}, {!r}\n".format(t, e, error[t]))


    def subtract(self, t):
        logging.info("original times {}-{}".format( min(self.times), max(self.times) ) )

        new_times = [time for time in self.times if time > t]

        for time in new_times:
            logging.info("redefinging correlator data for time {0} as C'({0}) = C({0}) - C({1})".format(time, t))
            for cfg in self.configs:
                self.data[cfg][time] = self.data[cfg][time] - self.data[cfg][t]
        self.asv = None
        self.jkasv = None
        self.average = None
        self.sums = None
        self.vevdata = None
        self.times = new_times

    def period_check(self, period):
        if period is None:
            if self.period:
                return self.period
            else:
                logging.warning("Period assmed to be number of times")
                T = len(self.times)
                if self.made_symmetric:
                    T *= 2
                return T
        else:
            if self.period is not None:
                if period != self.period:
                    logging.error("passed period does not agree with internal period")
                    raise RuntimeError("passed period does not agree with internal period")
        return period


    def check_symmetric(self, sigma=1.0, anti=False):
        antistring = "anti-" if anti else ""
        asymmetry = self.check_symmetric_sigma(anti=anti)
        if asymmetry > sigma:
            logging.error("correlator is not {}symmetric within {}sigma is {}".format(antistring, sigma,asymmetry))
            return False
        return True

    def check_symmetric_sigma(self, anti=False):
        antistring = "anti-" if anti else "     "
        asv = self.average_sub_vev()
        errors = self.jackknifed_errors()
        seperations = [t for t in sorted(self.times) if t>0]
        max_asymmetry = 0
        for tf, tb in zip(seperations, reversed(seperations)):
            if anti:
                asymmetry = abs(asv[tf] + asv[tb]) / (errors[tf]+errors[tb])
            else:
                asymmetry = abs(asv[tf] - asv[tb]) / (errors[tf]+errors[tb])
            logging.debug("asymmetry in correlator {}symmetry ({},{}): {}".format(antistring, tf,tb,asymmetry))
            max_asymmetry = max(asymmetry,max_asymmetry)
        logging.info("Max asymmetry in correlator {}symmetry: {}sigma".format(antistring, max_asymmetry))
        return max_asymmetry

    def make_symmetric(self, sigma=1.0):
        if self.made_symmetric:
            logging.error("Already made symmetric!!!")
            raise RuntimeError("Called make symmetric but correlator Already made symmetric!!!")

        corsym = self.determine_symmetry()
        if not corsym:
            logging.error("can not make symmetric, symmetry type undetermined")
            raise RuntimeError("Called make symmetric but correlator unknown symmetry!!!")
        anti = False
        if (corsym == 'anti-symmetric'):
            anti = True

        logging.info("making correlator {}symmetric".format("anti-" if anti else "     "))
        logging.info("original times {}-{}".format( min(self.times), max(self.times) ) )
        asv = self.average_sub_vev()
        errors = self.jackknifed_errors()

        removed_data = []
        # new_times = [t for t in self.times if asv[t] - 2.0 * errors[t] > 0.0]
        if self.period is not None:
            if self.period != len(self.times):
                raise RuntimeError("trying to make symmetric but period does not match times!")
        else:
            self.period = len(self.times)
        period = self.period
        seperations = list(range(1,period/2+1))
        biggest_change = 0.0
        for t in self.times:
            if t in seperations:
                logging.debug("averaging %d %d", t, period-t)
                for cfg in self.configs:
                    prevdata = self.data[cfg][t]
                    if anti:
                        newdata = (self.data[cfg][t] - self.data[cfg][period-t])/2.0
                    else:
                        newdata = (self.data[cfg][t] + self.data[cfg][period-t])/2.0

                    biggest_change = max(biggest_change, abs(prevdata - newdata)/prevdata)
                    self.data[cfg][t] = newdata # (self.data[cfg][t] + self.data[cfg][period-t])/2.0
            else:
                removed_data.append(t)
                for cfg in self.configs:
                    del self.data[cfg][t]
        logging.info("Removed data for t={}".format(repr(removed_data)))

        logging.info("Correlator made symetric, largest change was {}".format(biggest_change))
        self.times = seperations
        self.asv = None
        self.jkasv = None
        self.average = None
        self.sums = None
        self.vevdata = None
        self.made_symmetric = True

    def determine_symmetry(self, recheck=False):
        if self.made_symmetric:
            logging.error("Already made symmetric!!!")
            raise RuntimeError("Called determine symmetry but correlator Already made symmetric!!!")


        if self.symmetry:
            logging.info("correlator symmetry is {}".format(self.symmetry))
            if not recheck:
                return self.symmetry

        sym = self.check_symmetric_sigma()
        asym = self.check_symmetric_sigma(anti=True)

        if min(asym,sym) > 6.0:
            logging.info("No symmetry found within 6sigma!!")
        elif sym < asym:
            self.symmetry = "symmetric"
            logging.info("correlator found to be symmetric!")
        else:
            self.symmetry = "anti-symmetric"
            logging.info("correlator found to be anti-symmetric!")

        return self.symmetry

    def multiply_by_value_dict(self, d):


        logging.warn("Dividing correlator by {}!!!!!".format(d))

        for c in self.configs:
            for t in self.times:
                self.data[c][t] = self.data[c][t]/d[t]

        self.asv = None
        self.jkasv = None
        self.average = None
        self.sums = None
        self.vevdata = None
