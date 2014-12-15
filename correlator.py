import configtimeobj
import jackknife
import math
import vev
import logging


class Correlator(configtimeobj.Cfgtimeobj):

    vevdata = None

    # op1 = None
    # op2 = None

    asv = None
    jkasv = None

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
        self.vev1 = vev.Vev(vev1)
        self.vev2 = vev.Vev(vev2)
        super(Correlator, self).__init__(datadict)

    @classmethod
    def fromDataDicts(cls, corr, vev1, vev2):
        """ Create a correlator from a dictionaries for the correlator, and vevs
        """
        return cls(corr, vev1, vev2)

    def verify(self):
        logging.debug("verifying correlator")

        assert self.configs == self.vev1.configs
        assert self.configs == self.vev2.configs

        super(Correlator, self).verify()

    def writefullfile(self, filename, comp=False):
        logging.debug("writting vevs to %s", filename + ".vev1,2")
        self.vev1.writefullfile(filename + ".vev1")
        self.vev2.writefullfile(filename + ".vev2")
        logging.debug("writting correlator to %s", filename + ".cor")
        super(Correlator, self).writefullfile(filename + ".cor", comp=comp)

    def average_sub_vev(self):
        if not self.asv:
            vev1 = self.vev1.average()
            vev2 = self.vev2.average()
            self.asv = {t: corr - vev1 * vev2
                        for t, corr in self.average_over_configs().iteritems()}
        return self.asv

    def jackknife_average_sub_vev(self):
        if not self.jkasv:

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

    def cosh_effective_mass(self, dt):
        asv = self.average_sub_vev()
        emass = {}
        for t in self.times[dt:-dt]:
            try:
                emass[t] = (1.0 / float(dt))*math.acosh((asv[t+dt] + asv[t-dt])/(2.0*asv[t]))
            except ValueError:
                logging.debug("invalid argument to acosh, setting to zero")
                emass[t] = float('NaN')
            except KeyError:
                logging.error("index out of range")
            except ZeroDivisionError:
                logging.error("Div by zero either dt:{} or average value sub vev {}".format(dt,asv[t]))
                emass[t] = float('NaN')
        return emass

    def cosh_const_effective_mass(self, dt):
        asv = self.average_sub_vev()
        asvt = {t: asv[t+dt]-asv[t] for t in self.times[:-dt] }
        emass = {}
        for t in self.times[dt:-(dt+dt)]:
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

    def cosh_effective_mass_errors(self, dt):

        jkasv = self.jackknife_average_sub_vev()
        jkemass = {}
        for cfg in self.configs:
            asvc = jkasv[cfg]
            emass = {}
            for t in self.times[dt:-dt]:
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
        effmass_dt = self.cosh_effective_mass(dt)
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

    def writeemass(self, filename, dt=3, header=None, cosh=True):
        logging.info("Writing emass{} to {}".format(dt,filename))
        if cosh:
            emass = self.cosh_effective_mass(dt)
            error = self.cosh_effective_mass_errors(dt)
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
        self.sums = None
        self.times = new_times

    def check_symmetric(self, sigma=1.0, anti=False):
        asv = self.average_sub_vev()
        errors = self.jackknifed_errors()
        print asv
        print errors
        seperations = [t for t in sorted(self.times) if t>0]
        print seperations
        max_asymmetry = 0
        for tf, tb in zip(seperations, reversed(seperations)):
            if anti:
                asymmetry = abs(asv[tf] + asv[tb]) / (errors[tf]+errors[tb])
            else:
                asymmetry = abs(asv[tf] - asv[tb]) / (errors[tf]+errors[tb])
            logging.debug("asymmetry in correlator({},{}): {}".format(tf,tb,asymmetry))
            max_asymmetry = max(asymmetry,max_asymmetry)
            if max_asymmetry>sigma:
                logging.error("correlator is not symmetric within {}sigma".format(sigma))
                logging.error("C({}) - C({}) = {}, E({}) E({})".format(tf , tb, asv[tf] - asv[tb], errors[tf],errors[tb],))
                return False
        logging.info("Max asymmetry in correlator: {}sigma".format(max_asymmetry))
        return True

    def make_symmetric(self, sigma=1.0, anti=False):
        logging.info("original times {}-{}".format( min(self.times), max(self.times) ) )
        asv = self.average_sub_vev()
        errors = self.jackknifed_errors()

        # new_times = [t for t in self.times if asv[t] - 2.0 * errors[t] > 0.0]
        period = len(self.times)
        print period
        seperations = list(range(1,period/2+1))
        print seperations
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
                    self.data[cfg][t] = (self.data[cfg][t] + self.data[cfg][period-t])/2.0
            else:
                logging.info("Removing data for t={}".format(t))
                for cfg in self.configs:
                    del self.data[cfg][t]
        logging.info("Correlator made symetric, largest change was {}".format(biggest_change))
        self.times = seperations
        self.asv = None
        self.jkasv = None
