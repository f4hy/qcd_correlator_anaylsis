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

    def prune_invalid(self):
        asv = self.average_sub_vev()
        errors = self.jackknifed_errors()

        # new_times = [t for t in self.times if asv[t] - 2.0 * errors[t] > 0.0]
        new_times = []
        for t in self.times:
            if asv[t] - 2.0 * errors[t] > 0.0:
                new_times.append(t)
            else:
                break
        print new_times
        self.times = new_times
        self.asv = None
        self.jkasv = None
        print self.average_sub_vev()

    def effective_mass(self, dt):
        asv = self.average_sub_vev()
        emass = {}
        for t in self.times[:-dt]:
            try:
                emass[t] = (1.0 / float(dt)) * math.log(asv[t] / asv[t + dt])
            except ValueError:
                #logging.debug("invalid argument to log, setting to zero")
                emass[t] = 0.0
            except KeyError:
                logging.error("index out of range")
            except ZeroDivisionError:
                logging.error("Div by zero either dt:{} or average value sub vev {}".format(dt,asv[t + dt]))
                emass[t] = 0.0
        return emass

    def cosh_effective_mass(self, dt):
        asv = self.average_sub_vev()
        emass = {}
        for t in self.times[dt:-dt]:
            try:
                emass[t] = (1.0 / float(dt))*math.acosh((asv[t+dt] + asv[t-dt])/(2.0*asv[t]))
            except ValueError:
                logging.error("invalid argument to acosh, setting to zero")
                emass[t] = 0.0
            except KeyError:
                logging.error("index out of range")
            except ZeroDivisionError:
                logging.error("Div by zero either dt:{} or average value sub vev {}".format(dt,asv[t]))
                emass[t] = 0.0
        return emass

    def cosh_const_effective_mass(self, dt):
        asv = self.average_sub_vev()
        asvt = {t: asv[t+dt]-asv[t] for t in self.times[:-dt] }
        emass = {}
        for t in self.times[dt:-(dt+dt)]:
            try:
                emass[t] = (1.0 / float(dt))*math.acosh((asvt[t+dt] + asvt[t-dt])/(2.0*asvt[t]))
            except ValueError:
                logging.error("invalid argument to acosh, setting to nan")
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
