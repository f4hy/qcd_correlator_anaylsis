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
            dts = cls.times = opval1.times
        else:
            cls.times = dts
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
                acc = 0.0
                for t in times:
                    #acc += get1(config=cfg, time=((t+dt)%numtimes))*get2(config=cfg, time=t)
                    shifted_t = (t + dt) % numtimes
                    acc += opc1[shifted_t] * opc2[t]
                inerdata[dt] = acc / float(numtimes)
            data[cfg] = inerdata

        cls.vev1 = vev.Vev(opval1.average_over_times())
        cls.vev2 = vev.Vev(opval2.average_over_times())

        # data = {cfg:
        #         {dt:
        #          math.fsum(get1(config=cfg, time=((t+dt)%numtimes))*get2(config=cfg, time=t)
        #                    for t in times)/float(numtimes)
        #          for dt in dts }
        #          for cfg in configs}

        return cls(data)

    @classmethod
    def fromDataDicts(cls, corr, vev1, vev2):
        """ Create a correlator from a dictionaries for the correlator, and vevs
        """

        cls.vev1 = vev.Vev(vev1)
        cls.vev2 = vev.Vev(vev2)
        return cls(corr)

    def verify(self):
        logging.debug("verifying correlator")

        assert self.configs == self.vev1.configs
        assert self.configs == self.vev2.configs

        super(Correlator, self).verify()

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
