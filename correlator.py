import configtimeobj
import numpy as np
import jackknife
import math

class Correlator(configtimeobj.Cfgtimeobj):

    # config_sums = None

    vevdata = None

    op1 = None
    op2 = None
    
    # vevsubtracted = False

    asv = None
    jkasv = None
    #pass

    @classmethod
    def fromOpvalCTO(corr,opval1,opval2,dts = None):
        if not opval1.compatible(opval2):
            raise Exception("Incompatible opbjects")
        if dts is None:
            print "No dts given, using all"
            dts = corr.times = opval1.times
        else:
            corr.times = dts
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
                    #acc += get1(config=cfg,time=((t+dt)%numtimes))*get2(config=cfg,time=t)
                    shifted_t = (t+dt)%numtimes
                    acc += opc1[shifted_t]*opc2[t]
                inerdata[dt] = acc/float(numtimes)
            data[cfg] = inerdata
        
        corr.op1 = opval1
        corr.op2 = opval2

        # data = {cfg:
        #         {dt:
        #          math.fsum(get1(config=cfg,time=((t+dt)%numtimes))*get2(config=cfg,time=t) for t in times)/float(numtimes)
        #          for dt in dts }
        #          for cfg in configs}
        
        return corr(data)
        
    def verify(self):
        print "verifying correlator"

        if self.op1 is not None:
            self.op1.verify()
            self.op2.verify()
        else:
            print "verify vevs, not implemented"
        # print float(self[self.configs[0]][self.times[0]])
        # print len(self[self.configs[0]][self.times[0]])
        # if len(self[self.configs[0]][self.times[0]])
        # exit()
        super(Correlator, self).verify()


        
        
    # def subvev(self):
    #     if self.vevsubtracted:
    #         raise ValueError("vev already subtracted!")
    #     vevs = self.vevs()
    #     self.data = {cfg: {t: self.get(time=t,config=cfg)-vevs[cfg] for t in self.times} for cfg in self.configs}
    #     self.vevsubtracted = True
    #     self.sums = None
    #     return 0

    # def average_vev(self):
    #     return float(sum(self.vevs().values()))/float(len(self.vevs()))
        
    def average_sub_vev(self):
        if not self.asv:
            vev1 = self.op1.average_all()
            vev2 = self.op2.average_all()
            self.asv = {t: corr - vev1*vev2 for t,corr in self.average_over_configs().iteritems()}
        return self.asv

    def jackknife_average_sub_vev(self):
        if not self.jkasv:
            
            jkvev1 = self.op1.jackknifed_full_average()
            jkvev2 = self.op2.jackknifed_full_average()

            corrjk = self.jackknifed_averages()
            jk = configtimeobj.Cfgtimeobj.fromDataDict(self.jackknifed_averages(),silent=True)
            self.jkasv = {c: {t: jk.get(config=c,time=t) - jkvev1[c]*jkvev2[c] for t in self.times} for c in self.configs}
        return self.jkasv

    def jackknifed_errors(self):
        jk = configtimeobj.Cfgtimeobj.fromDataDict(self.jackknife_average_sub_vev(),silent=True)
        asv=self.average_sub_vev()
        return {t: jackknife.errorbars(asv[t],jk.get(time=t)) for t in self.times }


    def effective_mass(self,dt):
        asv=self.average_sub_vev()
        emass = {}
        for t in self.times[:-dt]:
            try:
                emass[t] = (1.0 / float(dt)) * math.log(asv[t]/asv[t+dt])
            except ValueError:
                emass[t] = 0.0
            except KeyError:
                print "out of range"
        return emass

    def effective_mass_errors(self,dt):

        jkasv = self.jackknife_average_sub_vev()
        jkemass = {}
        for cfg in self.configs:
            asvc = jkasv[cfg]
            emass = {}
            for t in self.times[:-dt]:
                try:
                    emass[t] = (1.0 / float(dt)) * math.log(asvc[t]/asvc[t+dt])
                except ValueError:
                    emass[t] = 0.0
                except KeyError:
                    print "out of range"
            jkemass[cfg] = emass
        jkemassobj = configtimeobj.Cfgtimeobj.fromDataDict(jkemass,silent=True)
        effmass_dt = self.effective_mass(dt)
        return {t: jackknife.errorbars(effmass_dt[t],jkemassobj.get(time=t)) for t in self.times[:-dt] }
        
        
        
    # def jackknifed_errors(self):
    #     jk = configtimeobj.Cfgtimeobj.fromDataDict(self.jackknifed_averages())
    #     #print jk.get(time=0)
    #     print self.average_over_configs()
    #     return {t: jackknife.errorbars(self.average_over_configs()[t],jk.get(time=t)) for t in self.times }
        


