import numpy as np
import math
DEBUG = False

class Cfgtimeobj(object):
    """ A class to handle all the objects that are indexed as
    x[config][time] and provide some methods to simplify accessing the
    data in nonstandard ways, such as a time slice of all configs. 
    """
    data = None
    sums = None
    average = None
    
    def __init__(self,datadict):
        #self.numberofvalues = len(self.data[self.configs[0]][self.times[0]])
        self.data = datadict

        self.configs = datadict.keys()
        self.times = datadict[self.configs[0]].keys()
        self.numconfigs = len(self.configs)
        self.numtimes = len(self.times)
        #print "configs",self.configs,"times",self.times
        print "numconfigs",self.numconfigs,"numtimes",self.numtimes

        #print self.indexes()
        
        dataitem = self[self.configs[0]][self.times[0]]
        
        self.datatype =  type(dataitem)
        self.scalar = np.isscalar(dataitem)
        
        self.verify()
        

        
    @classmethod
    def fromDataDict(cto,datadict):
        return cto(datadict)

    @classmethod
    def fromListTuple(cto,listtuple):
        configs = list(range(len(listtuple)))
        #print configs
        #cto.configs= configs
        #cto.times=list(range(len(listtuple[0])))

        data = {}
        for cfg in configs:
            rawcfgdata = listtuple[cfg]
            d = {}
            for rawtimedata in rawcfgdata:
                t = rawtimedata[0] # First element is the time
                # print type(rawtimedata)
                # print np.array(rawtimedata)
                arraydata = np.asarray(tuple(rawtimedata)[1:])
                #data = list(rawtimedata)[1]
                d[t]=arraydata # rest is data e.g. real,imag
            data[cfg]=d
        # print data
        # print blah
        #cto.data=data
        return cto(data)
        

    def verify(self):

        if not DEBUG:
            print("not verified")
            return True
        
        if not self.data:
            raise ValueError("data obejct empty or false")
        
        sizes = map(len,self.data.values())
        if (sizes.count(sizes[0]) != len(sizes)):
            raise ValueError("Object size is inconsistant")

        for cfg in self.configs:
            for time in self.times:
                if type(self.data[cfg][time]) != self.datatype:
                    print type(self.data[cfg][time])
                    raise TypeError("Not all data is the same type")
                if self.data[cfg][time] == None:
                    raise ValueError("indexed value is none")

        print("Cfg Time Object verified for consistancy")
        return True
        
    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key,value):
        self.data[key] = value

    
    def get(self, config = None, time = None):

        if config is None and time is None:
            return self.data
        if time is None:
            return self.data[config]
        if config is None:
            #timeslice = dict()
            return { cfg:v[time]   for cfg, v in self.data.iteritems() }
        else:
            return self.data[config][time]
        
        
    def compatible(self, otherobj):

        if self.configs == otherobj.configs and self.times == otherobj.times:
            return True
        else:
            return False

    
        
    def indexes(self):
        return (self.configs, self.times)
        
    # def vevs(self):
    #     self.vevdata = {config: sum(self.get(config=config).values())/float(self.numtimes) for config in self.configs}
    #     return self.vevdata

    # def average_vev(self):
    #     #print self.vevs()
    #     return float(sum(self.vevs().values()))/float(len(self.vevs()))

    def average_over_times(self):
        return {config: math.fsum(self.get(config=config).values())/float(self.numtimes) for config in self.configs}

    def average_all(self):
        if not self.average:
            aot = self.average_over_times()
            self.average = math.fsum(aot.values())/float(self.numconfigs)
        return self.average

        
    def average_over_configs(self):
        sums = self.sum_over_configs();
        return {time: sums[time]/float(self.numconfigs) for time in self.times}

    def sum_over_configs(self):
        if not self.sums:
            self.sums = {t: math.fsum(self.get(time=t).values()) for t in self.times}
        return self.sums
            
    def jackknifed_averages(self):
        sums = self.sum_over_configs()
        return {cfg: {t: (sums[t]-self.get(time=t,config=cfg))/(self.numconfigs-1) for t in self.times} for cfg in self.configs}
        
    def jackknifed_full_average(self):
        total = self.average_all()
        N = float(self.numconfigs)
        Njk = float(self.numconfigs-1)
        return {cfg: (N*total - single)/Njk for cfg,single in self.average_over_times().iteritems()}
        
    def writefullfile(self,filename):
        outfile = open(filename,'w')
        for config in self.configs:
            for time in self.times:
                outfile.write("%d   %f\n" %(time,self.get(config=config,time=time)))
                outfile.close()

    def writeeachconfig(self,filename):
        for config in self.configs:
            outfile = open(filename+'.'+str(config),'w')
            for time in self.times:
                outfile.write("%d   %f\n" %(time,self.get(config=config,time=time)))
                outfile.close()
