import numpy as np
import math
#DEBUG = False
DEBUG = True


class Cfgtimeobj(object):
    """ A class to handle all the objects that are indexed as
    x[config][time] and provide some methods to simplify accessing the
    data in nonstandard ways, such as a time slice of all configs.
    """
    data = None
    sums = None
    average = None

    def __init__(self, datadict, silent=False):
        self.data = datadict

        self.configs = datadict.keys()
        self.times = datadict[self.configs[0]].keys()
        self.numconfigs = len(self.configs)
        self.numtimes = len(self.times)
        if not silent:
            print "numconfigs", self.numconfigs, "numtimes", self.numtimes

        #print self.indexes()

        dataitem = self[self.configs[0]][self.times[0]]

        self.datatype = type(dataitem)
        self.scalar = np.isscalar(dataitem)

        if DEBUG:
            self.verify()

    @classmethod
    def fromDataDict(cto, datadict, silent=False):
        return cto(datadict, silent)

    @classmethod
    def fromListTuple(cto, listtuple):
        configs = list(range(len(listtuple)))

        data = {}
        for cfg in configs:
            rawcfgdata = listtuple[cfg]
            d = {}
            for rawtimedata in rawcfgdata:
                t = rawtimedata[0]  # First element is the time
                arraydata = np.asarray(tuple(rawtimedata)[1:])
                d[t] = arraydata  # rest is data e.g. real,imag
            data[cfg] = d
        return cto(data)

    def verify(self):

        if not DEBUG:
            return True

        if not self.data:
            raise ValueError("data obejct empty or false")

        sizes = map(len, self.data.values())
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

    def __setitem__(self, key, value):
        self.data[key] = value

    def get(self, config=None, time=None):

        if config is None and time is None:
            return self.data
        if time is None:
            return self.data[config]
        if config is None:
            return {cfg: v[time]   for cfg, v in self.data.iteritems()}
        else:
            return self.data[config][time]

    def compatible(self, otherobj):

        if self.configs == otherobj.configs and self.times == otherobj.times:
            return True
        else:
            return False

    def indexes(self):
        return (self.configs, self.times)

    def average_over_times(self):
        return {config: math.fsum(self.get(config=config).values()) / float(self.numtimes)
                for config in self.configs}

    def average_all(self):
        if not self.average:
            aot = self.average_over_times()
            self.average = math.fsum(aot.values()) / float(self.numconfigs)
        return self.average

    def average_over_configs(self):
        sums = self.sum_over_configs()
        return {time: sums[time] / float(self.numconfigs) for time in self.times}

    def sum_over_configs(self):
        if not self.sums:
            self.sums = {t: math.fsum(self.get(time=t).values())
                         for t in self.times}
        return self.sums

    def jackknifed_averages(self):
        sums = self.sum_over_configs()
        return {cfg: {t: (sums[t] - self.get(time=t, config=cfg)) / (self.numconfigs - 1)
                      for t in self.times} for cfg in self.configs}

    def jackknifed_full_average(self):
        total = self.average_all()
        N = float(self.numconfigs)
        Njk = float(self.numconfigs - 1)
        return {cfg: (N * total - single) / Njk
                for cfg, single in self.average_over_times().iteritems()}

    def writefullfile(self, filename):
        outfile = open(filename, 'w')
        for config in self.configs:
            for time in self.times:
                outfile.write("%d   %f\n" % (time, self.get(config=config, time=time)))
                outfile.close()

    def writeeachconfig(self, filename):
        for config in self.configs:
            outfile = open(filename + '.' + str(config), 'w')
            for time in self.times:
                outfile.write("%d   %f\n" % (time, self.get(config=config, time=time)))
                outfile.close()
