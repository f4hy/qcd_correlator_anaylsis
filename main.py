#!/usr/bin/env python

#import tools
import plot
#import readinput
#import indexes
import configtimeobj
import read_config_time_file
import correlator

#import cProfile

emass_dts = range(2,5)

srcop = "a1pp_0_optype0_op1"

snkop = "a1pp_0_optype0_op1"


name = "%s,%s"%(srcop,snkop)

folder = "/home/bfahy/data/michael/"

srcdata = read_config_time_file.read_config_time_data_real("%sopvals_%s.orthog_test1.dat"%(folder,srcop))

#snkdata = read_config_time_file.read_config_time_data_real("%sopvals_%s.orthog_test1.dat"%(folder,snkop))
snkdata = srcdata

#srcdata = readinput.reduce_to_first_eig_val(rawsnkdata.data)
#snkdata = readinput.reduce_to_first_eig_val(rawsrcdata.data)

corr = correlator.Correlator.fromOpvalCTO(srcdata,snkdata,dts=range(0,8))

# corr.writefullcorrfile("/home/bfahy/plots/testing/single")
# corr.writeeachconfig("/home/bfahy/plots/testing/correach")

#print corr.average_over_configs()
# corr.subvev()
# print "subtracting vevs"
# print corr.average_over_configs()
# print "vevs"
# print corr.vevs()[1]
# print corr.average_vev()
# print corr.op1.average_vev()
# print corr.op2.average_vev()

# print {t: c - corr.op1.average_vev()*corr.op2.average_vev() for t,c in corr.average_over_configs().iteritems()}

# print corr.average_sub_vev()
# print "JK"
# print corr.jackknife_average_sub_vev()
# print corr.jackknifed_errors()
#exit()

# cProfile.run('correlator.Correlator.fromOpvalCTO(srcdata,snkdata)',"profile_build.out")
# cProfile.run('corr.average_sub_vev()',"profile_avg_sub_vev.out")
# cProfile.run('corr.jackknifed_errors()',"profile_jack.out")
# exit()

avgcorr = corr.average_sub_vev()
corr_errors = corr.jackknifed_errors()

plot_corr = {"%s, \t error" %name : (avgcorr, corr_errors)}

prefix = "testing"
plot.plotwitherrorbarsnames("../plots/%s/corrilator" % prefix, plot_corr, avgcorr.keys() ,autoscale=True)

for dt in emass_dts:
    emass = corr.effective_mass(dt)
    emass_errors = corr.effective_mass_errors(dt)
    plot_emass = {"emass dt=%d, \t error" %dt : (emass, emass_errors)}
    plot.plotwitherrorbarsnames("../plots/%s/emass%d" % (prefix,dt), plot_emass, emass.keys() ,autoscale=True)


#cProfile.run('corr.average_sub_vev()',"profile_avg_sub_vev.out")




