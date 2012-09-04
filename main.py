#!/usr/bin/env python

import plot
import configtimeobj
import read_config_time_file as read
import correlator


emass_dts = range(1,4)

# srcop = "a1pp_0_optype0_op1"
# snkop = "a1pp_0_optype0_op1"

# srcop = "a1mp_0_optype6_op1"
# snkop = "a1mp_0_optype6_op1"


# if snkop == srcop:
#     name = "%s"%(srcop)
# else:
#     name = "%s.%s"%(srcop,snkop)

#ops = ["a1pp_0_optype0_op1","a1pp_0_optype10_op1","a1pp_0_optype1_op1"]
ops = ["a1pp_0_optype0_op1","a1pp_0_optype10_op1"]


data_folder = "/home/bfahy/r3/results/rngtest2_smear0.17/data/"

out_folder ="/home/bfahy/glueballstesting/rngtest2/"




def compute_and_plot(srcop,snkop,data_folder,out_folder,name):


    srcdata = read.read_config_time_data_real("%sopvals_%stest1.dat"%(data_folder,srcop))
    if srcop==snkop:
        snkdata = srcdata
    else:
        snkdata = read.read_config_time_data_real("%sopvals_%stest1.dat"%(data_folder,snkop))

    corr = correlator.Correlator.fromOpvalCTO(srcdata,snkdata,dts=list(range(9)))

    # corr.writefullfile("/home/bfahy/plots/testing/corr/%s" % name)
    
    avgcorr = corr.average_sub_vev()
    corr_errors = corr.jackknifed_errors()
    
    plot_corr = {"%s, \t error" %name : (avgcorr, corr_errors)}

    plot.plotwitherrorbarsnames("%scorrelator%s" % (out_folder,name),
                                plot_corr, avgcorr.keys() ,autoscale=True)

    for dt in emass_dts:
        emass = corr.effective_mass(dt)
        emass_errors = corr.effective_mass_errors(dt)
        plot_emass = {"%s emass dt=%d, \t error" %(name,dt) : ( emass, emass_errors)}
        plot.plotwitherrorbarsnames("%semass%d.%s" % (out_folder,dt,name), plot_emass,
                                    emass.keys() ,autoscale=True)

    print "done with %s %s to %s" % (srcop,snkop,out_folder)


for op in ops:
    name = op
    compute_and_plot(op,op,data_folder,out_folder,name)
