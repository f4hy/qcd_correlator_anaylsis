#!/usr/bin/env python


import correlator
import struct
import build_corr
import logging
import numpy as np
from main import plot_corr
from itertools import combinations

def diagonalize_using_coefs(coef_matrix, corr_matrix):
    "Takes correlators and coefficents and returns a correlator which is a linear combination using those coeffiecents."

    assert coef_matrix.shape == corr_matrix.shape
    print coef_matrix.shape, corr_matrix.shape

    for cor1,cor2 in combinations(corr_matrix.flat,2):
        assert  cor1.compatible(cor2)

    first_elm = corr_matrix.flat[0]
    configs = first_elm.configs
    times =  first_elm.times
        
    diaged = {}
    v = {}
    # for cfg in configs:
    #     matrix_of_cors[cfg] = {t: np.matrix(np.zeros((2,2))) for t in times}

    for cfg in configs:
        v[cfg] = np.matrix(np.array([c.vev1[cfg] for c in corr_matrix.diagonal().flat]))*np.transpose(coef_matrix)
        diaged[cfg] = [coef_matrix*np.matrix(np.array([c.get(config=cfg,time=t) for c in corr_matrix.flat]).reshape(2,2))*np.transpose(coef_matrix) for t in times]
        # for t in times:
        #     m= np.matrix(np.array([c.get(config=cfg,time=t) for c in corr_matrix.flat]).reshape(2,2))
        #     for x in (coef_matrix*m*np.transpose(coef_matrix)).flat:
        #         print x
        #     exit(0)
        #print diaged
        #exit(0)

    size = coef_matrix.shape[0]
    print coef_matrix
    print size

    # for i in range(size):
    #     for j in range(size):
    #         print (i,j)
    #         for cfg in configs:
    #             print "cfg",cfg
    #             for t in times:
    #                 print "{}, {!r}".format(t,diaged[cfg][t][i,j])

    print type(v)
    print type(diaged)
    print type(diaged.keys()[0])
    #print diaged
    
    #return diaged,v

    print diaged[configs[0]][times[0]]
    print diaged[configs[1]][times[1]]
    
    zero_zero = {cfg: {t: diaged[cfg][t][0,0] for t in times} for cfg in configs}
    vev_zero = {cfg: v[cfg][0,0] for cfg in configs}

    one_one = {cfg: {t: diaged[cfg][t][1,1] for t in times} for cfg in configs}
    vev_one = {cfg: v[cfg][0,1] for cfg in configs}
    
    
    print zero_zero[configs[0]][times[0]]
    print zero_zero[configs[1]][times[1]]

    print v[configs[0]]
    print v[configs[1]]
    print type(vev_zero[configs[0]])
    print vev_zero[configs[0]]
    print vev_zero[configs[1]]

    print vev_one[configs[0]]
    print vev_one[configs[1]]
        
    corr0 = correlator.Correlator.fromDataDicts(zero_zero,vev_zero,vev_zero)
    corr1 = correlator.Correlator.fromDataDicts(one_one,vev_one,vev_one)

    plot_corr(corr0, "/tmp/diag/", "diagzero")
    plot_corr(corr1, "/tmp/diag/", "diagone")
    
    #exit(0)



    
    
def get_diagonalize_coefs(corr_matrix,t,t0):
    first_elm = corr_matrix.flat[0]
    configs = first_elm.configs
    times =  first_elm.times
    vevmat = np.matrix(np.array([c.vev1.average()*c.vev2.average() for c in corr_matrix.flat]).reshape(2,2))
    A = np.matrix(np.array([c.average_over_configs()[t0] for c in corr_matrix.flat]).reshape(2,2))
    B = np.matrix(np.array([c.average_over_configs()[t] for c in corr_matrix.flat]).reshape(2,2))
    np.set_printoptions(precision=20)
    #print [c.sum_over_configs() for c in corr_matrix.flat]
    print "Vev",[repr(v) for v in vevmat.flat]
    print "a",[repr(a) for a in A.flat]
    print "b",[repr(b) for b in B.flat]
    print "av",[repr(a) for a in (A-vevmat).flat]
    print "bv",[repr(b) for b in (B-vevmat).flat]
    Av = A-vevmat
    Bv = B-vevmat
    invBv = np.linalg.inv(Bv)
    # invB = np.linalg.inv(B)
    # print "inv,b",[repr(b) for b in inv(B).flat]
    # #print np.linalg.inv(B)
    
    #print np.linalg.eigh(invBv*Av)
    evals,evecs = np.linalg.eig(invBv*Av)
    print evecs
    print evecs.flat[1]*22.6922
    print evecs.flat[3]*22.6922
    print evals
    return evecs
    #print evecs[0][0]*22.6922, evecs[0][1]*22.6922
    

def build_matrix_of_cors(corr_matrix):
    """Takes a matrix of corelator objects and returns an
    cfg_time_dict whose elements is a matrix of values

    """
    for cor1,cor2 in combinations(corr_matrix.flat,2):
        assert  cor1.compatible(cor2)

    first_elm = corr_matrix.flat[0]
    configs = first_elm.configs
    times =  first_elm.times
    
    matrix_of_cors = {}
    for cfg in configs:
        matrix_of_cors[cfg] = [np.matrix(np.array(
            [c.get(config=cfg,time=t)
             for c in corr_matrix.flat]).reshape(2,2)) for t in times]
                               

def split_matrix_of_cors(corr_matrix, vevs):
    """ Takes a dictionar whose elements are """
    pass
        
if __name__ == "__main__":
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

    exampledir = "/home/bfahy/r3/effectivemasses/meyer_binned/total/"
    corfile11 = exampledir +  "binned_500_a1pp_0_optype0_op1_a1pp_0_optype0_op1.cor"
    vev1file = exampledir + "binned_500_a1pp_0_optype0_op1_a1pp_0_optype0_op1.vev1"
    vev2file = exampledir + "binned_500_a1pp_0_optype0_op1_a1pp_0_optype0_op1.vev2"
    cor11 = build_corr.corr_and_vev_from_files(corfile11, vev1file, vev2file)

    corfile12 = exampledir + "binned_500_a1pp_0_optype0_op1_a1pp_0_optype10_op1.cor"
    vev1file = exampledir + "binned_500_a1pp_0_optype0_op1_a1pp_0_optype10_op1.vev1"
    vev2file = exampledir + "binned_500_a1pp_0_optype0_op1_a1pp_0_optype10_op1.vev2"
    cor12 = build_corr.corr_and_vev_from_files(corfile12, vev1file, vev2file)
    
    # corfile21 = exampledir + "binned_500_a1pp_0_optype10_op1_a1pp_0_optype0_op1.cor"
    # vev1file = exampledir + "binned_500_a1pp_0_optype10_op1_a1pp_0_optype0_op1.vev1"
    # vev2file = exampledir + "binned_500_a1pp_0_optype10_op1_a1pp_0_optype0_op1.vev2"
    # cor21 = build_corr.corr_and_vev_from_files(corfile21, vev1file, vev2file)

    cor21 = cor12

    corfile22 = exampledir + "binned_500_a1pp_0_optype10_op1_a1pp_0_optype10_op1.cor"
    vev1file = exampledir + "binned_500_a1pp_0_optype10_op1_a1pp_0_optype10_op1.vev1"
    vev2file = exampledir + "binned_500_a1pp_0_optype10_op1_a1pp_0_optype10_op1.vev2"
    cor22 = build_corr.corr_and_vev_from_files(corfile22, vev1file, vev2file)

    #sum_correlators_with_coeffs(zip([1.0,1.0,1.0],[cor11,cor12,cor22]))
    #sum_correlators_with_coeffs(zip([-7.708249,0.587736],[cor21,cor22]))
    #sum_correlators_with_coeffs(zip([22.427507,-1.557686],[cor21,cor22]))

    a = -7.708249
    b = 0.587736
    c = 22.427507
    d = -1.557686

    a = 22.63585114194867742
    b = -1.59877050648128449

    c = -4.74449486062705006
    d = 0.28418360796717246

    print np.matrix([[a,b],[c,d]])
    print np.matrix([[cor11,cor12],[cor21,cor22]])
    coefs = get_diagonalize_coefs(np.matrix([[cor11,cor12],[cor21,cor22]]),1,2)
    print coefs
    print np.matrix(coefs)
    #diagonalize_using_coefs(np.matrix(coefs),np.matrix([[cor11,cor12],[cor21,cor22]]))
    #exit(0)

    
    
    diagonalize_using_coefs(np.matrix([[a,b],[c,d]]),np.matrix([[cor11,cor12],[cor21,cor22]]))
    #sum_correlators_with_coeffs(zip([1.0],[cor11]))
