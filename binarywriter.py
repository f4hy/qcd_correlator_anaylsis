#!/usr/bin/env python

"""This contains methods to write out a correlator in binary format
which is expected by colin's analyis code. 

"""
import correlator
import struct
import build_corr
import logging
import math
from itertools import combinations



def write_binary_correlator(filename,cor):



    f = open(filename,'wb')
    NSTRINGS = 20
    f.write(struct.pack('i',NSTRINGS)) # Nstrings
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                           beta     \0") # beta at 27?
    f.write("                           %        \0") # a_{s}/a_{t}=?
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("               16              36   \0") # Lattice
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                               %s\0" %cor.numconfigs)
    f.write("                                               1 \0")
    f.write("                                                 \0")
    f.write("                                    \0")
    f.write(struct.pack('i',1)) # NCOR
    f.write(struct.pack('i',cor.numtimes-1)) # MAX_DT
    f.write(struct.pack('i',cor.numtimes)) # T_PERIOD
    f.write(struct.pack('i',cor.numconfigs)) # NBINS
    # for x in range(NBINS*(MAX_DT+1)*(NCOR*(NCOR+1))/2):
    #     f.write(struct.pack('d',float(x)))
    for cfg in cor.configs:
        for t in cor.times:
            f.write(struct.pack('d', cor.get(config=cfg,time=t)))
    for cfg in cor.configs:
        f.write(struct.pack('d', cor.vev1[cfg]))


def write_binary_correlator_matrix(filename,cor_matrix):

    assert all(c1.compatible(c2) for c1,c2 in combinations(cor_matrix.values(),2))

    size = int(math.sqrt(len(cor_matrix)))
        
    first = cor_matrix[(1,1)]
    
    f = open(filename,'wb')
    NSTRINGS = 20
    f.write(struct.pack('i',NSTRINGS)) # Nstrings
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                           beta     \0") # beta at 27?
    f.write("                           %        \0") # a_{s}/a_{t}=?
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("               16              36   \0") # Lattice
    f.write("                                    \0")
    f.write("                                    \0")
    f.write("                                               %s\0" %first.numconfigs)
    f.write("                                               1 \0")
    f.write("                                                 \0")
    f.write("                                    \0")
    f.write(struct.pack('i',size)) # NCOR
    f.write(struct.pack('i',first.numtimes-1)) # MAX_DT
    f.write(struct.pack('i',first.numtimes)) # T_PERIOD
    f.write(struct.pack('i',first.numconfigs)) # NBINS
    # for x in range(NBINS*(MAX_DT+1)*(NCOR*(NCOR+1))/2):
    #     f.write(struct.pack('d',float(x)))
    
    for cfg in first.configs:
        for t in first.times:
            for j in range(1,size+1):
                for i in range(1,j+1):
                    f.write(struct.pack('d', cor_matrix[(i,j)].get(config=cfg,time=t)))

    # for cfg in cor_list[0].configs:
    #     f.write(struct.pack('d', (cor_list[0]).vev1[cfg]))
            
    # for cfg in cor_list[-1].configs:
    #     f.write(struct.pack('d', cor_list[-1].vev1[cfg]))
    for cfg in first.configs:
        for diag in range(1,size+1):
            #print cor_matrix[(diag,diag)].vev1[cfg]
            f.write(struct.pack('d', cor_matrix[(diag,diag)].vev1[cfg]))
        
        
if __name__ == "__main__":
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

    # exampledir = "/home/bfahy/r3/effectivemasses/meyer_binned/total/"
    # corfile11 = exampledir +  "binned_500_a1pp_0_optype0_op1_a1pp_0_optype0_op1.cor"
    # vev1file = exampledir + "binned_500_a1pp_0_optype0_op1_a1pp_0_optype0_op1.vev1"
    # vev2file = exampledir + "binned_500_a1pp_0_optype0_op1_a1pp_0_optype0_op1.vev2"
    # cor11 = build_corr.corr_and_vev_from_files(corfile11, vev1file, vev2file)

    # corfile12 = exampledir + "binned_500_a1pp_0_optype0_op1_a1pp_0_optype10_op1.cor"
    # vev1file = exampledir + "binned_500_a1pp_0_optype0_op1_a1pp_0_optype10_op1.vev1"
    # vev2file = exampledir + "binned_500_a1pp_0_optype0_op1_a1pp_0_optype10_op1.vev2"
    # cor12 = build_corr.corr_and_vev_from_files(corfile12, vev1file, vev2file)
    
    # corfile21 = exampledir + "binned_500_a1pp_0_optype10_op1_a1pp_0_optype0_op1.cor"
    # vev1file = exampledir + "binned_500_a1pp_0_optype10_op1_a1pp_0_optype0_op1.vev1"
    # vev2file = exampledir + "binned_500_a1pp_0_optype10_op1_a1pp_0_optype0_op1.vev2"
    # cor21 = build_corr.corr_and_vev_from_files(corfile21, vev1file, vev2file)

    # corfile22 = exampledir + "binned_500_a1pp_0_optype10_op1_a1pp_0_optype10_op1.cor"
    # vev1file = exampledir + "binned_500_a1pp_0_optype10_op1_a1pp_0_optype10_op1.vev1"
    # vev2file = exampledir + "binned_500_a1pp_0_optype10_op1_a1pp_0_optype10_op1.vev2"
    # cor22 = build_corr.corr_and_vev_from_files(corfile22, vev1file, vev2file)

    
    

    # exampledir = "/home/bfahy/r3/meyer_results_moreops/stream_1/data/"
    # opfile1 = exampledir +  "opvals_a1pp_0_optype0_op1_diag.dat"
    # opfile2 = exampledir +  "opvals_a1pp_0_optype5_op1_diag.dat"
    # opfile3 = exampledir +  "opvals_a1pp_0_optype11_op1_diag.dat"

    #matrix1 = build_corr.matrix_from_opfiles([opfile1,opfile2,opfile3])

    #ops = ["a1pp_0_optype{}_op1".format(i) for i in [0,1,3,5,10,11]]
    ops = ["a1pp_0_optype{}_op1".format(i) for i in [5,11]]
    # ops = ["a1pp_0_optype{}_op1".format(i) for i in [0,10]]

    directory = "/home/bfahy/r3/effectivemasses/meyer_moreops/total/"
    
    matrix = build_corr.matrix_from_cor_and_vev(directory,"binned_500_{}_{}.cor","binned_500_{}_{}.vev1","binned_500_{}_{}.vev2", ops)

    


    # directory = "/home/bfahy/r3/diag_test_data/"
    # ops = ['Single_Site_0', 'Triply_Displaced_O_3', 'Triply_Displaced_U_0', 'Triply_Displaced_U_2']
    
    # matrix = build_corr.matrix_from_cor_and_vev(directory,"isovector_du_OhD_A1um_1_000_{}-isovector_du_OhD_A1um_1_000_{}.0.conn.dat","vev.dat","vev.dat", ops)
    
    
    
    # corfile11 = exampledir +  "isoscalar_OhD_A1gp_1_000_Single_Site_0-isoscalar_OhD_A1gp_1_000_Single_Site_0.A1gp.conn.dat"
    # corfile12 = exampledir +  "isoscalar_OhD_A1gp_1_000_Single_Site_0-isoscalar_OhD_A1gp_1_000_Singly_Displaced_2.A1gp.conn.dat"
    # corfile22 = exampledir +  "isoscalar_OhD_A1gp_1_000_Singly_Displaced_2-isoscalar_OhD_A1gp_1_000_Singly_Displaced_2.A1gp.conn.dat"
    # vev1file = exampledir + "vev.dat"
    # cor11 = build_corr.corr_and_vev_from_files(corfile11, vev1file, vev1file)
    # cor12 = build_corr.corr_and_vev_from_files(corfile12, vev1file, vev1file)
    # cor22 = build_corr.corr_and_vev_from_files(corfile22, vev1file, vev1file)

    #matrix = {(1,1) : cor11, (1,2) : cor12, (2,1) : cor21, (2,2): cor22}
    #matrix = {(1,1) : cor11}

    
    #write_binary_correlator("bin500.data1pp",cor11)
    write_binary_correlator_matrix("newop.data1pp",matrix)

    
