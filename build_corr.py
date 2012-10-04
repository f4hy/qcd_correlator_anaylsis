import read_config_time_file as read
import eigenvalues
import correlator
import configtimeobj
from itertools import product

def corr_and_vev_from_files(corrfile, srcvevfile, snkvevfile):

    corrdata = read.read_config_time_data_real(corrfile).data
    vevdata_src = read.read_config_vev(srcvevfile)
    vevdata_snk = read.read_config_vev(snkvevfile)
    return correlator.Correlator.fromDataDicts(corrdata, vevdata_src, vevdata_snk)


def from_opfiles(src_opfile, snk_opfile, N=8):
    srcdata = read.read_config_time_data_real(src_opfile)
    snkdata = read.read_config_time_data_real(snk_opfile)

    return correlator.Correlator.fromOpvalCTO(srcdata, snkdata, dts=list(range(N)))


def diag_from_opfiles(opfile, N=8):
    opdata = read.read_config_time_data_real(opfile)
    return correlator.Correlator.fromOpvalCTO(opdata, opdata, dts=list(range(N)))


def from_eigenvalue_24cubed_opfiles(dirname):
    rawdata = eigenvalues.readfile_neigenvalues(dirname, 112)
    data = configtimeobj.Cfgtimeobj.fromDataDict(eigenvalues.reduce_to_trace(rawdata))
    return correlator.Correlator.fromOpvalCTO(data, data)

def from_eigenvalue_32cubed_opfiles(dirname):
    rawdata = eigenvalues.read_configdir_timeorlevel_evalues(dirname, 264, recall=False, store=False)
    data = configtimeobj.Cfgtimeobj.fromDataDict(eigenvalues.reduce_to_trace(rawdata))
    return correlator.Correlator.fromOpvalCTO(data, data)


def matrix_from_opfiles(opfile_list):
    datas = [read.read_config_time_data_real(op) for op in opfile_list]

    
    
    matrix = {}
    # Product give all possible combinations
    for e1,e2 in product(enumerate(datas, start=1),repeat=2):
        print "e1,e2"
        print e1
        print e2
        index1,data1 = e1
        index2,data2 = e2
        print index1,data1,index2,data2
        matrix[(index1,index2)] = correlator.Correlator.fromOpvalCTO(data1,data2)

    print matrix
    return matrix
