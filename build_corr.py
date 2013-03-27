import read_config_time_file as read
import eigenvalues as eigen
import correlator
import configtimeobj
import logging
from itertools import product


def corr_and_vev_from_files(corrfile, srcvevfile, snkvevfile):

    corrdata = read.read_config_time_data_real(corrfile).data
    vevdata_src = read.read_config_vev(srcvevfile)
    vevdata_snk = read.read_config_vev(snkvevfile)
    return correlator.Correlator.fromDataDicts(corrdata, vevdata_src, vevdata_snk)


def from_opfiles(src_opfile, snk_opfile, N=None):
    srcdata = read.read_config_time_data_real(src_opfile)
    snkdata = read.read_config_time_data_real(snk_opfile)

    if N is None:
        return correlator.Correlator.fromOpvalCTO(srcdata, snkdata)
    else:
        return correlator.Correlator.fromOpvalCTO(srcdata, snkdata, dts=list(range(N)))


def diag_from_opfiles(opfile, N=8):
    opdata = read.read_config_time_data_real(opfile)
    if N is None:
        return correlator.Correlator.fromOpvalCTO(opdata, opdata)
    else:
        return correlator.Correlator.fromOpvalCTO(opdata, opdata, dts=list(range(N)))


def from_eigenvalue_24cubed_opfiles(dirname):
    rawdata = eigen.readfile_neigenvalues(dirname, 112)
    data = configtimeobj.Cfgtimeobj.fromDataDict(eigen.reduce_to_weighted_trace(rawdata))
    data.writeeachconfig(dirname[:-1]+"_weighted/weighted")
    return correlator.Correlator.fromOpvalCTO(data, data)


def from_eigenvalue_32cubed_opfiles(dirname):
    rawdata = eigen.read_configdir_timeorlevel_evalues(dirname, 264, recall=False, store=False)
    data = configtimeobj.Cfgtimeobj.fromDataDict(eigen.reduce_to_trace(rawdata))
    return correlator.Correlator.fromOpvalCTO(data, data)


def matrix_from_opfiles(opfile_list):
    logging.debug("building matrix of correlators using %s", str(opfile_list))
    datas = [read.read_config_time_data_real(op) for op in opfile_list]

    matrix = {}
    # Product give all possible combinations
    for e1, e2 in product(enumerate(datas, start=1), repeat=2):
        index1, data1 = e1
        index2, data2 = e2
        matrix[(index1, index2)] = correlator.Correlator.fromOpvalCTO(data1, data2)

    return matrix


def matrix_from_cor_and_vev(directory, cortemplate, vevtemplate, operator_list):
    logging.debug("building matrix of correlators using %s", str(operator_list))
    matrix = {}
    for e1, e2 in product(enumerate(operator_list, start=1), repeat=2):
        index1, op1 = e1
        index2, op2 = e2
        corfile = (directory + cortemplate).format(op1, op2)
        srcvevfile = (directory + vevtemplate).format(op1)
        snkvevfile = (directory + vevtemplate).format(op2)
        matrix[(index1, index2)] = corr_and_vev_from_files(corfile, srcvevfile, snkvevfile)

    return matrix
