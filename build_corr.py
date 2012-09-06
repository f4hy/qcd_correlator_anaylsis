import read_config_time_file as read
import correlator


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