#!/usr/bin/env python
import re
import logging


def translate(name):
    logging.debug("Translating {}".format(name))
    new = translate_partclenames(name)
    new = translate_momentum(new)
    new = translate_irrep(new)
    new = prune_isotriplet(new)
    return new


def translate_partclenames(name):
    new = name.replace("pion","$\pi$")
    new = new.replace("rho","$\rho$")
    new = new.replace("eta","$\eta$")
    new = new.replace("phi","$\phi$")
    new = new.replace("kaon","$K$")
    new = new.replace("kbar","${K}$")
    return new


def translate_momentum(name):
    new = name.replace("p000","_AR_")
    new = new.replace("p001","_OA_")
    new = new.replace("p00-1","_OA_")
    new = new.replace("p010","_OA_")
    new = new.replace("p0-10","_OA_")
    new = new.replace("p100","_OA_")
    new = new.replace("p-100","_OA_")
    new = new.replace("p011","_PD_")
    new = new.replace("p101","_PD_")
    new = new.replace("p110","_PD_")
    new = new.replace("p0-1-1","_PD_")
    new = new.replace("p-10-1","_PD_")
    new = new.replace("p-1-10","_PD_")
    new = new.replace("p111","_CD_")
    new = new.replace("p-1-1-1","_CD_")

    new = new.replace("p002","_OA(2)_")
    new = new.replace("p00-2","_OA(2)_")
    new = new.replace("p020","_OA(2)_")
    new = new.replace("p0-20","_OA(2)_")
    new = new.replace("p200","_OA(2)_")
    new = new.replace("p-200","_OA(2)_")
    new = new.replace("p022","_PD(2)_")
    new = new.replace("p202","_PD(2)_")
    new = new.replace("p220","_PD(2)_")
    new = new.replace("p0-2-2","_PD(2)_")
    new = new.replace("p-20-2","_PD(2)_")
    new = new.replace("p-2-20","_PD(2)_")
    new = new.replace("p222","_CD(2)_")
    new = new.replace("p-2-2-2","_CD(2)_")
    return new


def translate_irrep(name):
    new = name
    for rep in ["A1", "A2", "E", "T1", "T2", "B1", "B2"]:
        baserep = rep[0]
        subrep = rep[1] if len(rep) > 1 else ""
        for gparity in ["g", "u", ""]:
            for p, newp in [("p", "+"), ("m", "-"), ("", "")]:
                if "{}{}{}".format(rep, gparity, p) == "E":
                    continue
                if "{}{}{}".format(rep, gparity, p) in new:
                    new = new.replace("{}{}{}".format(rep, gparity, p), "${}_{{{}{}}}^{{{}}}$".format(baserep, subrep, gparity, newp))

    return new


def prune_isotriplet(name):
    isomatch = "iso.*?\-"
    return re.sub(isomatch, "", name)

if __name__ == "__main__":
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

    # print translate("pionp000TDO15T1up1")
    # print translate("pionp000TDO15A2m")
    # print translate("isotripletT1up1-etap100LSD3A2m-pionp-100SS1A2m")
    # print translate("isotripletT1up1-pionp100SS1Ep-pionp-100SS1Ep")
    print translate("isotripletT1up1-etap101SS0A2p-pionp-10-1SS1B1p")
