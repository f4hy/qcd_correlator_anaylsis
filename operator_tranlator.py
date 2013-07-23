#!/usr/bin/env python
import re

def translate(name):
    new = translate_partclenames(name)
    print new
    new = translate_momentum(new)
    print new
    new = translate_irrep(new)
    print new
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
    new = new.replace("p001","__OA__")
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
    return new

def translate_irrep(name):
    new = name
    for rep in ["A1", "A2", "E", "T1", "T2"]:
        baserep = rep[0]
        subrep = rep[1] if len(rep) > 1 else ""
        for gparity in ["g", "u"]:
            for p,newp in [("p","+"), ("m","-"), ("","")]:
                print rep,gparity,p
                print "{}{}{}".format(rep,gparity,p), "${}_{}^{}$".format(rep,gparity,newp)
                new = new.replace("{}{}{}".format(rep,gparity,p), "${}_{{{}{}}}^{{{}}}$".format(baserep, subrep,gparity,newp) )
    return new

def prune_isotriplet(name):
    isomatch = "isotriplet.*?\-"
    return re.sub(isomatch, "", name)
