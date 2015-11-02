##This is basic script for using dadi
#Import libraries
import dadi
import numpy as np
import itertools as itt
import pandas as pd
import sys, os

#Prettier printing of arrays
pd.set_option('precision',4)
pd.set_option('chop_threshold', 0.001)

#Plotting libraries
import matplotlib.pyplot as plt
import seaborn as sns
#matplotlib inline

###CONVERT SNP DATA TO FREQUENCY SPECTRUM###
#DEFINE a function to select most common element in a list
def most_common(L):
    return max(itt.groupby(sorted(L)),
        key=lambda(x, v):(len(list(v)),
        -L.index(x)))[0]

## a function to get alleles from ambiguous bases
def unstruct(amb):
    " returns bases from ambiguity code"
    D = {"R":["G","A"],
         "K":["G","T"],
         "S":["G","C"],
         "Y":["T","C"],
         "W":["T","A"],
         "M":["C","A"]}
    if amb in D:
        return D.get(amb)
    else:
        return [amb,amb]

##parse the loci
locifile = open("/home/mwinston/GBS/ALL_LANES/outfiles/all_lanes_clusters_GBS_10000L.loci")
loci = locifile.read().strip().split("|")

## define focal populations and the outgroups
burchellii_A        = ['sample227','sample236','sample269','sample270','sample275','sample276','sample277','sample278','sample280','sample297','sample321','sample339','sample344','sample345','sample357','sample358','sample360','sample361','sample364','sample365','sample366','sample369','sample370','sample374','sample375','sample376','sample377','sample378','sample413','sample422','sample424','sample426','sample427','sample428','sample430','sample431','sample432','sample54']
burchellii_B    = ['sample415','sample224','sample226','sample52','sample55','sample56']
burchellii_C  = ['sample115','sample127','sample162','sample165','sample168','sample175','sample182','sample19','sample194','sample198','sample209','sample211','sample213','sample244','sample261','sample284','sample326','sample327','sample330','sample34','sample37','sample387','sample394','sample400','sample416','sample418','sample60','sample63','sample77','sample79','sample84','sample86','sample97','sample99']
outgroup    = ['sample26','sample28','sample31','sample38','sample39','sample49','sample72']

## BC is a composite of two 'species' (B + C)
burchellii_BC = burchellii_B + burchellii_C

###PROJECTION OF DATA

## minimum samples required to use the locus
proj = [10,6,10]

## only examine loci w/ at least 1 outgroup sample
## and at least two inds from each focal populations
Floci = []   ## filtered locus list
for loc in loci:
    names = [i.strip().split(" ")[0] for i in loc.strip().split("\n")[:-1]]
    if len(set([">"+i for i in outgroup]).intersection(names)) >= 1:
        if len(set([">"+i for i in burchellii_A]).intersection(names)) >= proj[0]/2:
            if len(set([">"+i for i in burchellii_B]).intersection(names)) >= proj[1]/2:
                if len(set([">"+i for i in burchellii_C]).intersection(names)) >= proj[2]/2:
                    Floci.append(loc)

##
## outfile location
outfile = open("burch.dadi.snps", 'w')

## print header
print >>outfile, "\t".join(["REF","OUT",
                            "Allele1","burchellii_A","burchellii_B","burchellii_C",
                            "Allele2","burchellii_A","burchellii_B","burchellii_C",
                            "locus","site"])

## recording values
uncalled = 0
calledSNP = 0
uncalledSNP = 0

## iterate over loci
for locN, loc in enumerate(Floci):
    ## to filter a single SNP per locus
    singleSNP = 0

    ## separate names, loci, and SNP identifier line
    names = [i.strip().split(" ")[0] for i in loc.strip().split("\n")[:-1]]
    dat = [i.strip().split(" ")[-1] for i in loc.strip().split("\n")[:-1]]
    snps = loc.strip().split("\n")[-1][14:]

    ## select ingroups v. outgroups
    ihits = [names.index(">"+i) for i in burchellii_A+burchellii_B+burchellii_C if ">"+i in names]
    ohits = [names.index(">"+i) for i in outgroup if ">"+i in names]

    ## select a single variable SNP from each locus
    ## but also record how many are skipped for calculating
    ## effective sequence length later
    for site,char in enumerate(snps):
        if char in ["-","*"]:
            ## get common (reference) alleles
            i2 = ["-"+dat[k][site]+"-" for k in ihits if dat[k][site] not in ["N","-"]]
            o2 = ["-"+dat[k][site]+"-" for k in ohits if dat[k][site] not in ["N","-"]]
            ## filter out uninformative sites
            b1 = any([i[1] not in ['N','-'] for i in i2])
            b2 = any([i[1] not in ['N','-'] for i in o2])
            if not (b1 and b2):
                uncalled += 1
            else:
                ## if site is variable in the ingroups
                if len(set(i2)) > 1:
                    ## get the segregating alleles
                    alleles = list(itt.chain(*[unstruct(i[1]) for i in i2]))
                    allele1 = most_common(alleles)
                    allele2 = most_common([i for i in alleles if i!=allele1])
                    outg = most_common(list(itt.chain(*[unstruct(i[1]) for i in o2])))
                    ## burchellii_A
                    bA = [names.index(">"+z) for z in burchellii_A if ">"+z in names]
                    bA_dat = list(itt.chain(*[unstruct(dat[z][site]) for z in bA]))
                    bA1,bA2 = bA_dat.count(allele1), bA_dat.count(allele2)
                    ## burchellii_B
                    bB = [names.index(">"+z) for z in burchellii_B if ">"+z in names]
                    bB_dat = list(itt.chain(*[unstruct(dat[z][site]) for z in bB]))
                    bB1,bB2 = bB_dat.count(allele1), bB_dat.count(allele2)
                    ## burchellii_C
                    bC = [names.index(">"+z) for z in burchellii_C if ">"+z in names]
                    bC_dat = list(itt.chain(*[unstruct(dat[z][site]) for z in bC]))
                    bC1,bC2 = bC_dat.count(allele1), bC_dat.count(allele2)

                    if not singleSNP:
                        calledSNP += 1
                        print >>outfile, "\t".join(map(str,["-"+allele1+"-",
                                                            "-"+outg+"-",
                                                            allele1, bA1, bB1, bC1,
                                                            allele2, bA2, bB2, bC2,
                                                            str(locN), str(site)]))
                        singleSNP = 1
                    else:
                        uncalledSNP += 1

                        
outfile.close()

Nloci = len(Floci)
print "loci\t",Nloci
print "bp\t",Nloci*90
L = Nloci*90-uncalled
print "called bp\t",L
print 'called SNP\t', calledSNP
print 'uncalledSNP\t', uncalledSNP
propcalled = calledSNP/float(calledSNP+uncalledSNP)
print 'prop\t', propcalled                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
