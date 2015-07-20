# FET
#Need to make it so that for each Pfam term, make a contingency table, perform FET, output significantly over- or under- represented terms.

#Make a dictionary with genes:lists of their Pfam terms
PfamDict = {}
f=open("/home/laura/Dropbox/IRT3/PfamScanOutput/Dgal_okay_PfamScanNH.txt", 'r')
for line in f.readlines()[1:]:
    line = line.strip().split("\t")
    if line[0] not in PfamDict:
        PfamDict[line[0]] = [line[5]]
    else:
        PfamDict[line[0]].append(line[5])

f.close()
#How many genes (of all genes in proteome) have annotated domains (all)?  len(PfamDict)   16239

#Need a subset, e.g. SBG, so I could just subset the dictionary by its keys.
subsetFile="/home/laura/Dropbox/IRT3/DgalSexBiased.txt"
sf = open(subsetFile).readlines()
subsetGeneNames = []
for line in sf:
    gene = line.strip()
    subsetGeneNames.append(gene)

#len(subsetGeneNames) #5838
for g in subsetGeneNames:
    if g not in PfamDict.keys():
        subsetGeneNames.remove(g)

#len(subsetGeneNames) #4092
#Make a subset dictionary:
subsetDict={}
for k in PfamDict.keys():
    if k in subsetGeneNames:
        subsetDict[k] = PfamDict[k]

#For each Pfam term, how many genes have this term?
PfamScanResults = "/home/laura/Dropbox/IRT3/PfamScanOutput/Dgal_okay_PfamScanNH.txt"
lines=open(PfamScanResults).readlines()[1:]
PFNrs = []
PFlines = []
for line in lines:
    line = line.strip().split('\t')
    PFlines.append(line)
    geneName = line[0]
    PFNrs.append(line[5])

'''
allGeneNames = []
proteome = "/home/laura/Dropbox/IRT3/DgalAllGeneNames.txt"
l = open(proteome).readlines()
for line in l:
    gene = line.strip()
    allGeneNames.append(gene)
'''

def countOccurrencesOfPfamTerm(myPFnr, PFlines):
    count = 0
    for l in PFlines:
        if l[5] == myPFnr:
            count +=1
    return count

for p in PFNrs:
    countOccurrencesOfPfamTerm(p, PFlines)

#Now for the subset, how many genes have this term?
#Subset the PfamScan output file:
lines=open(PfamScanResults).readlines()[1:]
subsetPFlines = []
subsetPFNrs = []
for line in lines:
    line = line.strip().split('\t')
    if line[0] in subsetGeneNames:
        subsetPFlines.append(line)
        subsetPFNrs.append(line[5])

#Call the function for the subset
for p in subsetPFNrs:
    countOccurrencesOfPfamTerm(p, subsetPFlines)


#Make contingency tables like [[subset genes with domain x, subset genes],[genes with domain x, genes]]
import numpy
#How many genes have all annotated domains?      len(PfamDict)  #16239
#How many genes in the subset have all annoated domains?    len(subsetDict)   #3160
#How many genes have domain x?

def continTableOfPfamTerm(myPFnr, PFlines):
    bl = 0
    for l in PFlines:
        if l[5] == myPFnr:
            bl +=1
    tl = 0
    for l in subsetPFlines:
        if l[5] == myPFnr:
            tl += 1
    tr = len(subsetDict)
    br = len(PfamDict)
    return numpy.matrix([[tl,tr],[bl,br]])


#import scipy stats and perform FET for each term.
from scipy import stats
#example:  stats.fisher_exact(numpy.matrix([[10,3160],[43,16239]]))
FETresults={}
for p in PFNrs:
    #term,oddsratio,pvalue = p, stats.fisher_exact(continTableOfPfamTerm(p,PFlines))
    FETresults[p] =stats.fisher_exact(continTableOfPfamTerm(p,PFlines))

#output all the significant genes.
outfile= open("/home/laura/Dropbox/IRT3/SignificantPfamTerms_FET", 'w')
outfile.write(subsetFile)
outfile.write("\n")
for p in FETresults.keys():
    if FETresults[p][1] < 0.05:
        outfile.write("%s" % p)
        outfile.write("\t")
        outfile.write("%s" % FETresults[p][1])
        outfile.write("\n")

outfile.close()

for p in FETresults.keys():
    if FETresults[p][1] < 0.05:
        print p
        print FETresults[p]
