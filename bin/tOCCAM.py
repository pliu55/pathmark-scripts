#!/usr/bin/env python
"""tOCCAM.py: 

Usage:
  tOCCAM.py [options] phenotype data

Options:
  -q            run quietly
"""
## Written By: Sam Ng
import getopt, math, os, re, sys
from mPATHMARK import *

verbose = True

def usage(code = 0):
    print __doc__
    if code != None:
        sys.exit(code)

def log(msg, die = False):
    if verbose:
        sys.stderr.write(msg)
    if die:
        sys.exit(1)

def syscmd(cmd):
    log("running:\n\t"+cmd+"\n")
    exitstatus = os.system(cmd)
    if exitstatus != 0:
        print "Failed with exit status %i" % (exitstatus)
        sys.exit(10)
    log("... done\n")

def ttest(values0, values1, alpha = 0.05):
    (mean0, std0) = mean_std(values0)
    (mean1, std1) = mean_std(values1)
    try:
        tval = (mean1-mean0)/(math.sqrt((std1**2)/len(values1)+(std0**2)/len(values0))+alpha)
    except ZeroDivisionError:
        tval = "NA"
    except TypeError:
        tval = "NA"
    return(tval)

def main(args):
    ## parse arguments
    try:
        opts, args = getopt.getopt(args, "q")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) != 2:
        log("ERROR: incorrect number of arguments", die = True)
    
    phenotypeFile = args[0]
    dataFile = args[1]
    
    global verbose
    for o, a in opts:
        if o == "-q":
            verbose = False
    
    ## execute
    phenotypeName = re.split("/", phenotypeFile)[-1]
    dataName = re.split("/", dataFile)[-1]
    outputDir = "OCCAM__%s__%s" % (phenotypeName, dataName)
    syscmd("mkdir %s" % (outputDir))
    (phenData, phenFeatures, phenSamples) = rCRSData(phenotypeFile, retFeatures = True)
    
    ## compute stats
    o = open("%s/results.tab" % (outputDir), "w")
    o.write("id\t%s\n" % ("\t".join(phenFeatures)))
    f = open(dataFile, "r")
    dataSamples = f.readline().rstrip().split("\t")[1:]
    sampleIndex = {}
    for i, sample in enumerate(dataSamples):
        sampleIndex[sample] = i
    for line in f:
        if line.isspace():
            continue
        dataValues = line.rstrip().split("\t")
        dataFeature = dataValues.pop(0)
        o.write("%s" % (dataFeature))
        for phenotype in phenFeatures:
            posSamples = []
            negSamples = []
            for sample in phenSamples:
                if sample not in dataSamples:
                    continue
                if phenData[phenotype][sample] == 1:
                    posSamples.append(sample)
                elif phenData[phenotype][sample] == 0:
                    negSamples.append(sample)
            val = ttest([float(dataValues[sampleIndex[i]]) for i in negSamples], [float(dataValues[sampleIndex[i]]) for i in posSamples])
            o.write("\t%s" % (val))
        o.write("\n")
    f.close()
    o.close()

if __name__ == "__main__":
    main(sys.argv[1:])
