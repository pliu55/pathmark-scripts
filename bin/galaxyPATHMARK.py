#!/usr/bin/env python
"""
galaxyPATHMARK.py
    by Sam Ng
"""
import logging, math, os, random, re, shutil, sys, types, zipfile
from copy import deepcopy

from optparse import OptionParser

## logger
logging.basicConfig(filename = "galaxy-pathmark.log", level = logging.INFO)

## executables
bin_dir = os.path.dirname(os.path.abspath(__file__))
signature_exec = os.path.join(bin_dir, "signature.py")
pathmark_exec = os.path.join(bin_dir, "PATHMARK.py")

## functions
def zipDirectory(directory, zip):
    for root, dirs, files in os.walk(directory):
        for file in files:
            zip.write(os.path.join(root, file))

def main():
    ## check for fresh run
    if os.path.exists(".jobTree"):
        logging.warning("WARNING: '.jobTree' directory found, remove it first to start a fresh run\n")
    
    ## parse arguments
    parser = OptionParser(usage = "%prog [options] data_matrix phenotype_matrix pathway_file")
    parser.add_option("-b", "--bootstrap", dest = "bootstrap_size", default = 0,
                      help = "")
    parser.add_option("-n", "--null", dest = "null_size", default = 0,
                      help = "")
    parser.add_option("-p", "--permute", dest = "null_permute", default = "paradigm",
                      help = "")
    parser.add_option("-m", "--method", dest = "signature_method", default = "sam",
                      help = "")
    parser.add_option("-f", "--filter", dest = "filter_parameters", default = "0.0;0.0",
                      help = "")
    parser.add_option("-t", "--heat", dest = "heat_diffusion", default = "0",
                      help = "")
    parser.add_option("-u", "--hub", dest = "hub_filter", action = "store_true", default = False,
                      help = "")
    parser.add_option("-o", "--output", dest = "output_file", default = None,
                      help = "")
    parser.add_option("-z", "--seed", dest = "seed", default = None,
                      help = "")
    options, args = parser.parse_args()
    logging.info("options: %s" % (str(options)))
    
    work_dir = os.path.abspath("./")
    
    if len(args) != 3:
        logging.error("ERROR: incorrect number of arguments\n")
        sys.exit(1)
    data_file = os.path.abspath(args[0])
    phenotype_file = os.path.abspath(args[1])
    pathway_file = os.path.abspath(args[2])
    
    ## run signature.py
    cmd = "%s %s" % (sys.executable, signature_exec)
    cmd += " -b %s" % (options.bootstrap_size)
    cmd += " -n %s" % (options.null_size)
    cmd += " -p %s" % (options.null_permute)
    cmd += " -m %s" % (options.signature_method)
    if options.seed is not None:
        cmd += " -z %s" % (options.seed)
    cmd += " %s %s" % (data_file, phenotype_file)
    os.system(cmd)
    logging.info("system: %s" % (cmd))
    
    ## run PATHMARK.py
    cmd = "%s %s" % (sys.executable, pathmark_exec)
    if os.path.exists("null_signature.tab"):
        cmd += " -n %s" % ("null_signature.tab")
    if os.path.exists("bootstrap_signature.tab"):
        cmd += " -b %s" % ("bootstrap_signature.tab")
    cmd += " -f \"%s\"" % (options.filter_parameters)
    cmd += " -t %s" % (options.heat_diffusion)
    if options.hub_filter:
        cmd += " -u"
    cmd += " signature.tab %s" % (pathway_file)
    os.system(cmd)
    logging.info("system: %s" % (cmd))
    
    ## prepare outputs
    if options.output_file is not None:
        zip_file = zipfile.ZipFile("report.zip", "w")
        zipDirectory("report", zip_file)
        zip_file.close()
        shutil.copy(os.path.join(work_dir, "report.zip"), options.output_file)

if __name__ == "__main__":
    main()
