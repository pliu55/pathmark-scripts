#!/usr/bin/env python
"""
galaxyPATHMARK.py
    by Sam Ng
"""
import logging, math, os, random, re, shutil, sys, types, zipfile
from copy import deepcopy

from optparse import OptionParser

## logger
logging.basicConfig(filename="galaxy-pathmark.log", level=logging.INFO)

## executables
bin_dir = os.path.dirname(os.path.abspath(__file__))
signature_exec = os.path.join(bin_dir, "signature.py")
pathmark_exec = os.path.join(bin_dir, "PATHMARK.py")

def main():
    ## check for fresh run
    if os.path.exists(".jobTree"):
        logging.warning("WARNING: '.jobTree' directory found, remove it first to start a fresh run\n")
    
    ## parse arguments
    parser = OptionParser(usage = "%prog [options] data_matrix phenotype_matrix pathway_file")
    parser.add_option("-f", "--filter", dest = "filter_parameters", default = "0.0;0.0",
                      help = "")
    parser.add_option("-t", "--heat", dest = "heat_diffusion", default = "0",
                      help = "")
    parser.add_option("-u", "--hub", dest = "hub_filter", action = "store_true", default = False,
                      help = "")
    parser.add_option("-o", "--output", dest = "output_file", default = None,
                      help = "")
    options, args = parser.parse_args()
    
    work_dir = os.path.abspath("./")
    
    if len(args) != 3:
        logging.error("ERROR: incorrect number of arguments\n")
        sys.exit(1)
    data_file = os.path.abspath(args[0])
    phenotype_file = os.path.abspath(args[1])
    pathway_file = os.path.abspath(args[2])
    
    ## run signature.py
    cmd = "%s %s %s" % (signature_exec, data_file, phenotype_file)
    os.system(cmd)
    
    ## run PATHMARK.py
    cmd = "%s" % (pathmark_exec)
    if os.path.exists("null_signature.tab"):
        cmd += " -n %s" % ("null_signature.tab")
    if os.path.exists("bootstrap_signature.tab"):
        cmd += " -b %s" % ("bootstrap_signature.tab")
    if options.output_file is not None:
        cmd += " -o %s" % (options.output_file)
    cmd += " -f \"%s\"" % (options.filter_parameters)
    cmd += " -t %s" % (options.heat_diffusion)
    if options.hub_filter:
        cmd += " -u"
    cmd += " signature.tab %s" % (pathway_file)
    os.system(cmd)
    
    ## prepare outputs
    if options.output_file is not None:
        zip_file = zipfile.ZipFile("report.zip", "w")
        zipDirectory("report", zip_file)
        zip_file.close()
        shutil.copy(os.path.join(work_dir, "report.zip"), options.output_file)


if __name__ == "__main__":
    main()
