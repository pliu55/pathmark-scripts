#!/usr/bin/env python
"""
galaxyPATHMARK.py

Author: Sam Ng
Last Updated: 2014-06-13
"""
import math, os, random, re, shutil, sys, types
from copy import deepcopy
from optparse import OptionParser

## default variables
base_directory = os.path.dirname(os.path.abspath(__file__))
signature_executable = 'signature.py'
pathmark_executable = 'PATHMARK.py'

def main():
    ## check for fresh run
    if os.path.exists('.jobTree'):
        logger('WARNING: .jobTree directory found, remove it first to start a fresh run\n', die = True)
    
    ## parse arguments
    parser = OptionParser(usage = '%prog [options]')
    parser.add_option('-d', '--data', dest='data_file', default=None)
    parser.add_option('-p', '--phenotype', dest='phenotype_file', default=None)
    parser.add_option('-n', '--pathway', dest='pathway_file', default=None)
    parser.add_option('-o', '--output', dest='output_file', default=None)
    parser.add_option('-f', '--filter', dest='filter_parameters', default='0.0;0.0')
    parser.add_option('-t', '--heat', dest='heat_diffusion', default='0')
    parser.add_option('-u', '--hub', dest="hub_filter", action='store_true', default=False)
    parser.add_option('-g', '--galaxy', dest='galaxy_run', action='store_true', default=False)
    options, args = parser.parse_args()
    
    assert(len(args) == 0)
    
    ## set Galaxy Executables
    global signature_executable, pathmark_executable
    if options.galaxy_run:
        signature_executable = os.path.join(base_directory, 'signature.py')
        pathmark_executable = os.path.join(base_directory, 'PATHMARK.py')
    
    ## run signature.py
    cmd = '%s' % (signature_executable)
    if options.data_file is not None:
        cmd += ' -d %s' % (options.data_file)
    if options.phenotype_file is not None:
        cmd += ' -p %s' % (options.phenotype_file)
    # print cmd
    os.system(cmd)
    
    ## run PATHMARK.py
    cmd = '%s' % (pathmark_executable)
    if os.path.exists('signature.tab'):
        cmd += ' -s %s' % ('signature.tab')
    if os.path.exists('null_signature.tab'):
        cmd += ' -n %s' % ('null_signature.tab')
    if os.path.exists('bootstrap_signature.tab'):
        cmd += ' -b %s' % ('bootstrap_signature.tab')
    if options.output_file is not None:
        cmd += ' -o %s' % (options.output_file)
    cmd += ' -f "%s"' % (options.filter_parameters)
    cmd += ' -t %s' % (options.heat_diffusion)
    if options.hub_filter:
        cmd += ' -u'
    if options.galaxy_run:
        cmd += ' -g'
    # print cmd
    os.system(cmd)
    
if __name__ == '__main__':
    main()
