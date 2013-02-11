#!/usr/bin/env python
"""layout-cytoscape.py: Generate cytoscape sessions from .sif and
   .NA files

Usage:
  layout-cytoscape.py [options] feature

Options:
  -c path   path to cytoscape
  -v path   path to vizmap
  -p path   path to plugins
  -l str    cytoscape layout
  -q        run quietly
  
Notes:
  Remember to set up your path (-c, -v, -p)
"""
## Written By: Sam Ng
## Last Updated: 5/17/11
import os, os.path, sys, getopt, re

verbose = True

if os.path.exists("/projects/sysbio/apps/java/Cytoscape/cytoscape-2.8.0"):
    cysPath = "/projects/sysbio/apps/java/Cytoscape/cytoscape-2.8.0/cytoscape.sh"
    vizPath = "/projects/sysbio/apps/java/Cytoscape/cytoscape-2.8.0/vizmap/vizmap_nodrug.props"
    pluginPath = "/projects/sysbio/apps/java/Cytoscape/cytoscape-2.8.0/plugins"
else:
    cysPath = "/Applications/Cytoscape_v2.8.2/cytoscape.sh"
    vizPath = "/Users/OrosTheAvenger/Desktop/Dropbox/My_Research/bin/pathmark-scripts/vizmaps/vizmap_nodrug.props"
    pluginPath = "/Applications/Cytoscape_v2.8.2/plugins"

layoutSpec = 'layout.default="force-directed" defaultVisualStyle="Local-Red-Blue-On-White"'
netExtension = ".sif"

def usage(code = 0):
    print __doc__
    if code != None: sys.exit(code)

def log(msg):
    if (verbose):
        sys.stderr.write(msg)

def main(args):
    ## Parse arguments
    try:
        opts, args = getopt.getopt(args, "c:v:p:l:q")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) != 1:
        print "incorrect number of arguments"
        usage(1)
    
    feature = args[0]
    
    global verbose, cysPath, vizPath, pluginPath, layoutSpec
    for o, a in opts:
        if o == "-c":
            cysPath = a
        elif o == "-v":
            vizPath = a
        elif o == "-p":
            pluginPath = a
        elif o == "-l":
            layoutSpec = a
        elif o == "-q":
            verbose = False
    
    ## Check structure
    assert os.path.exists("%s/LABEL.NA" % (feature))
    assert os.path.exists("%s/TYPE.NA" % (feature))
    assert os.path.exists("%s/SCORE.NA" % (feature))
    assert os.path.exists("%s" % (feature))
    
    naFiles = ["%s/LABEL.NA" % (feature), "%s/TYPE.NA" % (feature), "%s/SCORE.NA" % (feature)]
    if os.path.exists("%s/BORDER.NA" % (feature)):
        naFiles.append("%s/BORDER.NA" % (feature))
    
    ## Identify nets with feature
    sifFiles = list()
    for i in os.listdir("%s" % (feature)):
        if i.endswith(netExtension):
            sifFiles.append(feature+"/"+i)
                
    ## Launch cytoscape
    cmd = "%s -N %s -n %s -V %s -p %s -P %s" % (cysPath, " ".join(sifFiles), " ".join(naFiles), vizPath, pluginPath, layoutSpec)
    log(cmd+"\n")
    os.system(cmd)

if __name__ == "__main__":
    main(sys.argv[1:])
