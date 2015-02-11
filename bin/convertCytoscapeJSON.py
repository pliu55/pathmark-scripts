#!/usr/bin/env python
"""
convertCytoscapeJSON.py
    by Sam Ng
"""
import logging, math, os, re, shutil, sys, types
from copy import deepcopy

import xmltodict
import simplejson as json

from optparse import OptionParser

def convertCytoscapeXGMMLtoCytoscapeJSON(output_file, input_file):
    f = open(input_file, 'r')
    xgmml_object = xmltodict.parse(f.read())
    f.close()
    json_object = {}
    json_object["nodes"] = []
    for node_object in xgmml_object["graph"]["node"]:
        index = int(node_object['@id'])
        node_element = {"data" : {}}
        node_element["data"]["id"] = str(index)
        for attribute_object in node_object["att"]:
            if attribute_object["@type"] == "real":
                node_element["data"][str(attribute_object["@name"])] = float(attribute_object["@value"])
            else:
                node_element["data"][str(attribute_object["@name"])] = str(attribute_object["@value"])
        json_object["nodes"].append(node_element)
    json_object["edges"] = []
    for edge_object in xgmml_object["graph"]["edge"]:
        edge_element = {"data" : {}}
        edge_element["data"]["source"] = str(edge_object['@source'])
        edge_element["data"]["target"] = str(edge_object['@target'])
        for attribute_object in edge_object["att"]:
            if attribute_object["@type"] == "real":
                edge_element["data"][str(attribute_object["@name"])] = float(attribute_object["@value"])
            else:
                edge_element["data"][str(attribute_object["@name"])] = str(attribute_object["@value"])
        json_object["edges"].append(edge_element)
    o = open(output_file, "w")
    o.write(json.dumps(json_object, indent = "  "))
    o.close()

def main():
    ## parse arguments
    parser = OptionParser(usage = "%prog [options] input_file output_file")
    options, args = parser.parse_args()
    logging.info("options: %s" % (str(options)))
    
    if len(args) != 2:
        logging.error("ERROR: incorrect number of arguments\n")
        sys.exit(1)
    input_file = os.path.abspath(args[0])
    output_file = args[1]
    
    convertCytoscapeXGMMLtoCytoscapeJSON(output_file, input_file)

if __name__ == "__main__":
    main()