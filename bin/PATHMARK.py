#!/usr/bin/env python
"""PATHMARK.py: identifies subnets in a score-matrix (feature x phenotype)

Usage:
  PATHMARK.py [options] scoref [scoref ...]

Options:
  -p  str                     pathway file to superimpose scores onto (default: ./*_pathway.tab)
  -f  str[,str,...]           features to run PATHMARK on (default: all)
  -s  flt;flt[,flt;flt,...]   force mean and standard deviation statistics for score files
  -b  flt;flt                 filter lower and upper boundary parameters (default: 0;0)
  -d  str                     output to user-specified directory (default: feature/)
  -l  str                     log file
  -u                          use good hub inclusion filter
  -q                          run quietly
"""
## Written By: Sam Ng
## Last Updated: 7/23/11
import os, os.path, sys, getopt, re
import networkx

from xml.dom.minidom import Document
from mPATHMARK import rCRSData
from mPATHMARK import median_mad
from mPATHMARK import mean_std
from copy import deepcopy


verbose = True
outputCleaned = True

logFile = None

def usage(code = 0):
    print __doc__
    if code != None: sys.exit(code)

def log(msg, file = None, die = False):
    if (verbose):
        if file is None:
            sys.stderr.write(msg)
        else:
            l = open(file, "a")
            l.write(msg)
            l.close()
    if (die):
        sys.exit(1)

def syscmd(cmd):
    log("running: "+cmd)
    exitstatus = os.system(cmd)
    if exitstatus != 0:
        print "Failed with exit status %i" % exitstatus
        sys.exit(10)
    log(" ... done\n")

class FormatException(Exception):
    pass

def read_paradigm_graph(handle):
    gr = networkx.MultiDiGraph()
    for line in handle:
        tmp = line.rstrip().split("\t")
        if len(tmp) == 2:
            if tmp[1] in gr.node:
                continue              
            gr.add_node( tmp[1], type=tmp[0] )
        elif len(tmp) == 3:
            if tmp[0] not in gr.node:
                raise FormatException("Missing element declaration for : %s" % (tmp[0]))
            if tmp[1] not in gr.node:
                raise FormatException("Missing element declaration for : %s" % (tmp[1]))
            gr.add_edge(tmp[0], tmp[1], interaction=tmp[2])
        else:
            raise FormatException("Bad line: %s" % (line))
    return(gr)

def reverseGraph(gr):
    rgr = networkx.MultiDiGraph()
    for node in gr.node:
        rgr.add_node(node, type = gr.node[node]["type"])
    for source in gr.edge:
        for target in gr.edge[source]:
            for edge in gr.edge[source][target]:
                rgr.add_edge(target, source,
                             interaction = gr.edge[source][target][edge]["interaction"])
    return(rgr)

def write_sif(gr, handle, edge_type_field='interaction'):
    for (e1,e2,data) in gr.edges_iter(data=True):
        interaction = data.get(edge_type_field, "pp")
        handle.write("%s\t%s\t%s\n" % (e1, interaction, e2))

def write_cytoscape_xgmml(gr, handle, scoreMap = {}):
    doc = Document()
    graph_node = doc.createElement('graph')
    graph_node.setAttribute("xmlns", "http://www.cs.rpi.edu/XGMML")


    graph_node.setAttribute("xmlns:dc", "http://purl.org/dc/elements/1.1/")
    graph_node.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink" )
    graph_node.setAttribute("xmlns:rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#" )
    graph_node.setAttribute("xmlns:cy", "http://www.cytoscape.org" )
    graph_node.setAttribute("directed", "1")
    doc.appendChild(graph_node)

    name_map = {}
    for i, n in enumerate(gr.node):
        name_map[n] = i
        node = doc.createElement('node')
        node.setAttribute('label', str(n))
        node.setAttribute('id', str(i))

        att_node = doc.createElement('att')
        att_node.setAttribute('name', "TYPE")
        att_node.setAttribute('value', str(gr.node[n]["type"]))
        att_node.setAttribute('type', "string")
        node.appendChild(att_node)
        
        att_node = doc.createElement('att')
        att_node.setAttribute('name', "LABEL")
        if gr.node[n]["type"] == "protein":
            att_node.setAttribute('value', str(n))
        else:
            att_node.setAttribute('value', "")
        att_node.setAttribute('type', "string")
        node.appendChild(att_node)
        
        att_node = doc.createElement('att')
        att_node.setAttribute('name', "SCORE")
        if str(n) in scoreMap:
            att_node.setAttribute('value', str(scoreMap[str(n)]))
        else:
            att_node.setAttribute('value', "NA")
        att_node.setAttribute('type', "real")
        node.appendChild(att_node)
        graph_node.appendChild(node)

    for source in gr.edge:
        for target in gr.edge[source]:
            for edge in gr.edge[source][target]:
                edge_node = doc.createElement("edge")
                edge_node.setAttribute("label", "%s - %s" % (source, target))
                edge_node.setAttribute("source", str(name_map[source]))
                edge_node.setAttribute("target", str(name_map[target]))
                for key, value in gr.edge[source][target][edge].items():
                    att_node = doc.createElement('att')
                    att_node.setAttribute('name', key)
                    att_node.setAttribute('value', str(value))
                    if type(value) == float:
                        att_node.setAttribute('type', "real")
                    elif type(value) == int:
                        att_node.setAttribute('type', "integer")
                    else:
                        att_node.setAttribute('type', "string")
                    edge_node.appendChild(att_node)
                graph_node.appendChild(edge_node)

    doc.writexml(handle, addindent=" ", newl="\n")

def selectLink(feature, source, target, unsignedData, filterStats, filterBounds,
               selectionRule = "OR"):
    linkScore = []
    srcScore = []
    trgScore = []
    for i in range(len(unsignedData.keys())):
        linkScore.append([unsignedData[i][feature][source],
                          unsignedData[i][feature][target]])
    for i in range(len(unsignedData.keys())):
        if linkScore[i][0] > filterStats[i][0]+filterBounds[1]*filterStats[i][1]:
            srcScore.append(2)
        elif linkScore[i][0] > filterStats[i][0]+filterBounds[0]*filterStats[i][1]:
            srcScore.append(1)
        else:
            srcScore.append(0)
        if linkScore[i][1] > filterStats[i][0]+filterBounds[1]*filterStats[i][1]:
            trgScore.append(2)
        elif linkScore[i][1] > filterStats[i][0]+filterBounds[0]*filterStats[i][1]:
            trgScore.append(1)
        else:
            trgScore.append(0)
    
    ## selection rule
    if selectionRule == "OR":
        if max(srcScore)+max(trgScore) >= 3:
            return(True)
    elif selectionRule == "AND":
        votes = 0
        for i in range(len(unsignedData.keys())):
            if srcScore[i]+trgScore[i] >= 3:
                votes += 0
        if votes == len(unsignedData.keys()):
            return(True)
    elif selectionRule == "MAIN":
        if srcScore[0]+trgScore[0] >= 3:
            return(True)
    return(False)

def getGroupSize(node, rgr):
    if node in rgr.edge:
        groupList = []
        for source in rgr.edge[node]:
            edgePass = False
            for edge in rgr.edge[node][source]:
                if rgr.edge[node][source][edge]["interaction"] in ["component>", "member>"]:
                    edgePass = True
            if edgePass:
                groupList.append(source)
        return(len(groupList))
    else:
        return(0)

def filterGeneGroups(subnetGraph, globalGraph, cutoff = 0.5):
    """clean filter"""
    outputGraph = deepcopy(subnetGraph)
    routputGraph = reverseGraph(outputGraph)
    rglobalGraph = reverseGraph(globalGraph)
    groupNodes = []
    for node in subnetGraph.node:
        if subnetGraph.node[node]["type"] in ["complex", "family"]:
            groupNodes.append(node)
    deleteNodes = ["-"]
    while (len(deleteNodes) > 0):
        deleteNodes = []
        for node in groupNodes:
            top = getGroupSize(node, routputGraph)
            bottom = getGroupSize(node, rglobalGraph)
            if bottom == 0:
                pass
            elif float(top)/float(bottom) < cutoff:
                deleteNodes.append(node)
        groupNodes = list(set(groupNodes) - set(deleteNodes))
        for node in deleteNodes:
            outputGraph.remove_node(node)
            routputGraph.remove_node(node)
    return(outputGraph)

def filterGoodHub(subnetGraph, globalGraph, nodeScore, nodeThreshold,
                  childrenTypes = ["protein"], childrenThreshold = 3,
                  includeThreshold = 0.5):
    ## initialize output variables
    oNodes = deepcopy(sNodes)
    oInteractions = deepcopy(sInteractions)
    
    ## build reverse
    rgInteractions = reverseInteractions(gInteractions)
    
    ## identify all hubs
    badHubs  = set()
    goodHubs = set()
    childrenMap = {}
    for node in gNodes:
        childrenMap[node] = set()
        if node in gInteractions:
            for child in gInteractions[node]:
                if gNodes[child] in childrenTypes:
                    childrenMap[node].update([child])
        if node not in oNodes:
            if len(childrenMap[node]) >= childrenThreshold:
                badHubs.update([node])
    
    ## iteratively add good hubs
    hubScore = {}
    currentCount = 1
    currentNodes = []
    while (currentCount > 0):
        currentCount = 0
        currentNodes = []
        currentHubs = deepcopy(badHubs)
        ## identify hubs to add
        for hub in currentHubs:
            goodChildren = set()
            for child in childrenMap[hub]:
                if child in nodeScore:
                    if nodeScore[child] >= nodeThreshold:
                        goodChildren.update([child])
                if child in goodHubs:
                    goodChildren.update([child])
            if float(len(goodChildren))/float(len(childrenMap[hub])) >= includeThreshold:
                log("> %s\n" % (hub))
                goodHubs.update([hub])
                badHubs.remove(hub)
                currentCount += 1
                currentNodes.append(hub)
                for child in childrenMap[hub]:
                    if child not in oNodes:
                        currentNodes.append(child)
        ## connect new good hubs
        while len(currentNodes) > 0:
            node = currentNodes.pop(0)
            oNodes[node] = gNodes[node]
            if node in gInteractions:
                for target in gInteractions[node]:
                    if target in oNodes:
                        if node not in oInteractions:
                            oInteractions[node] = {}
                        oInteractions[node][target] = gInteractions[node][target]
            if node in rgInteractions:
                for source in rgInteractions[node]:
                    if source in oNodes:
                        if source not in oInteractions:
                            oInteractions[source] = {}
                        oInteractions[source][node] = rgInteractions[node][source]
    return(oNodes, oInteractions)

def PATHMARK(files, globalPathway, pathmarkFeatures = None, forcedStats = None,
             filterBounds = [0, 0], outputDirectory = None, selectionRule = "OR",
             topDisconnected = 100, applyGoodHubFilter = False):
    filterString = "%s_%s" % (filterBounds[0], filterBounds[1])
    
    ## read global pathway
    f = open(globalPathway, "r")
    globalGraph = read_paradigm_graph(f)
    f.close()
    
    ## read scores
    signedData = {}
    unsignedData = {}
    for i in range(len(files)):
        signedData[i] = rCRSData(files[i])
        unsignedData[i] = {}
        for j in signedData[i].keys():
            unsignedData[i][j] = {}
            for k in signedData[i][j].keys():
                try:
                    unsignedData[i][j][k] = abs(float(signedData[i][j][k]))
                except ValueError:
                    unsignedData[i][j][k] = "NA"
    
    ## iterate features
    for feature in unsignedData[0].keys():
        if pathmarkFeatures is not None:
            if feature not in pathmarkFeatures:
                continue
        ## initiate pathmark graph
        pathmarkGraph = networkx.MultiDiGraph()
        
        ## compute score statistics
        filterStats = []
        if forcedStats is None:
            for i in range(len(unsignedData.keys())):
                filterStats.append(median_mad(unsignedData[i][feature].values()))
        else:
            for i in forcedStats.split(","):
                (v1, v2) = i.split(";")
                filterStats.append((float(v1), float(v2)))
        
        ## log statistics
        log("%s\t%s;%s" % (feature, filterStats[0][0], filterStats[0][1]), file = logFile)
        for i in range(1, len(filterStats)):
            log(",%s;%s" % (filterStats[i][0], filterStats[i][1]), file = logFile)
        log("\n", file = logFile)
        
        ## add top scoring links
        for source in globalGraph.edge:
            if source not in unsignedData[0][feature]:
                continue
            elif unsignedData[0][feature][source] == "NA":
                continue
            for target in globalGraph.edge[source]:
                if target not in unsignedData[0][feature]:
                    continue
                elif unsignedData[0][feature][target] == "NA":
                    continue
                
                if selectLink(feature, source, target, unsignedData, filterStats,
                              filterBounds, selectionRule = selectionRule):
                    if source not in pathmarkGraph.node:
                        pathmarkGraph.add_node(source,
                                               type = globalGraph.node[source]["type"])
                    if target not in pathmarkGraph.node:
                        pathmarkGraph.add_node(target,
                                               type = globalGraph.node[target]["type"])
                    for edge in globalGraph.edge[source][target]:
                        pathmarkGraph.add_edge(source, target, interaction =
                                    globalGraph.edge[source][target][edge]["interaction"])
        
        ## add top scoring disconnected nodes
        sortedTop = []
        for node in unsignedData[0][feature]:
            if node not in globalGraph.node:
                continue
            if globalGraph.node[node] in ["protein"]:
                sortedTop.append(node)
        if len(sortedTop) > 0:
            sortedTop.sort(lambda x, y: cmp(unsignedData[0][feature][y],unsignedData[0][feature][x]))
            while (unsignedData[0][feature][sortedTop[0]] == "NA"):
                sortedTop.pop(0)
                if len(sortedTop) == 0:
                    break
            for i in range(topDisconnected):
                if i > len(sortedTop)-1:
                    break
                if unsignedData[0][feature][sortedTop[i]] < filterStats[0][0]+filterBounds[0]*filterStats[0][1]:
                    break
                if sortedTop[i] not in pathmarkGraph.node:
                    if "__DISCONNECTED__" not in pathmarkGraph.node:
                        pathmarkGraph.add_node("__DISCONNECTED__", type = "abstract")
                    pathmarkGraph.add_node(sortedTop[i],
                                       type = globalGraph.node[sortedTop[i]]["type"])
                    pathmarkGraph.add_edge(sortedTop[i], "__DISCONNECTED__",
                                           interaction = "-disconnected-")
        
        ## apply filters
        if applyGoodHubFilter:
            pass
        #    (pNodes, pInteractions) = filterGoodHub(pNodes, pInteractions, gNodes, gInteractions,
        #                                            unsignedData[0][feature], filterStats[0][0],
        #                                            childrenTypes = ["protein"],
        #                                            childrenThreshold = 3,
        #                                            includeThreshold = 0.5)
        pathmarkGraph = filterGeneGroups(pathmarkGraph, globalGraph)
        
        ## output networks
        if outputDirectory == None:
            currentDirectory = feature
        else:
            currentDirectory = outputDirectory
        if not os.path.exists(currentDirectory):
            syscmd("mkdir %s" % (currentDirectory))
        f = open("%s/%s_%s_nodrug.xgmml" % (currentDirectory, feature, filterString), "w")
        write_cytoscape_xgmml(pathmarkGraph, f, signedData[0][feature])
        f.close()
        
        f = open("%s/%s_%s_nodrug.xgmml" % (currentDirectory, feature, filterString), "w")
        write_sif(pathmarkGraph, f)        
        f.close()

if __name__ == "__main__":
    ## parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], "p:f:s:b:d:l:uq")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) < 1:
        print "incorrect number of arguments"
        usage(1)
    
    globalPathway = None
    pathmarkFeatures = None
    forcedStats = None
    filterBounds = [0, 0]
    outputDirectory = None
    selectionRule = "OR"
    applyGoodHubFilter = False
    for o, a in opts:
        if o == "-p":
            globalPathway = a
        elif o == "-f":
            pathmarkFeatures = a.split(";")
        elif o == "-s":
            forcedStats = a
            if os.path.exists(forcedStats):
                f = open(forcedStats, "r")
                forcedStats = f.readline().rstrip().split("\t")[1]
                f.close()
        elif o == "-b":
            (v1, v2) = a.split(";")
            filterBounds = [float(v1), float(v2)]
            filterBounds.sort()
        elif o == "-d":
            outputDirectory = a.rstrip("/")
        elif o == "-l":
            logFile = a
        elif o == "-u":
            applyGoodHubFilter = True
        elif o == "-q":
            verbose = False
    if globalPathway is None:
        for file in os.listdir("."):
            if file.endswith("_pathway.tab"):
                globalPathway = file
    assert globalPathway is not None
    
    ## run
    PATHMARK(args, globalPathway, pathmarkFeatures = pathmarkFeatures,
             forcedStats = forcedStats,  filterBounds = filterBounds,
             outputDirectory = outputDirectory, selectionRule = selectionRule,
             topDisconnected = 100, applyGoodHubFilter = applyGoodHubFilter)
