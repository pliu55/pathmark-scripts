#!/usr/bin/env python
"""
PATHMARK.py
    by Sam Ng
"""
import logging, math, os, re, sys, types
from copy import deepcopy

import pandas
import networkx
import matplotlib.pyplot as plt
from xml.dom.minidom import Document

from optparse import OptionParser
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

#### NOTE BLOCK
#### - Add in option to collapse the input pathway

## logger
logging.basicConfig(filename = "pathmark.log", level = logging.INFO)

## pm classes
class Parameters:
    """
    Stores parameters used for this [PATHMARK.py specific]
    """
    def __init__(self, pathway_file, filter_parameters, diffusion_time, hub_filter):
        self.pathway_file = pathway_file
        self.filter_parameters = [float(i) for i in filter_parameters.split(";")]
        self.filter_parameters.sort()
        self.top_disconnected = 100
        self.heat_diffusion_time = float(diffusion_time)
        self.hub_filter = hub_filter
        self.parameters_string = "%s_%s" % (self.filter_parameters[0], self.filter_parameters[1])

class Pathway:
    """
    Paradigm compatible pathway class [2014-6-9]
    Dependencies: logger
    """
    def __init__(self, input, pid = None):
        self.nodes = {}
        self.interactions = {}
        self.pid = pid
        if type(input) == types.TupleType:
            self.nodes = deepcopy(input[0])
            self.interactions = deepcopy(input[1])
        elif type(input) == types.StringType:
            (self.nodes, self.interactions) = self.readSPF(input)
        else:
            logger("ERROR: invalid input for pathway import (%s)\n" % (input), die = True)
    def readSPF(self, input_file):
        nodes = {}
        interactions = {}
        f = open(input_file, "r")
        for line in f:
            if line.isspace():
                continue
            pline = line.rstrip().split("\t")
            if len(pline) == 2:
                nodes[pline[1]] = pline[0]
            elif len(pline) == 3:
                if pline[0] not in interactions:
                    interactions[pline[0]] = {}
                if pline[1] not in interactions[pline[0]]:
                    interactions[pline[0]][pline[1]] = pline[2]
                else:
                    interactions[pline[0]][pline[1]] += ";" + pline[2]
            else:
                logger("ERROR: line length not 2 or 3 (%s)\n" % (line), die = True)
        f.close()
        return(nodes, interactions)
    def writeSPF(self, output_file, reverse = False):
        o = open(output_file, "w")
        for node in self.nodes:
            o.write("%s\t%s\n" % (self.nodes[node], node))
        for source in self.interactions:
            for target in self.interactions[source]:
                for interaction in self.interactions[source][target].split(";"):
                    o.write("%s\t%s\t%s\n" % (source, target, interaction))
        o.close()
    def writeSIF(self, output_file):
        o = open(output_file, "w")
        for source in self.interactions:
            for target in self.interactions[source]:
                for interaction in self.interactions[source][target].split(";"):
                    o.write("%s\t%s\t%s\n" % (source, interaction, target))
        o.close()
    def networkx(self):
        networkx_graph = networkx.MultiDiGraph()
        for node in self.nodes:
            networkx_graph.add_node(node, type = self.nodes[node])
        for source in self.interactions:
            for target in self.interactions[source]:
                for interaction in self.interactions[source][target].split(";"):
                    networkx_graph.add_edge(source, target, interaction = interaction)
        return(networkx_graph)
    def reverse(self):
        reversed_interactions = {}
        for source in self.interactions:
            for target in self.interactions[source]:
                if target not in reversed_interactions:
                    reversed_interactions[target] = {}
                reversed_interactions[target][source] = self.interactions[source][target]
        return(reversed_interactions)
    def getShortestPaths(self, source, target, max_distance = None):
        shortest_paths = []
        all_walks = [[source]]
        while len(shortest_paths) == 0:
            next_walks = []
            while len(all_walks) > 0:
                current_walk = all_walks.pop()
                if current_walk[-1] not in self.interactions:
                    continue
                for intermediate in self.interactions[current_walk[-1]]:
                    if intermediate in current_walk:
                        continue
                    if intermediate == target:
                        complete_walk = current_walk + [intermediate]
                        shortest_paths.append([(complete_walk[index],
                                              self.interactions[complete_walk[index]][complete_walk[index + 1]],
                                              complete_walk[index + 1])
                                              for index in range(len(complete_walk) - 1)])
                    next_walks.append(current_walk + [intermediate])
            if len(next_walks) == 0:
                break
            if max_distance is not None:
                if len(next_walks[0]) >= max_distance + 1:
                    break
            all_walks = deepcopy(next_walks)
        return(shortest_paths)
    def getAllPaths(self, source, target, max_distance):
        all_paths = []
        all_walks = [[source]]
        for distance in range(1, max_distance + 1):
            next_walks = []
            while len(all_walks) > 0:
                current_walk = all_walks.pop()
                if current_walk[-1] not in self.interactions:
                    continue
                for intermediate in self.interactions[current_walk[-1]]:
                    if intermediate in current_walk:
                        continue
                    if intermediate == target:
                        complete_walk = current_walk + [intermediate]
                        all_paths.append([(complete_walk[index], 
                                         self.interactions[complete_walk[index]][complete_walk[index + 1]],
                                         complete_walk[index + 1])
                                         for index in range(len(complete_walk) - 1)])
                    next_walks.append(current_walk + [intermediate])
            if len(next_walks) == 0:
                break
            all_walks = deepcopy(next_walks)
        return(all_paths)
    def getNeighbors(self, node, max_distance):
        reversed_interactions = self.reverse()
        seen_nodes = set([node])
        border_nodes = [node]
        frontier_nodes = []
        for distance in range(max_distance):
            while len(border_nodes) > 0:
                current_node = border_nodes.pop()
                if current_node in self.interactions:
                    for target in self.interactions[current_node]:
                        if target not in seen_nodes:
                            seen_nodes.update([target])
                            frontier_nodes.append(target)
                if current_node in reversed_interactions:
                    for source in reversed_interactions[current_node]:
                        if source not in seen_nodes:
                            seen_nodes.update([source])
                            frontier_nodes.append(source)
            border_nodes = deepcopy(frontier_nodes)
            frontier_nodes = []
        return(list(seen_nodes))
    def subsetPathway(self, subsetted_nodes):
        subsetted_pathway = Pathway( ({}, {}) )
        for source in subsetted_nodes:
            if source not in subsetted_pathway.nodes:
                subsetted_pathway.nodes[source] = self.nodes[source]
            if source in self.interactions:
                for target in self.interactions[source]:
                    if target not in subsetted_nodes:
                        continue
                    if target not in subsetted_pathway.nodes:
                        subsetted_pathway.nodes[target] = self.nodes[target]
                    if source not in subsetted_pathway.interactions:
                        subsetted_pathway.interactions[source] = {}
                    subsetted_pathway.interactions[source][target] = self.interactions[source][target]
        return(subsetted_pathway)
    def appendPathway(self, append_pathway):
        for source in append_pathway.interactions:
            if source not in self.nodes:
                self.nodes[source] = append_pathway.nodes[source]
            for target in append_pathway.interactions[source]:
                if target not in self.nodes:
                    self.nodes[target] = append_pathway.nodes[target]
                if source not in self.interactions:
                    self.interactions[source] = {}
                if target not in self.interactions[source]:
                    self.interactions[source][target] = append_pathway.interactions[source][target]
                else:
                    self.interactions[source][target] = ";".join(list(set(self.interactions[source][target].split(";")) |
                                                                      set(append_pathway.interactions[source][target].split(";"))))

## pm functions
def logger(message, file = None, die = False):
    """
    Writes messages to standard error [2014-3-1]
    """
    if file is None:
        sys.stderr.write(message)
    else:
        o = open(file, "a")
        o.write(message)
        o.close()
    if die:
        sys.exit(1)

def diffuseHeat(signature_frame, global_graph, parameters):
    """
    Refactoring of https://github.com/ucscCancer/pathway_tools/blob/master/scripts/diffuse.py [613be7d5baad339e8ddc852be6e10baff0cf8f9c]
    """
    from array import array
    from numpy import genfromtxt, dot
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import expm
    
    class SciPYKernel:
        def __init__(self):
            self.laplacian = None
            self.index2node = None
            self.kernel = None
            self.labels = None
        def readKernel(self, kernel_file):
            self.kernel = coo_matrix(genfromtxt(kernel_file, delimiter = "\t")[1:, 1:])
            f = open(kernel_file, "r")
            self.labels = f.readline().rstrip().split("\t")[1:]
            f.close()
        def makeKernel(self, networkx_graph, diffusion_time = 0.1):
            ## parse the network, build the graph laplacian
            edges, nodes, node_out_degrees = self.parseGraph(networkx_graph)
            num_nodes = len(nodes)
            node_order = list(nodes)
            index2node = {}
            node2index = {}
            for i in range(0, num_nodes):
                index2node[i] = node_order[i]   
                node2index[node_order[i]] = i
            ## construct the diagonals
            row = array("i")
            col = array("i")
            data = array("f")
            for i in range(0, num_nodes):
                ## diag entries: out degree
                degree = 0 
                if index2node[i] in node_out_degrees:
                    degree = node_out_degrees[index2node[i]]
                ## append to the end
                data.insert(len(data), degree)  
                row.insert(len(row), i) 
                col.insert(len(col), i) 
            ## add edges
            for i in range(0, num_nodes):
                for j in range(0, num_nodes):
                    if i == j:
                        continue
                    if (index2node[i], index2node[j]) not in edges:
                        continue
                    ## append index to i-th row, j-th column
                    row.insert(len(row), i)
                    col.insert(len(col), j)
                    ## -1 for laplacian edges
                    data.insert(len(data), -1)
            ## graph laplacian
            L = coo_matrix((data,(row, col)), shape=(num_nodes,num_nodes)).tocsc()
            self.laplacian = L
            self.index2node = index2node
            self.kernel = expm(-1*diffusion_time*L)
            self.labels = node_order
        def printLaplacian(self):
            cx = self.laplacian.tocoo()
            for i,j,v in zip(cx.row, cx.col, cx.data):
                a = self.index2node[i]
                b = self.index2node[j]
                print "\t".join([a,b,str(v)])
        def parseGraph(self, networkx_graph):
            edges = set()
            nodes = set()
            degrees = {}
            for source in networkx_graph.edge:
                for target in networkx_graph.edge[source]:
                    if (source, target) in edges:
                        continue
                    edges.add((source, target))
                    edges.add((target, source))
                    nodes.add(source)
                    nodes.add(target)
                    if source not in degrees:
                        degrees[source] = 0
                    if target not in degrees:
                        degrees[target] = 0
                    degrees[source] += 1
                    degrees[target] += 1
            return (edges, nodes, degrees)
        def writeKernel(self, output_file):
            o = open(output_file, "w")
            cx = self.kernel.tocoo()
            edges = {}
            for i, j, v in zip(cx.row, cx.col, cx.data):
                a = self.index2node[i]
                b = self.index2node[j]
                edges[(a, b)] = str(v)
            ## iterate through rows
            ## sort labels in alphabetical order
            o.write("Key\t" + "\t".join(sorted(self.labels)) + "\n")
            for nodeA in sorted(self.labels):
                printstr = nodeA
                # through columns
                for nodeB in sorted(self.labels):
                    if (nodeA, nodeB) in edges:
                        printstr += "\t" + edges[(nodeA, nodeB)]	
                    else:
                        printstr += "\t0"
                o.write(printstr + "\n")
            o.close()
        def kernelMultiplyOne(self, vector):
            array = []
            ## loop over gene names in the network kernel: add the starting value if 
            ## it's present in the supplied input vector
            for label in self.labels:
                if label in vector:
                    array.append(vector[label])
                else:
                    array.append(0)
            ## take the dot product
            value = self.kernel*array
            return_vec = {}
            idx = 0
            for label in self.labels:
                return_vec[label] = float(value[idx])
                idx += 1
            return return_vec
        @staticmethod
        def getAngle(v1, v2):
            arry1 = []
            arry2 = []
            for key in v1:
                arry1.append(float(v1[key]))
                arry2.append(float(v2[key]))
            mag_1 = math.sqrt(dot(arry1,arry1))
            mag_2 = math.sqrt(dot(arry2,arry2))
            cos_theta = dot(arry1,arry2)/(mag_1*mag_2)
            return math.acos(cos_theta)
        def diffuse(self, vector, reverse=False):
            ## reverse is not used: heat diffusion is undirected
            diffused_vector = self.kernelMultiplyOne(vector)
            return diffused_vector
    
    def collapseNetworkX(networkx_graph):
        collapsed_map = {}
        collapsed_graph = networkx.MultiDiGraph()
        reversed_networkx_graph = reverseNetworkX(networkx_graph)
        for target in reversed_networkx_graph.edge:
            for source in reversed_networkx_graph.edge[target]:
                if pathmark_graph.node[node]["type"] in ["complex", "family"]:
                    ####
                    pass
                for edge in reversed_networkx_graph.edge[target][source]:
                    ####
                    pass
        return(collapsed_map, collapsed_graph)
    
    def expandNetworkX(diffused_heats, collapsed_map):
        ## watch out for complex loops
        return(0)
     
    ## https://github.com/ucscCancer/pathway_tools/blob/master/scripts/diffuse.py
    kernel_file = "analysis/kernels/%s.kernel" % (".".join(parameters.pathway_file.split("/")[-1].split(".")[:-1]))
    diffuser = SciPYKernel()
    if os.path.exists(kernel_file):
        diffuser.readKernel(kernel_file)
    else:
        diffuser.makeKernel(global_graph, diffusion_time = parameters.heat_diffusion_time)
        if not os.path.exists("analysis/kernels"):
            os.mkdir("analysis/kernels")
        diffuser.writeKernel(kernel_file)
    input_heats = abs(signature_frame.icol(0))
    diffused_heats = pandas.DataFrame(pandas.Series(diffuser.diffuse(input_heats, reverse = False), name = signature_frame.columns[0]))
    return(diffused_heats)

def convertListToFloat(input_list):
    """
    Converts a list of strings to floats and removes non-float values [2014-6-7]
    """
    float_list = []
    for value in input_list:
        try:
            float_value = float(value)
            if float_value != float_value:
                raise ValueError
            float_list.append(float_value)
        except ValueError:
            continue
    return(float_list)

def computeMean(input_list, null = "NA", return_sd = False, sample_sd = True):
    """
    Computes mean and optionally the sd [2014-6-7]
    Dependencies: convertListToFloat
    """
    float_list = convertListToFloat(input_list)
    if len(float_list) == 0:
        mean = null
        sd = null
    else:
        mean = sum(float_list)/float(len(float_list))
        sd = 0.0
        for value in float_list:
            sd += (value - mean)**2
        if len(float_list) > 1:
            if sample_sd:
                sd = math.sqrt(sd/(len(float_list) - 1))
            else:
                sd = math.sqrt(sd/len(float_list))
        else:
            sd = 0.0
    if return_sd:
        return(mean, sd)
    else:
        return(mean)

def computeMedian(input_list, null = "NA", return_mad = False):
    """
    Computes median and optionally the mad [2014-6-7]
    Dependencies: convertListToFloat
    """
    float_list = convertListToFloat(input_list)
    float_list.sort()
    if len(float_list) == 0:
        median = null
        mad = null
    else:
        float_list_size = len(float_list)
        if float_list_size%2:
            median = float_list[float_list_size/2]
        else:
            median = (float_list[(float_list_size/2) - 1] + float_list[float_list_size/2])/2.0
        if not return_mad:
            return(median)
        else:
            ad_list = []
            for value in float_list:
                ad_list.append(abs(value - median))
            mad = computeMedian(ad_list)
            return(median, mad)

def computeStatistics(score_frame, method = "median"):
    """
    Passthrough function for computing statistics [PATHMARK.py specific]
    Dependencies: computeMean, computeMedian, convertListToFloat
    """
    values = [value[0] for value in score_frame.values.tolist()]
    if method == "mean":
        score_statistics = computeMean(values, return_sd = True)
    elif method == "median":
        score_statistics = computeMedian(values, return_mad = True)
    return(score_statistics)

def getPathmark(score_frame, global_graph, parameters, forced_statistics = None):
    """
    Selects the PATHMARK and returns it as a networkx object [PATHMARK.py specific]
    Dependencies: computeStatistics, computeMean, computeMedian, convertListToFloat
    """
    ## compute score_statistics
    if forced_statistics != None:
        assert(len(forced_statistics) == 2)
        score_statistics = forced_statistics
    else:
        score_statistics = computeStatistics(score_frame)
    score_thresholds = (score_statistics[0] + parameters.filter_parameters[0]*score_statistics[1],
                        score_statistics[0] + parameters.filter_parameters[1]*score_statistics[1])
    
    ## initiate pathmark_graph
    pathmark_graph = networkx.MultiDiGraph()
    
    ## store edges that pass selection parameters
    for source in global_graph.edge:
        if source not in score_frame.index:
            continue
        source_value = float(score_frame.icol(0)[source])
        if source_value != source_value:
            continue
        for target in global_graph.edge[source]:
            if target == source:
                continue
            if target not in score_frame.index:
                continue
            target_value = float(score_frame.icol(0)[target])
            if target_value != target_value:
                continue
            edge_values = [source_value, target_value]
            edge_values.sort()
            if abs(edge_values[0]) >= score_thresholds[0] and abs(edge_values[1]) >= score_thresholds[1]:
                if source not in pathmark_graph.node:
                    pathmark_graph.add_node(source,
                                            type = global_graph.node[source]["type"])
                if target not in pathmark_graph.node:
                    pathmark_graph.add_node(target,
                                            type = global_graph.node[target]["type"])
                for edge in global_graph.edge[source][target]:
                    pathmark_graph.add_edge(source, target, interaction =
                                            global_graph.edge[source][target][edge]["interaction"])
    
    ## store high scoring isolated nodes
    sorted_nodes = []
    for node in score_frame.index:
        if node not in global_graph.node:
            continue
        node_value = float(score_frame.icol(0)[node])
        if global_graph.node[node]["type"] in ["protein"] and node_value == node_value:
            sorted_nodes.append(node)
    if len(sorted_nodes) > 0:
        sorted_nodes.sort(lambda x, y: cmp(float(score_frame.icol(0)[y]), float(score_frame.icol(0)[x])))
        for index in range(parameters.top_disconnected):
            if index > len(sorted_nodes) - 1:
                break
            if abs(float(score_frame.icol(0)[sorted_nodes[index]])) < score_thresholds[1]:
                break
            if sorted_nodes[index] not in pathmark_graph.node:
                if "__DISCONNECTED__" not in pathmark_graph.node:
                    pathmark_graph.add_node("__DISCONNECTED__", type = "abstract")
                pathmark_graph.add_node(sorted_nodes[index],
                                        type = global_graph.node[sorted_nodes[index]]["type"])
                pathmark_graph.add_edge(sorted_nodes[index], "__DISCONNECTED__",
                                        interaction = "-disconnected-")
    
    ## doing a pass through all nodes to catch self-interactions
    for node in pathmark_graph.node:
        if node in global_graph.edge:
            if node in global_graph.edge[node]:
                for edge in global_graph.edge[node][node]:
                    pathmark_graph.add_edge(node, node, interaction =
                                            global_graph.edge[node][node][edge]["interaction"])
    ## note that edges will be missed if both nodes were between the two thresholds
    return(pathmark_graph)

def reverseNetworkX(networkx_graph):
    """
    Flips the direction of a networkx graph [PATHMARK.py specific]
    """
    reversed_networkx_graph = networkx.MultiDiGraph()
    for node in networkx_graph.node:
        reversed_networkx_graph.add_node(node, type = networkx_graph.node[node]["type"])
    for source in networkx_graph.edge:
        for target in networkx_graph.edge[source]:
            for edge in networkx_graph.edge[source][target]:
                reversed_networkx_graph.add_edge(target, source,
                             interaction = networkx_graph.edge[source][target][edge]["interaction"])
    return(reversed_networkx_graph)

def filterGeneGroups(pathmark_graph, global_graph, include_threshold = 0.5):
    """
    Applies a filter to remove complexes and families with less than the threshold number of members [PATHMARK.py specific]
    Dependencies: reverseNetworkX
    """
    def getGroupSize(node, reversed_graph):
        if node in reversed_graph.edge:
            group_list = []
            for source in reversed_graph.edge[node]:
                edge_valid = False
                for edge in reversed_graph.edge[node][source]:
                    if reversed_graph.edge[node][source][edge]["interaction"] in ["component>", "member>"]:
                        edge_valid = True
                if edge_valid:
                    group_list.append(source)
            return(len(group_list))
        else:
            return(0)
    
    filtered_graph = deepcopy(pathmark_graph)
    reversed_filtered_graph = reverseNetworkX(filtered_graph)
    reversed_global_graph = reverseNetworkX(global_graph)
    group_nodes = []
    for node in pathmark_graph.node:
        if pathmark_graph.node[node]["type"] in ["complex", "family"]:
            group_nodes.append(node)
    remove_nodes = [""]
    while len(remove_nodes) > 0:
        remove_nodes = []
        for node in group_nodes:
            top_size = getGroupSize(node, reversed_filtered_graph)
            bottom_size = getGroupSize(node, reversed_global_graph)
            if bottom_size == 0:
                pass
            elif float(top_size)/float(bottom_size) < include_threshold:
                remove_nodes.append(node)
        group_nodes = list(set(group_nodes) - set(remove_nodes))
        for node in remove_nodes:
            filtered_graph.remove_node(node)
            reversed_filtered_graph.remove_node(node)
    return(filtered_graph)

def filterGoodHub(pathmark_graph, global_graph, score_map, node_threshold, children_types = ["protein"], children_threshold = 3, include_threshold = 0.66):
    """
    Applies a filter to include hubs with greater than the threshold number of members [PATHMARK.py specific]
    Dependencies: reverseNetworkX
    """
    l = open('filter_hub.log', 'w')
    
    filtered_graph = deepcopy(pathmark_graph)
    reversed_global_graph = reverseNetworkX(global_graph)
    
    ## identify all hubs
    hub_bad = set()
    hub_good = set()
    children_map = {}
    for node in global_graph.node:
        children_map[node] = set()
        for child in global_graph.edge[node]:
            if global_graph.node[child]["type"] in children_types:
                children_map[node].update([child])
        if node not in filtered_graph.node:
            if len(children_map[node]) >= children_threshold:
                hub_bad.update([node])
    
    ## iteratively add good hubs
    hub_score = {}
    current_count = 1
    current_nodes = []
    while current_count > 0:
        current_count = 0
        current_nodes = []
        current_hubs = deepcopy(hub_bad)
        ## identify hubs to add
        for hub in current_hubs:
            children_good = set()
            for child in children_map[hub]:
                if child in score_map:
                    if score_map[child] >= node_threshold:
                        children_good.update([child])
                if child in hub_good:
                    children_good.update([child])
            if float(len(children_good))/float(len(children_map[hub])) >= include_threshold:
                l.write("> %s\t%s\t%s\n" % (hub, len(children_good), len(children_map[hub])))
                hub_good.update([hub])
                hub_bad.remove(hub)
                current_count += 1
                current_nodes.append(hub)
                for child in children_good:
                    if child in score_map:
                        l.write("  %s\t%s\n" % (child, score_map[child]))
                    else:
                        l.write("  %s\t-\n" % (child))
                    if child not in filtered_graph.node:
                        current_nodes.append(child)
        ## connect new good hubs
        while len(current_nodes) > 0:
            node = current_nodes.pop(0)
            filtered_graph.add_node(node, type = global_graph.node[node]["type"])
            for target in global_graph.edge[node]:
                if target in filtered_graph.node:
                    for edge in global_graph.edge[node][target]:
                        filtered_graph.add_edge(node, target, interaction =
                                                global_graph.edge[node][target][edge]["interaction"])
            for source in reversed_global_graph.edge[node]:
                if source in filtered_graph.node:
                    for edge in reversed_global_graph.edge[node][source]:
                        filtered_graph.add_edge(source, node, interaction =
                                                reversed_global_graph.edge[node][source][edge]["interaction"])
    l.close()
    return(filtered_graph)

def getLargestSubgraph(pathmark_graph):
    subgraphs = networkx.weakly_connected_component_subgraphs(pathmark_graph)
    if len(subgraphs) == 0:
        largest_subgraph = networkx.MultiDiGraph()
    elif "__DISCONNECTED__" in subgraphs[0].node:
        if len(subgraphs) == 1:
            largest_subgraph = networkx.MultiDiGraph()
        else:
            largest_subgraph = subgraphs[1]
    else:
        largest_subgraph = subgraphs[0]
    return(largest_subgraph)

def readNetworkXSIF(input_file, node_type_map = None, default_node_type = "abstract"):
    """
    Read in a SIF format file into a networkx object [PATHMARK.py specific]
    """
    networkx_graph = networkx.MultiDiGraph()
    f = open(input_file, "r")
    for line in f:
        if line.isspace():
            continue
        (source, interaction, target) = line.rstrip().split("\t")
        if source not in networkx_graph.node:
            if node_type_map is not None:
                if source in node_type_map:
                    networkx_graph.add_node(source, type = node_type_map[source])
                else:
                    networkx_graph.add_node(source, type = default_node_type)
            else:
                networkx_graph.add_node(source, type = default_node_type)
        if target not in networkx_graph.node:
            if node_type_map is not None:
                if target in node_type_map:
                    networkx_graph.add_node(target, type = node_type_map[target])
                else:
                    networkx_graph.add_node(target, type = default_node_type)
            else:
                networkx_graph.add_node(target, type = default_node_type)
        networkx_graph.add_edge(source, target, interaction = interaction)
    f.close()
    return(networkx_graph)

def writeNetworkXSIF(output_file, networkx_graph, edge_type_field = "interaction"):
    """
    Write a SIF format file from a networkx object [PATHMARK.py specific]
    """
    o = open(output_file, "w")
    for (source, target, data) in networkx_graph.edges_iter(data = True):
        interaction = data.get(edge_type_field, "pp")
        o.write("%s\t%s\t%s\n" % (source, interaction, target))
    o.close()

def writeNetworkXSPF(output_file, networkx_graph, edge_type_field = "interaction"):
    """
    Write a SPF format file from a networkx object [PATHMARK.py specific]
    """
    o = open(output_file, "w")
    for node in networkx_graph.node:
        o.write("%s\t%s\n" % (str(networkx_graph.node[node]["type"]), node))
    for (source, target, data) in networkx_graph.edges_iter(data = True):
        interaction = data.get(edge_type_field, "pp")
        o.write("%s\t%s\t%s\n" % (source, interaction, target))
    o.close()

def writeNetworkXCytoscapeXGMML(output_file, networkx_graph, signature_map, score_map, bootstrap_map):
    """
    Write a Cytoscape XGMML format file from a networkx object [PATHMARK.py specific]
    """
    doc = Document()
    graph_node = doc.createElement("graph")
    graph_node.setAttribute("xmlns", "http://www.cs.rpi.edu/XGMML")
    graph_node.setAttribute("xmlns:dc", "http://purl.org/dc/elements/1.1/")
    graph_node.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink")
    graph_node.setAttribute("xmlns:rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#")
    graph_node.setAttribute("xmlns:cy", "http://www.cytoscape.org")
    graph_node.setAttribute("directed", "1")
    doc.appendChild(graph_node)
    name_map = {}
    for i, n in enumerate(networkx_graph.node):
        name_map[n] = i
        node = doc.createElement("node")
        node.setAttribute("label", str(n))
        node.setAttribute("id", str(i))
        ## set type string
        att_node = doc.createElement("att")
        att_node.setAttribute("name", "TYPE")
        att_node.setAttribute("value", str(networkx_graph.node[n]["type"]))
        att_node.setAttribute("type", "string")
        node.appendChild(att_node)
        ## set label string
        att_node = doc.createElement("att")
        att_node.setAttribute("name", "LABEL")
        if networkx_graph.node[n]["type"] == "protein":
            att_node.setAttribute("value", str(n))
        else:
            att_node.setAttribute("value", "")
        att_node.setAttribute("type", "string")
        node.appendChild(att_node)
        ## set signature values
        att_node = doc.createElement("att")
        att_node.setAttribute("name", "SIGNATURE")
        if str(n) in signature_map:
            att_node.setAttribute("value", str(signature_map[str(n)]))
        else:
            att_node.setAttribute("value", "0")
        att_node.setAttribute("type", "real")
        node.appendChild(att_node)
        ## set score values
        att_node = doc.createElement("att")
        att_node.setAttribute("name", "SCORE")
        if str(n) in score_map:
            att_node.setAttribute("value", str(score_map[str(n)]))
        else:
            att_node.setAttribute("value", "0")
        att_node.setAttribute("type", "real")
        node.appendChild(att_node)
        ## set bootstrap values
        att_node = doc.createElement("att")
        att_node.setAttribute("name", "BOOTSTRAP")
        if str(n) in bootstrap_map:
            att_node.setAttribute("value", str(bootstrap_map[str(n)]))
        else:
            att_node.setAttribute("value", "0")
        att_node.setAttribute("type", "real")
        node.appendChild(att_node)
        graph_node.appendChild(node)

    for source in networkx_graph.edge:
        for target in networkx_graph.edge[source]:
            for edge in networkx_graph.edge[source][target]:
                edge_node = doc.createElement("edge")
                edge_node.setAttribute("label", "%s - %s" % (source, target))
                edge_node.setAttribute("source", str(name_map[source]))
                edge_node.setAttribute("target", str(name_map[target]))
                ## set edge attributes
                for key, value in networkx_graph.edge[source][target][edge].items():
                    att_node = doc.createElement("att")
                    att_node.setAttribute("name", key)
                    att_node.setAttribute("value", str(value))
                    if type(value) == float:
                        att_node.setAttribute("type", "real")
                    elif type(value) == int:
                        att_node.setAttribute("type", "integer")
                    else:
                        att_node.setAttribute("type", "string")
                    edge_node.appendChild(att_node)
                ## set bootstrap values
                att_node = doc.createElement("att")
                att_node.setAttribute("name", "BOOTSTRAP")
                if (str(source), str(target)) in bootstrap_map:
                    att_node.setAttribute("value", str(bootstrap_map[(str(source), str(target))]))
                else:
                    att_node.setAttribute("value", "0")
                att_node.setAttribute("type", "real")
                edge_node.appendChild(att_node)
                graph_node.appendChild(edge_node)
    o = open(output_file, "w")
    doc.writexml(o, addindent = " ", newl = "\n")
    o.close()

## jt classes
class branchPATHMARK(Target):
    def __init__(self, signature_file, null_file, bootstrap_file, global_pathway, parameters, directory):
        Target.__init__(self, time=10000)
        self.signature_file = signature_file
        self.null_file = null_file
        self.bootstrap_file = bootstrap_file
        self.global_pathway = global_pathway
        self.parameters = parameters
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        if not os.path.exists("analysis"):
            os.mkdir("analysis")
        
        ## convert global_pathway to global_graph
        global_graph = self.global_pathway.networkx()
        
        ## read input files
        real_frame = pandas.read_csv(self.signature_file, sep = "\t", index_col = 0)
        real_signatures = list(real_frame.columns)
        null_frame = None
        null_signatures = []
        if self.null_file:
            null_frame = pandas.read_csv(self.null_file, sep = "\t", index_col = 0)
            null_signatures = list(null_frame.columns)
        bootstrap_frame = None
        bootstrap_signatures = []
        if self.bootstrap_file:
            bootstrap_frame = pandas.read_csv(self.bootstrap_file, sep = "\t", index_col = 0)
            bootstrap_signatures = list(bootstrap_frame.columns)
        
        ## iterate over signatures
        signature_map = {}
        score_map = {}
        for signature in real_signatures:
            assert(not signature.startswith("_"))
            signature_directory = "analysis/%s" % (signature)
            signature_sif = "%s_%s.sif" % (signature, self.parameters.parameters_string)
            if not os.path.exists(signature_directory):
                os.mkdir(signature_directory)
            else:
                continue
            ## generate real PATHMARK
            signature_frame = real_frame[[signature]]
            if self.parameters.heat_diffusion_time > 0.0:
                score_frame = diffuseHeat(signature_frame,
                                          global_graph,
                                          self.parameters)
            else:
                score_frame = abs(signature_frame)
            score_statistics = computeStatistics(score_frame)
            signature_graph = getPathmark(score_frame,
                                          global_graph,
                                          self.parameters,
                                          forced_statistics = score_statistics)
            signature_graph = filterGeneGroups(signature_graph, global_graph)
            if self.parameters.hub_filter:
                node_threshold = score_statistics[0] + self.parameters.filter_parameters[0]*score_statistics[1]
                signature_graph = filterGoodHub(signature_graph, global_graph, score_frame.icol(0), node_threshold)
            writeNetworkXSIF("%s/%s" % (signature_directory, signature_sif), signature_graph)
            signature_map[signature] = deepcopy(signature_frame.icol(0))
            score_map[signature] = deepcopy(score_frame.icol(0))
            ## generate bootstrap PATHMARKs
            bootstraps = filter(lambda x: x.startswith(signature), bootstrap_signatures)
            for bootstrap in bootstraps:
                bootstrap_index = bootstrap.split(":")[-1]
                self.addChildTarget(queuePATHMARK("%s/_b%s_%s" % (signature_directory, bootstrap_index, signature_sif),
                                                  bootstrap_frame[[bootstrap]],
                                                  global_graph,
                                                  self.parameters,
                                                  self.directory,
                                                  forced_statistics = None))
            ## generate null PATHMARKs
            nulls = filter(lambda x: x.startswith(signature), null_signatures)
            for null in nulls:
                null_index = null.split(":")[-1]
                self.addChildTarget(queuePATHMARK("%s/_n%s_%s" % (signature_directory, null_index, signature_sif),
                                                  null_frame[[null]],
                                                  global_graph,
                                                  self.parameters,
                                                  self.directory,
                                                  forced_statistics = score_statistics))
        self.setFollowOnTarget(analyzeGraphStatistics(real_signatures,
                                                      signature_map,
                                                      score_map,
                                                      self.global_pathway,
                                                      self.parameters,
                                                      self.directory))

class queuePATHMARK(Target):
    def __init__(self, output_sif, signature_frame, global_graph, parameters, directory, forced_statistics = None):
        Target.__init__(self, time=10000)
        self.output_sif = output_sif
        self.signature_frame = signature_frame
        self.global_graph = global_graph
        self.parameters = parameters
        self.directory = directory
        self.forced_statistics = forced_statistics
    def run(self):
        os.chdir(self.directory)
        
        if self.parameters.heat_diffusion_time > 0.0:
            score_frame = diffuseHeat(self.signature_frame,
                                      self.global_graph,
                                      self.parameters)
        else:
            score_frame = abs(self.signature_frame)
        score_statistics = computeStatistics(score_frame)
        signature_graph = getPathmark(score_frame,
                                      self.global_graph,
                                      self.parameters,
                                      forced_statistics = self.forced_statistics)
        signature_graph = filterGeneGroups(signature_graph, self.global_graph)
        if self.parameters.hub_filter:
            node_threshold = score_statistics[0] + self.parameters.filter_parameters[0]*score_statistics[1]
            signature_graph = filterGoodHub(signature_graph, self.global_graph, score_frame.icol(0), node_threshold)
        writeNetworkXSIF(self.output_sif, signature_graph)

class analyzeGraphStatistics(Target):
    def __init__(self, signatures, signature_map, score_map, global_pathway, parameters, directory):
        Target.__init__(self, time=10000)
        self.signatures = signatures
        self.signature_map = signature_map
        self.score_map = score_map
        self.global_pathway = global_pathway
        self.parameters = parameters
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## iterate over signatures
        graph_map = {}
        bootstrap_map = {}
        null_map = {}
        for signature in self.signatures:
            ## store real_graph
            real_file = "%s_%s.sif" % (signature, self.parameters.parameters_string)
            real_graph = readNetworkXSIF("analysis/%s/%s" % (signature, real_file), node_type_map = self.global_pathway.nodes)
            graph_map[signature] = deepcopy(real_graph)
            ## compute node_occurrences
            bootstrap_counts = {}
            bootstrap_occurences = {}
            bootstrap_files = filter(lambda x: x.startswith("_b"), os.listdir("analysis/%s" % (signature)))
            for bootstrap_file in bootstrap_files:
                bootstrap_graph = readNetworkXSIF("analysis/%s/%s" % (signature, bootstrap_file), node_type_map = self.global_pathway.nodes)
                bootstrap_nodes = set(bootstrap_graph.node)
                for node in bootstrap_nodes:
                    if node not in bootstrap_counts:
                        bootstrap_counts[node] = 1
                    else:
                        bootstrap_counts[node] += 1
                bootstrap_edges = set([(i[1], i[0]) for i in set(bootstrap_graph.edges_iter())]) | set(bootstrap_graph.edges_iter())
                for edge in bootstrap_edges:
                    if edge not in bootstrap_counts:
                        bootstrap_counts[edge] = 1
                    else:
                        bootstrap_counts[edge] += 1
            for element in bootstrap_counts:
                bootstrap_occurences[element] = float(bootstrap_counts[element])/float(len(bootstrap_files))
            if len(bootstrap_files) > 0:
                bootstrap_map[signature] = deepcopy(bootstrap_occurences)
            else:
                bootstrap_map[signature] = {}
            ## compute null_statistics
            null_statistics = []
            null_files = filter(lambda x: x.startswith("_n"), os.listdir("analysis/%s" % (signature)))
            for null_file in null_files:
                null_graph = readNetworkXSIF("analysis/%s/%s" % (signature, null_file), node_type_map = self.global_pathway.nodes)
                null_nodes = set(null_graph.node)
                null_node_count = len(null_nodes)
                null_edges = set([(i[1], i[0]) for i in set(null_graph.edges_iter())]) | set(null_graph.edges_iter())
                null_edge_count = len(null_edges)/2
                largest_null_subgraph = getLargestSubgraph(null_graph)
                largest_null_nodes = set(largest_null_subgraph.node)
                largest_null_node_count = len(largest_null_nodes)
                largest_null_edges = set([(i[1], i[0]) for i in set(largest_null_subgraph.edges_iter())]) | set(largest_null_subgraph.edges_iter())
                largest_null_edge_count = len(largest_null_edges)/2
                null_statistics.append( (null_node_count, null_edge_count, largest_null_node_count, largest_null_edge_count) )
            if len(null_files) > 0:
                null_map[signature] = deepcopy(zip(*null_statistics))
            else:
                null_map[signature] = None
        self.setFollowOnTarget(generateOutput(self.signatures,
                                              graph_map,
                                              self.signature_map,
                                              self.score_map,
                                              bootstrap_map,
                                              null_map,
                                              self.parameters,
                                              self.directory))

class generateOutput(Target):
    def __init__(self, signatures, graph_map, signature_map, score_map, bootstrap_map, null_map, parameters, directory):
        Target.__init__(self, time=10000)
        self.signatures = signatures
        self.graph_map = graph_map
        self.signature_map = signature_map
        self.score_map = score_map
        self.bootstrap_map = bootstrap_map
        self.null_map = null_map
        self.parameters = parameters
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        if not os.path.exists("report"):
            os.mkdir("report")
        
        ## iterate over signatures
        for signature in self.signatures:
            ## output xgmml, sif, and table
            writeNetworkXCytoscapeXGMML("report/%s_%s.xgmml" % (signature, self.parameters.parameters_string),
                                        self.graph_map[signature],
                                        self.signature_map[signature],
                                        self.score_map[signature],
                                        self.bootstrap_map[signature])
            writeNetworkXSIF("report/%s_%s.sif" % (signature, self.parameters.parameters_string),
                             self.graph_map[signature])
            if len(self.bootstrap_map[signature]) > 0:
                o = open("report/%s_%s.bootstrap.tab" % (signature, self.parameters.parameters_string), "w")
                for node in self.bootstrap_map[signature]:
                    o.write("%s\t%s\n" % (node, self.bootstrap_map[signature][node]))
                o.close()
                pandas.DataFrame([self.signature_map[signature], self.score_map[signature]],
                                 index = ["Signature","Score"]).transpose().to_csv("report/%s_%s.tab" % (signature, self.parameters.parameters_string), sep = "\t", na_rep = "nan")
            else:
                pandas.DataFrame([self.signature_map[signature], self.score_map[signature]],
                                 index = ["Signature","Score"]).transpose().to_csv("report/%s_%s.tab" % (signature, self.parameters.parameters_string), sep = "\t", na_rep = "nan")
            
            ## output significance plots
            if self.null_map[signature] is not None:
                real_graph = self.graph_map[signature]
                real_nodes = set(real_graph.node)
                real_node_count = len(real_nodes)
                real_edges = set([(i[1], i[0]) for i in set(real_graph.edges_iter())]) | set(real_graph.edges_iter())
                real_edge_count = len(real_edges)/2
                largest_real_subgraph = getLargestSubgraph(real_graph)
                largest_real_nodes = set(largest_real_subgraph.node)
                largest_real_node_count = len(largest_real_nodes)
                largest_real_edges = set([(i[1], i[0]) for i in set(largest_real_subgraph.edges_iter())]) | set(largest_real_subgraph.edges_iter())
                largest_real_edge_count = len(largest_real_edges)/2
                real_statistics = [real_node_count, real_edge_count, largest_real_node_count, largest_real_edge_count]
                for index in range(4):
                    xmin = max(0, min([real_statistics[index]] + list(self.null_map[signature][index])) - 50)
                    xmax = max([real_statistics[index]] + list(self.null_map[signature][index])) + 50
                    (n, bins, patches) = plt.hist(list(self.null_map[signature][index]), bins = 20, range = (xmin, xmax), histtype = "stepfilled")
                    plt.axvline(x = real_statistics[index], linestyle = "--", linewidth = 2, color = "r")
                    plt.savefig("report/%s_%s_%s.png" % (signature, self.parameters.parameters_string, index))
                    plt.clf()

def main():
    ## check for fresh run
    if os.path.exists(".jobTree"):
        logging.warning("WARNING: '.jobTree' directory found, remove it first to start a fresh run\n")
    
    ## parse arguments
    parser = OptionParser(usage = "%prog [options] signature_file pathway_file")
    Stack.addJobTreeOptions(parser)
    parser.add_option("--jobFile",
                      help = "Add as child of jobFile rather than new jobTree")
    parser.add_option("-b", "--bootstrap", dest = "bootstrap_file", default = None,
                      help = "")
    parser.add_option("-n", "--null", dest = "null_file", default = None,
                      help = "")
    parser.add_option("-f", "--filter", dest = "filter_parameters", default = "0.0;0.0",
                      help = "")
    parser.add_option("-t", "--heat", dest = "heat_diffusion", default = "0",
                      help = "")
    parser.add_option("-u", "--hub", dest = "hub_filter", action = "store_true", default = False,
                      help = "")
    options, args = parser.parse_args()
    logging.info("options: %s" % (str(options)))
    print "Using Batch System '%s'" % (options.batchSystem)
    
    if len(args) != 2:
        logging.error("ERROR: incorrect number of arguments\n")
        sys.exit(1)
    signature_file = os.path.abspath(args[0])
    pathway_file = os.path.abspath(args[1])
    
    ## set pathway
    global_pathway = Pathway(pathway_file)
    
    ## set parameters
    parameters = Parameters(pathway_file,
                            options.filter_parameters,
                            options.heat_diffusion,
                            options.hub_filter)
    
    ## run
    s = Stack(branchPATHMARK(signature_file,
                             options.null_file,
                             options.bootstrap_file,
                             global_pathway,
                             parameters,
                             os.getcwd().rstrip("/")))
    if options.jobFile:
        s.addToJobFile(options.jobFile)
    else:
        if options.jobTree == None:
            options.jobTree = "./.jobTree"
        
        jobtree_dir = options.jobTree
        lasttree_dir = jobtree_dir.replace("job", "last")
        assert(jobtree_dir != lasttree_dir)
        
        failed = s.startJobTree(options)
        if failed:
            logging.warning("WARNING: %d jobs failed" % (failed))
        else:
            logging.info("Run complete!")
            if os.path.exists(lasttree_dir):
                shutil.rmtree(lasttree_dir)
            if os.path.exists(jobtree_dir):
                shutil.move(jobtree_dir, lasttree_dir)

if __name__ == "__main__":
    from PATHMARK import *
    main()
