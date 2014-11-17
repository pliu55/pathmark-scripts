#!/usr/bin/env python
"""
shazamCytoscapeWeb.py
    by Sam Ng
"""
import math, os, random, re, shutil, sys, types
from copy import deepcopy

import numpy
import pandas
from PATHMARK import *
from signature import *

from optparse import OptionParser

def generatePathways(summary_list, pathway_directory, name_map):
    for pathway_name in summary_list:
        if pathway_name in name_map:
            pathway_pid = name_map[pathway_name]
            matches = filter(lambda x: pathway_pid == x.lstrip("pid_").rstrip("_pathway.tab"), os.listdir(pathway_directory))
            if len(matches) == 1:
                pathway_file = matches[0]
                yield(pathway_name, pathway_pid, pathway_file)

def filterDataFrame(data_frame, pathway_pid, pathway_features):
    features = filter(lambda x : x.startswith(pathway_pid + "_"), data_frame.index)
    if len(features) == 0:
        features = list(set(data_frame.index) & set(pathway_features))
        filtered_frame = data_frame.loc[features].copy()
    else:
        filtered_frame = data_frame.loc[features].copy()
        filtered_frame.index = [i.lstrip(pathway_pid + "_") for i in features]
    return(filtered_frame)

def shazamCytoscapeWeb():
    ## parse arguments
    parser = OptionParser(usage = "%prog [options]")
    parser.add_option("-t", "--summary", dest = "summary_tsv", default = None,
                      help = "")
    parser.add_option("-n", "--visualize", type = "int", dest = "visualize_n", default = 5,
                      help = "")
    parser.add_option("-d", "--data", dest = "data_files", default = None,
                      help = "")
    parser.add_option("-s", "--signature", type = "int", dest = "signature_index", default = 0,
                      help = "")
    parser.add_option("-p", "--phenotype", dest = "phenotype_file", default = None,
                      help = "")
    parser.add_option("-i", "--pathway", dest = "pathway_source", default = None,
                      help = "")
    parser.add_option("-c", "--circle", dest = "plot_circle", action = "store_true", default = False,
                      help = "")
    parser.add_option("-l", "--legend", dest = "legend_file", default = None,
                      help = "")
    parser.add_option("-o", "--output", dest = "output_directory", default = None,
                      help = "")
    options, args = parser.parse_args()
    
    assert(len(args) == 0)
    
    ## create signature_files directory
    if not os.path.exists("analysis"):
        os.mkdir("analysis")
    assert(not os.path.exists("analysis/layout_files"))
    os.mkdir("analysis/layout_files")
    
    ## set pathway_directory and name_map
    assert(options.pathway_source is not None)
    assert(os.path.exists(options.pathway_source))
    if options.pathway_source.endswith(".zip"):
        pathway_directory = options.pathway_source.rstrip(".zip")
    else:
        pathway_directory = options.pathway_source
    assert(os.path.exists("%s/names.tab" % (pathway_directory)))
    name_map = {}
    f = open("%s/names.tab" % (pathway_directory), "r")
    f.readline()
    for line in f:
        if line.isspace():
            continue
        pline = line.rstrip().split("\t")
        name_map[pline[1]] = pline[0]
    f.close()
    
    ## read summary_tsv for the list of pathways to plot and pick the top visualize_n
    summary_map = {}
    summary_list = []
    f = open(options.summary_tsv, "r")
    f.readline()
    for line in f:
        if line.isspace():
            continue
        pline = line.rstrip().split("\t")
        try:
            fvalue = float(pline[8])
            summary_map[pline[0]] = fvalue
        except ValueError:
            pass
    summary_features = summary_map.keys()
    summary_features.sort(lambda x, y : cmp(summary_map[y], summary_map[x]))
    for index in range(options.visualize_n):
        if index >= len(summary_features):
            break
        summary_list.append(summary_features[index])
    
    ## create signatures
    data_files = options.data_files.split(",")
    assert(len(data_files) - 1 >= options.signature_index)
    assert(options.signature_index >= 0)
    data_frames = {}
    for data_file in data_files:
        data_frames[data_file] = pandas.read_csv(data_file, sep = "\t", index_col = 0)
    signature_file = data_files[options.signature_index]
    signature_samples = list(data_frames[signature_file].columns)
    phenotype_frame = pandas.read_csv(options.phenotype_file, sep = "\t", index_col = 0)
    phenotypes = list(phenotype_frame.columns)
    for phenotype in phenotypes:
        assert(not phenotype.startswith("_"))
        for signature_name, positive_samples, negative_samples in generateDichotomies(phenotype, phenotype_frame, signature_samples):
            siggenes = rpy2SignatureGenes()
            signature_frame = siggenes.calculate("sam", data_frames[signature_file], positive_samples, negative_samples)[["Score"]]
            signature_frame.columns = [signature_name]
            for pathway_name, pathway_pid, pathway_file in generatePathways(summary_list, pathway_directory, name_map):
                pathway_features = Pathway("%s/%s" % (pathway_directory, pathway_file)).nodes.keys()
                if os.path.exists("%s/pid_%s_%s" % (options.output_directory, pathway_pid, signature_name)):
                    break
                os.mkdir("analysis/layout_files/pid_%s_%s" % (pathway_pid, signature_name))
                ## generate signature file
                filtered_signature_frame = filterDataFrame(signature_frame, pathway_pid, pathway_features)
                filtered_signature_frame.to_csv("analysis/layout_files/pid_%s_%s/%s.tab" % (pathway_pid, signature_name, signature_name), sep = "\t", index_label = "id")
                ## generate circleplots and layout
                if options.plot_circle:
                    l = open("analysis/layout_files/pid_%s_%s/os.log" % (pathway_pid, signature_name), "w")
                    os.mkdir("analysis/layout_files/pid_%s_%s/img" % (pathway_pid, signature_name))
                    o = open("analysis/layout_files/pid_%s_%s/include.samples" % (pathway_pid, signature_name), "w")
                    o.write("%s\n" % ("\n".join(positive_samples + negative_samples)))
                    o.close()
                    o = open("analysis/layout_files/pid_%s_%s/include.features" % (pathway_pid, signature_name), "w")
                    o.write("%s\n" % ("\n".join(pathway_features + ["GLOBAL"])))
                    o.close()
                    o = open("analysis/layout_files/pid_%s_%s/signature.circle" % (pathway_pid, signature_name), "w")
                    o.write("id\t%s\n" % ("\t".join(positive_samples + negative_samples)))
                    o.write("*\t%s\n" % ("\t".join(["1" for sample in positive_samples] + ["0" for sample in negative_samples])))
                    o.close()
                    o = open("analysis/layout_files/pid_%s_%s/color.map" % (pathway_pid, signature_name), "w")
                    o.write("> 1\n0\t255.255.255\n1\t0.0.0\n")
                    o.close()
                    for data_file in data_files:
                        filtered_data_frame = filterDataFrame(data_frames[data_file], pathway_pid, pathway_features)
                        if "*" in filtered_data_frame.index:
                            output_data_frame = filtered_data_frame
                        else:
                            global_data_frame = pandas.DataFrame([{sample : sum([abs(filtered_signature_frame[signature_name][feature])*filtered_data_frame[sample][feature] for feature in list(set(filtered_data_frame.index) & set(filtered_signature_frame.index) - set(filtered_signature_frame.index[filtered_signature_frame[signature_name].apply(numpy.isnan)]))]) for sample in filtered_data_frame.columns}], index = ["GLOBAL"])
                            output_data_frame = filtered_data_frame.append(global_data_frame)
                        output_data_frame.to_csv("analysis/layout_files/pid_%s_%s/%s" % (pathway_pid, signature_name, data_file.split("/")[-1]), sep = "\t", index_label = "id")
                    os.chdir("analysis/layout_files/pid_%s_%s" % (pathway_pid, signature_name))
                    sort_files = ["signature.circle", data_files[-1].split("/")[-1]]
                    cmd = "circlePlot.py -s include.samples -f include.features -m color.map -o \"%s;%s\" img/ signature.circle %s" % ("GLOBAL", ",".join(sort_files), " ".join([data_file.split("/")[-1] for data_file in data_files]))
                    l.write(cmd + "\n")
                    os.system(cmd)
                    os.chdir("../../..")
                    cmd = "layoutCytoscapeWeb.py -s %s -p %s -i %s -o %s" % ("analysis/layout_files/pid_%s_%s/%s.tab" % (pathway_pid, signature_name, signature_name), "%s/%s" % (pathway_directory, pathway_file), "analysis/layout_files/pid_%s_%s/img/" % (pathway_pid, signature_name), options.output_directory)
                    if options.legend_file:
                        cmd += " -l %s" % (options.legend_file)
                    l.write(cmd + "\n")
                    os.system(cmd)
                    l.close()
                else:
                    cmd = "python layoutCytoscapeWeb.py -s %s -p %s -o %s" % ("analysis/layout_files/pid_%s_%s/%s.tab" % (pathway_pid, signature_name, signature_name), "%s/%s" % (pathway_directory, pathway_file), options.output_directory)
                    os.system(cmd)
    # shutil.rmtree("analysis/layout_files")

if __name__ == "__main__":
    shazamCytoscapeWeb()
