#!/usr/bin/env python
"""
signature.py
    by Sam Ng
"""
import logging, math, os, random, re, shutil, sys, types
from copy import deepcopy

import pandas
import pandas.rpy.common as common
from rpy2 import robjects
from rpy2.robjects.packages import importr
from scipy import stats

from optparse import OptionParser
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

#### NOTE BLOCK
#### - Convert LIMMA code to rpy2
#### - Add method for dealing with either median or mean dichotomies of continuous data

## logger
logging.basicConfig(filename = "signature.log", level = logging.INFO)

## default variables
null_prefixes = ["na_", "nw_"]          ## valid null samples must start with one of these
                                        ## prefixes and end in a sample name

## pm classes
class Parameters:
    """
    Stores parameters used for this [signature.py specific]
    """
    def __init__(self, random_seed = None, bootstrap_size = 10, null_method = "paradigm", null_size = 0, signature_method = "sam"):
        if random_seed is None:
            self.random_seed = random.randint(0, 999999999)
        else:
            if os.path.exists(random_seed):
                f = open(random_seed, 'r')
                self.random_seed = int(f.readline().rstrip())
                f.close()
            else:
                self.random_seed = int(random_seed)
        self.bootstrap_size = bootstrap_size
        self.null_method = null_method
        self.null_size = null_size
        self.signature_method = signature_method
        self.signature_file = "signature.tab"

class rpy2SignatureGenes:
    def __init__(self):
        self.siggenes = importr("siggenes")
    def calculate(self, method, data_frame, positive_samples, negative_samples):
        ## construct matrix_r
        r = robjects.r
        samples = []
        for sample in data_frame.axes[1]:
            if sample in positive_samples + negative_samples:
                samples.append(sample)
        samples.sort()
        features = list(data_frame.axes[0])
        features.sort()
        matrix = data_frame[samples].loc[features]
        matrix_r = common.convert_to_r_dataframe(matrix)
        ## construct cls_r
        cls = {}
        for sample in samples:
            if sample in positive_samples:
                cls[sample] = 1
            elif sample in negative_samples:
                cls[sample] = 0
        cls_r = common.convert_to_r_matrix( pandas.DataFrame( [cls] ) )
        ## generate signature with method
        sam_out = self.siggenes.sam(matrix_r, r.c(cls_r))
        sam_att = r.cbind(
            r.c(r.attributes(sam_out).rx2("d")),
            r.c(r.attributes(sam_out).rx2("vec.false")),
            r.c(r.attributes(sam_out).rx2("q.value")),
            r.c(r.attributes(sam_out).rx2("p.value")),
            r.c(r.attributes(sam_out).rx2("s"))
        )
        ## return results as a data_frame
        ocols = ["Score", "FalseCalls", "Q-value", "P-value", "StdDev"]
        output = {}
        for j, col in enumerate(ocols):
            row = {}
            for i, n in enumerate(features):
                # print n, col, sam_att.rx(i + 1, j + 1)[0]
                row[n] = sam_att.rx(i + 1, j + 1)[0]
            # print row
            output[col] = row
        return(pandas.DataFrame(output))

class scipySignatureGenes:
    def __init__(self):
        pass
    def calculate(self, method, data_frame, positive_samples, negative_samples):
        t_map = {}
        for feature, values in data_frame.iterrows():
            t_map[feature] = stats.ttest_ind(values[positive_samples], values[negative_samples], equal_var = False)[0]
        return(pandas.DataFrame(pandas.Series(t_map, name = "T-statistic")))

class rpy2LIMMA:
    def __init__(self):
        self.edger = importr("edgeR")
        self.limma = importr("limma")
    def calculate(self, data_frame, positive_samples, negative_samples, output_file):
        ## construct matrix_r
        # r = robjects.r
        samples = []
        for sample in data_frame.axes[1]:
            if sample in positive_samples + negative_samples:
                samples.append(sample)
        samples.sort()
        features = list(data_frame.axes[0])
        features.sort()
        matrix = data_frame[samples].loc[features]
        # matrix_r = common.convert_to_r_matrix(matrix)
        matrix_file = "%s.input" % (output_file)
        matrix.to_csv(matrix_file, sep = "\t", index_label = "id")
        
        ## construct cls_r
        # cls = []
        # for sample in samples:
        #     if sample in positive_samples:
        #         cls.append(1)
        #     elif sample in negative_samples:
        #         cls.append(0)
        # cls_r = robjects.IntVector(cls)
        cls_file = "%s.contrast" % (output_file)
        o = open(cls_file, "w")
        o.write("sample\tcluster\n")
        for sample in samples:
            if sample in positive_samples:
                o.write("%s\t+\n" % (sample))
            elif sample in negative_samples:
                o.write("%s\t-\n" % (sample))
        o.close()
        
        ## generate signature with limma
        # dge = self.edger.DGEList(counts=matrix_r)
        
        # isexpr = r.rowSums(edger.cpm(dge) > 10) >= 2 #### comparison not right for R
        # flt = dge[isexpr,]
        # tmm = calcNormFactors(flt)
        # design = model.matrix(~ contrast)
        # y = voom(tmm, design, plot=FALSE)
        # fit = eBayes(lmFit(y, design))
        # cn = sprintf("contrast%s", as.character(levels(contrast)[2]))
        # tt = topTable(fit, coef=cn, number=200)
        # limma_out = sol = list(design=design, y=y, fit=fit, tt=tt)
        
        # sam_att = r.cbind(
        #     r.c(r.attributes(sam_out).rx2("d")),
        #     r.c(r.attributes(sam_out).rx2("vec.false")),
        #     r.c(r.attributes(sam_out).rx2("q.value")),
        #     r.c(r.attributes(sam_out).rx2("p.value")),
        #     r.c(r.attributes(sam_out).rx2("s"))
        # )
        
        ## return results as a data_frame
        # ocols = ["Score", "FalseCalls", "Q-value", "P-value", "StdDev"]
        # output = {}
        # for j, col in enumerate(ocols):
        #     row = {}
        #     for i, n in enumerate(features):
        #         # print n, col, sam_att.rx(i + 1, j + 1)[0]
        #         row[n] = sam_att.rx(i + 1, j + 1)[0]
        #     # print row
        #     output[col] = row
        # return(pandas.DataFrame(output))
        output_name = output_file.split("/")[-1]
        result_file = "%s.result" % (output_file)
        top_file = "%s.top_100.tab" % (output_name)
        plot_file = "%s.plot.pdf" % (output_name)
        cmd = "Rscript limma_ng.R %s %s 100 %s %s %s" % (matrix_file, cls_file, result_file, top_file, plot_file)
        os.system(cmd)
        result_frame = pandas.read_csv(result_file, sep = "\t", index_col = 0)
        os.remove(matrix_file)
        os.remove(cls_file)
        os.remove(result_file)
        if output_name.startswith("_n") or output_name.startswith("_b"):
            os.remove(top_file)
            os.remove(plot_file)
        return(pandas.DataFrame(result_frame.icol(1)))
        
## signature functions
def generateDichotomies(phenotype, phenotype_frame, data_samples, reverse_sort = False):
    """
    Generates dichotomies for a given column in the phenotype file [signature.py specific]
    """
    phenotype_map = phenotype_frame[phenotype]
    phenotype_is_categorical = True
    for element in set(phenotype_map.values):
        try:
            float_element = float(element)
            if not float_element.is_integer():
                phenotype_is_categorical = False
        except ValueError:
            pass
    if phenotype_is_categorical:
        phenotype_categories = list(set(phenotype_map.values))
        phenotype_categories.sort(reverse = reverse_sort)
        for index in range(len(phenotype_categories)):
            if index == 1 and len(phenotype_categories) == 2:
                break
            signature_name = "%s=%s_SIGNATURE" % (phenotype, phenotype_categories[index])
            positive_samples = list(set(phenotype_map[phenotype_map == phenotype_categories[index]].index) & set(data_samples))
            negative_samples = list(set(phenotype_map[phenotype_map != phenotype_categories[index]].index) & set(data_samples))
            yield(signature_name, positive_samples, negative_samples)
    #### add method for dealing with either median or mean splits of dichotomies
    #### what if we want to split by median or mean for integer counts?

## jt classes
class branchSignatures(Target):
    def __init__(self, data_file, phenotype_file, parameters, directory):
        Target.__init__(self, time=10000)
        self.data_file = data_file
        self.phenotype_file = phenotype_file
        self.parameters = parameters
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        if not os.path.exists("analysis"):
            os.mkdir("analysis")
        assert(not os.path.exists("analysis/signature_files"))
        os.mkdir("analysis/signature_files")
        
        ## read in data and phenotype matrices
        data_frame = pandas.read_csv(self.data_file, sep = "\t", index_col = 0)
        data_samples = list(data_frame.columns)
        phenotype_frame = pandas.read_csv(self.phenotype_file, sep = "\t", index_col = 0)
        phenotypes = list(phenotype_frame.columns)
        for phenotype in phenotypes:
            assert(not phenotype.startswith("_"))
            for signature_name, positive_samples, negative_samples in generateDichotomies(phenotype, phenotype_frame, data_samples):
                self.addChildTarget(generateBatches(data_frame,
                                                    signature_name,
                                                    positive_samples,
                                                    negative_samples,
                                                    self.parameters,
                                                    self.directory))
        self.setFollowOnTarget(mergeAllSignatures(self.parameters,
                                                  self.directory))

class generateBatches(Target):
    def __init__(self, data_frame, signature_name, positive_samples, negative_samples, parameters, directory):
        Target.__init__(self, time=10000)
        self.data_frame = data_frame
        self.signature_name = signature_name
        self.positive_samples = positive_samples
        self.negative_samples = negative_samples
        self.parameters = parameters
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        random.seed(self.parameters.random_seed)
        
        ## queue real
        self.addChildTarget(computeSignatures("analysis/signature_files/%s" % (self.signature_name),
                                              self.data_frame,
                                              self.signature_name,
                                              self.positive_samples,
                                              self.negative_samples,
                                              self.parameters,
                                              self.directory))
        
        ## generate and queue bootstraps
        for bootstrap_index in range(self.parameters.bootstrap_size):
            bootstrap_positive_data_frame = pandas.DataFrame(index = self.data_frame.index)
            for sample in self.positive_samples:
                bootstrap_positive_data_frame[sample] = self.data_frame[random.sample(self.positive_samples, 1)[0]]
            bootstrap_negative_data_frame = pandas.DataFrame(index = self.data_frame.index)
            for sample in self.negative_samples:
                bootstrap_negative_data_frame[sample] = self.data_frame[random.sample(self.negative_samples, 1)[0]]
            bootstrap_data_frame = bootstrap_positive_data_frame.join(bootstrap_negative_data_frame, how = "outer")
            self.addChildTarget(computeSignatures("analysis/signature_files/_b%s_%s" % (bootstrap_index + 1, self.signature_name),
                                                  bootstrap_data_frame,
                                                  "%s:%s" % (self.signature_name, bootstrap_index + 1),
                                                  self.positive_samples,
                                                  self.negative_samples,
                                                  self.parameters,
                                                  self.directory))
        
        ## generate and queue nulls
        for null_index in range(self.parameters.null_size):
            if self.parameters.null_method == "samples":
                null_data_frame = self.data_frame[self.positive_samples + self.negative_samples].copy()
                null_data_frame.columns = random.sample(null_data_frame.columns, len(null_data_frame.columns))
            elif self.parameters.null_method == "features":
                null_data_frame = self.data_frame[self.positive_samples + self.negative_samples].copy()
                null_data_frame.index = random.sample(null_data_frame.index, len(null_data_frame.index))
            elif self.parameters.null_method == "paradigm":
                null_data_frame = pandas.DataFrame(index = self.data_frame.index)
                for sample in self.positive_samples + self.negative_samples:
                    null_samples = filter(lambda x: x.endswith(sample), self.data_frame.columns)
                    null_samples = filter(lambda x: filter(x.startswith, null_prefixes), null_samples)
                    null_data_frame[sample] = self.data_frame[random.sample(null_samples, 1)[0]]
            self.addChildTarget(computeSignatures("analysis/signature_files/_n%s_%s" % (null_index + 1, self.signature_name),
                                                  null_data_frame,
                                                  "%s:%s" % (self.signature_name, null_index + 1),
                                                  self.positive_samples,
                                                  self.negative_samples,
                                                  self.parameters,
                                                  self.directory))

class computeSignatures(Target):
    def __init__(self, output_file, data_frame, signature_name, positive_samples, negative_samples, parameters, directory):
        Target.__init__(self, time=10000)
        self.output_file = output_file
        self.data_frame = data_frame
        self.signature_name = signature_name
        self.positive_samples = positive_samples
        self.negative_samples = negative_samples
        self.parameters = parameters
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        if self.parameters.signature_method == "sam":
            siggenes = rpy2SignatureGenes()
            signature_frame = siggenes.calculate(self.parameters.signature_method, self.data_frame, self.positive_samples, self.negative_samples)[["Score"]]
            signature_frame.columns = [self.signature_name]
            signature_frame.to_csv(self.output_file, sep = "\t", index_label = "id")
        elif self.parameters.signature_method == "ttest":
            scipy = scipySignatureGenes()
            signature_frame = scipy.calculate(self.parameters.signature_method, self.data_frame, self.positive_samples, self.negative_samples)
            signature_frame.columns = [self.signature_name]
            signature_frame.to_csv(self.output_file, sep = "\t", index_label = "id")
        elif self.parameters.signature_method == "limma":
            limma = rpy2LIMMA()
            signature_frame = limma.calculate(self.data_frame, self.positive_samples, self.negative_samples, self.output_file)
            signature_frame.columns = [self.signature_name]
            signature_frame.to_csv(self.output_file, sep = "\t", index_label = "id")

class mergeAllSignatures(Target):
    def __init__(self, parameters, directory):
        Target.__init__(self, time=10000)
        self.parameters = parameters
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## identify files
        all_files = os.listdir("analysis/signature_files")
        bootstrap_files = filter(lambda x: x.startswith("_b"), all_files)
        null_files = filter(lambda x: x.startswith("_n"), all_files)
        real_files = list(set(all_files) - set(bootstrap_files) - set(null_files))

        ## merge data
        real_signature_frame = pandas.DataFrame()
        for file in real_files:
            join_frame = pandas.read_csv("analysis/signature_files/%s" % (file), index_col = 0 , sep = "\t")
            real_signature_frame = real_signature_frame.join(join_frame, how = "outer")
        if len(real_files) > 0:
            real_signature_frame.to_csv(self.parameters.signature_file, sep = "\t", index_label = "id")
        bootstrap_signature_frame = pandas.DataFrame()
        for file in bootstrap_files:
            join_frame = pandas.read_csv("analysis/signature_files/%s" % (file), index_col = 0 , sep = "\t")
            bootstrap_signature_frame = bootstrap_signature_frame.join(join_frame, how = "outer")
        if len(bootstrap_files) > 0:
            bootstrap_signature_frame.to_csv("bootstrap_%s" % (self.parameters.signature_file), sep = "\t", index_label = "id")
        null_signature_frame = pandas.DataFrame()
        for file in null_files:
            join_frame = pandas.read_csv("analysis/signature_files/%s" % (file), index_col = 0 , sep = "\t")
            null_signature_frame = null_signature_frame.join(join_frame, how = "outer")
        if len(null_files) > 0:
            null_signature_frame.to_csv("null_%s" % (self.parameters.signature_file), sep = "\t", index_label = "id")
        
        ## remove signature_files
        # shutil.rmtree("analysis/signature_files")
        
def main():
    ## check for fresh run
    if os.path.exists(".jobTree"):
        logging.warning("WARNING: '.jobTree' directory found, remove it first to start a fresh run\n")
    
    ## parse arguments
    parser = OptionParser(usage = "%prog [options] data_matrix phenotype_matrix")
    Stack.addJobTreeOptions(parser)
    parser.add_option("--jobFile",
                      help = "Add as child of jobFile rather than new jobTree")
    parser.add_option("-b", "--bootstrap", dest = "bootstrap_size", default = 0,
                      help = "")
    parser.add_option("-n", "--null", dest = "null_size", default = 0,
                      help = "")
    parser.add_option("-p", "--permute", dest = "null_permute", default = "paradigm",
                      help = "")
    parser.add_option("-m", "--method", dest = "signature_method", default = "sam",
                      help = "")
    parser.add_option("-z", "--seed", dest = "seed", default = None,
                      help = "")
    options, args = parser.parse_args()
    logging.info("options: %s" % (str(options)))
    print "Using Batch System '%s'" % (options.batchSystem)
    
    if len(args) != 2:
        logging.error("ERROR: incorrect number of arguments\n")
        sys.exit(1)
    data_file = os.path.abspath(args[0])
    phenotype_file = os.path.abspath(args[1])
    
    ## set parameters
    parameters = Parameters(signature_method = options.signature_method,
                            bootstrap_size = int(options.bootstrap_size),
                            null_size = int(options.null_size),
                            null_method = options.null_permute,
                            random_seed = options.seed)
    
    ## run
    s = Stack(branchSignatures(data_file,
                               phenotype_file,
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
    from signature import *
    main()
