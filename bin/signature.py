#!/usr/bin/env python
"""
signature.py

Author: Sam Ng
Last Updated: 2014-06-23
"""
import math, os, random, re, shutil, sys, types
from copy import deepcopy

import pandas
import pandas.rpy.common as common
from rpy2 import robjects
from rpy2.robjects.packages import importr
from scipy import stats

from optparse import OptionParser
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

## default variables
null_prefixes = ['na_', 'nw_']          ## valid null samples must start with one of these
                                        ## prefixes and end in a sample name

#### NOTE BLOCK
#### - Test LIMMA and add in data counts check
#### - Add in permute and paradigm null_method

## pm classes
class Parameters:
    """
    Stores parameters used for this [signature.py specific]
    """
    def __init__(self, random_seed = 0, bootstrap_size = 10, bootstrap_proportion = 0.85, bootstrap_replacement = False, null_method = 'labels', null_size = 0, signature_method = 'sam'):
        self.random_seed = random_seed
        self.bootstrap_size = bootstrap_size
        self.bootstrap_proportion = bootstrap_proportion
        self.bootstrap_replacement = bootstrap_replacement
        self.null_method = null_method
        self.null_size = null_size
        self.signature_method = signature_method
        self.signature_file = 'signature.tab'

class rpy2SignatureGenes:
    def __init__(self):
        self.siggenes = importr('siggenes')
    def calculate(self, method, data_frame, positive_samples, negative_samples):
        ## construct matrix_r
        r = robjects.r
        samples = []
        for sample in data_frame.axes[1]:
            if sample in positive_samples + negative_samples:
                samples.append(sample)
        features = data_frame.axes[0]
        matrix = data_frame[samples]
        matrix_r = common.convert_to_r_matrix(matrix)
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
            r.c(r.attributes(sam_out).rx2('d')),
            r.c(r.attributes(sam_out).rx2('vec.false')),
            r.c(r.attributes(sam_out).rx2('q.value')),
            r.c(r.attributes(sam_out).rx2('p.value')),
            r.c(r.attributes(sam_out).rx2('s'))
        )
        ## return results as a data_frame
        ocols = ['Score', 'FalseCalls', 'Q-value', 'P-value', 'StdDev']
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
        return(pandas.DataFrame(pandas.Series(t_map, name = 'T-statistic')))

class rpy2LIMMA:
    def __init__(self):
        self.siggenes = importr('siggenes')
    def calculate(self, method, data_frame, positive_samples, negative_samples):
        ## construct matrix_r
        r = robjects.r
        samples = []
        for sample in data_frame.axes[1]:
            if sample in positive_samples + negative_samples:
                samples.append(sample)
        features = data_frame.axes[0]
        matrix = data_frame[samples]
        matrix_r = common.convert_to_r_matrix(matrix)
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
            r.c(r.attributes(sam_out).rx2('d')),
            r.c(r.attributes(sam_out).rx2('vec.false')),
            r.c(r.attributes(sam_out).rx2('q.value')),
            r.c(r.attributes(sam_out).rx2('p.value')),
            r.c(r.attributes(sam_out).rx2('s'))
        )
        ## return results as a data_frame
        ocols = ['Score', 'FalseCalls', 'Q-value', 'P-value', 'StdDev']
        output = {}
        for j, col in enumerate(ocols):
            row = {}
            for i, n in enumerate(features):
                # print n, col, sam_att.rx(i + 1, j + 1)[0]
                row[n] = sam_att.rx(i + 1, j + 1)[0]
            # print row
            output[col] = row
        return(pandas.DataFrame(output))

## pm functions
def logger(message, file = None, die = False):
    """
    Writes messages to standard error [2014-3-1]
    """
    if file is None:
        sys.stderr.write(message)
    else:
        o = open(file, 'a')
        o.write(message)
        o.close()
    if die:
        sys.exit(1)

def generateDichotomies(phenotype, phenotype_frame, data_samples, reverse = True):
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
        phenotype_categories.sort(reverse = reverse)
        for index in range(len(phenotype_categories)):
            if index == 1 and len(phenotype_categories) == 2:
                break
            signature_name = '%s=%s_SIGNATURE' % (phenotype, phenotype_categories[index])
            positive_samples = list(set(phenotype_map[phenotype_map == phenotype_categories[index]].index) & set(data_samples))
            negative_samples = list(set(phenotype_map[phenotype_map != phenotype_categories[index]].index) & set(data_samples))
            yield(signature_name, positive_samples, negative_samples)

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
        if not os.path.exists('analysis'):
            os.mkdir('analysis')
        assert(not os.path.exists('analysis/signature_files'))
        os.mkdir('analysis/signature_files')
        
        ## read in data and phenotype matrices
        data_frame = pandas.read_csv(self.data_file, sep = '\t', index_col = 0)
        data_samples = list(data_frame.columns)
        phenotype_frame = pandas.read_csv(self.phenotype_file, sep = '\t', index_col = 0)
        phenotypes = list(phenotype_frame.columns)
        for phenotype in phenotypes:
            assert(not phenotype.startswith('_'))
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
        random.seed(self.parameters.random_seed+12321)
        
        ## queue real
        self.addChildTarget(computeSignatures('analysis/signature_files/%s' % (self.signature_name),
                                              self.data_frame,
                                              self.signature_name,
                                              self.positive_samples,
                                              self.negative_samples,
                                              self.parameters,
                                              self.directory))
        
        ## generate and queue bootstraps
        bootstrap_positive_count = int(round(self.parameters.bootstrap_proportion*len(self.positive_samples)))
        bootstrap_negative_count = int(round(self.parameters.bootstrap_proportion*len(self.negative_samples)))
        for bootstrap_index in range(self.parameters.bootstrap_size):
            bootstrap_positive_data_frame = self.data_frame[random.sample(self.positive_samples, bootstrap_positive_count)].copy()
            bootstrap_negative_data_frame = self.data_frame[random.sample(self.negative_samples, bootstrap_negative_count)].copy()
            bootstrap_data_frame = bootstrap_positive_data_frame.join(bootstrap_negative_data_frame, how = 'outer')
            self.addChildTarget(computeSignatures('analysis/signature_files/_b%s_%s' % (bootstrap_index + 1, self.signature_name),
                                                  bootstrap_data_frame,
                                                  '%s:%s' % (self.signature_name, bootstrap_index + 1),
                                                  self.positive_samples,
                                                  self.negative_samples,
                                                  self.parameters,
                                                  self.directory))
        
        ## generate and queue nulls
        for null_index in range(self.parameters.null_size):
            if self.parameters.null_method == 'labels':
                null_data_frame = self.data_frame.copy()
                null_data_frame.columns = random.sample(self.data_frame.columns, len(self.data_frame.columns))
            self.addChildTarget(computeSignatures('analysis/signature_files/_n%s_%s' % (null_index + 1, self.signature_name),
                                                  null_data_frame,
                                                  '%s:%s' % (self.signature_name, null_index + 1),
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
        
        if self.parameters.signature_method == 'sam':
            siggenes = rpy2SignatureGenes()
            signature_frame = siggenes.calculate(self.parameters.signature_method, self.data_frame, self.positive_samples, self.negative_samples)[['Score']]
            signature_frame.columns = [self.signature_name]
            signature_frame.to_csv(self.output_file, sep = '\t', index_label = 'id')
        elif self.parameters.signature_method == 'ttest':
            siggenes = scipySignatureGenes()
            signature_frame = siggenes.calculate(self.parameters.signature_method, self.data_frame, self.positive_samples, self.negative_samples)
            signature_frame.columns = [self.signature_name]
            signature_frame.to_csv(self.output_file, sep = '\t', index_label = 'id')
        elif self.parameters.signature_method == 'limma':
            siggenes = rpy2SignatureGenes()
            signature_frame = siggenes.calculate(self.parameters.signature_method, self.data_frame, self.positive_samples, self.negative_samples)[['Score']]
            signature_frame.columns = [self.signature_name]
            signature_frame.to_csv(self.output_file, sep = '\t', index_label = 'id')

class mergeAllSignatures(Target):
    def __init__(self, parameters, directory):
        Target.__init__(self, time=10000)
        self.parameters = parameters
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## identify files
        all_files = os.listdir('analysis/signature_files')
        bootstrap_files = filter(lambda x: x.startswith('_b'), all_files)
        null_files = filter(lambda x: x.startswith('_n'), all_files)
        real_files = list(set(all_files) - set(bootstrap_files) - set(null_files))

        ## merge data
        real_signature_frame = pandas.DataFrame()
        for file in real_files:
            join_frame = pandas.read_csv('analysis/signature_files/%s' % (file), index_col = 0 , sep = '\t')
            real_signature_frame = real_signature_frame.join(join_frame, how = 'outer')
        if len(real_files) > 0:
            real_signature_frame.to_csv(self.parameters.signature_file, sep = '\t', index_label = 'id')
        bootstrap_signature_frame = pandas.DataFrame()
        for file in bootstrap_files:
            join_frame = pandas.read_csv('analysis/signature_files/%s' % (file), index_col = 0 , sep = '\t')
            bootstrap_signature_frame = bootstrap_signature_frame.join(join_frame, how = 'outer')
        if len(bootstrap_files) > 0:
            bootstrap_signature_frame.to_csv('bootstrap_%s' % (self.parameters.signature_file), sep = '\t', index_label = 'id')
        null_signature_frame = pandas.DataFrame()
        for file in null_files:
            join_frame = pandas.read_csv('analysis/signature_files/%s' % (file), index_col = 0 , sep = '\t')
            null_signature_frame = null_signature_frame.join(join_frame, how = 'outer')
        if len(null_files) > 0:
            null_signature_frame.to_csv('null_%s' % (self.parameters.signature_file), sep = '\t', index_label = 'id')
        
        ## remove signature_files
        shutil.rmtree('analysis/signature_files')
        
def main():
    ## check for fresh run
    if os.path.exists('.jobTree'):
        logger('WARNING: .jobTree directory found, remove it first to start a fresh run\n')
    
    ## parse arguments
    parser = OptionParser(usage = '%prog [options]')
    Stack.addJobTreeOptions(parser)
    parser.add_option('--jobFile', help='Add as child of jobFile rather than new jobTree')
    parser.add_option('-d', '--data', dest='data_file', default=None)
    parser.add_option('-p', '--phenotype', dest='phenotype_file', default=None)
    parser.add_option('-m', '--method', dest='signature_method', default='sam')
    parser.add_option('-n', '--null', dest='null_size', default=0)
    parser.add_option('-b', '--bootstrap', dest='bootstrap_size', default=0)
    options, args = parser.parse_args()
    logger('Using Batch System : %s\n' % (options.batchSystem))
    
    if len(args) == 1:
        if args[0] == 'clean':
            command = 'rm -rf .jobTree analysis'
            logger(command)
            os.system(command)
            sys.exit(0)
    
    assert(len(args) == 0)
    
    ## set parameters
    parameters = Parameters(signature_method = options.signature_method,
                            null_size = options.null_size,
                            bootstrap_size = options.bootstrap_size)
    
    ## run
    s = Stack(branchSignatures(options.data_file,
                               options.phenotype_file,
                               parameters,
                               os.getcwd().rstrip('/')))
    if options.jobFile:
        s.addToJobFile(options.jobFile)
    else:
        if options.jobTree == None:
            options.jobTree = "./.jobTree"
        
        failed = s.startJobTree(options)
        if failed:
            logger('%d jobs failed\n' % failed)
        else:
            os.system('rm -rf .lastjobTree')
            os.system('mv .jobTree .lastjobTree')

if __name__ == '__main__':
    from signature import *
    main()
