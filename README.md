PathMark: Identification of Pathway Markers of Interest
========

Current Version 
--------

2.0

Authors
--------

Sam Ng and Joshua M. Stuart


Requirements
--------

- [python](http://www.python.org/) >= 2.7
   - [scipy](http://www.scipy.org/) >= 0.12.0
   - [numpy](http://numpy.scipy.org/)
   - [pandas](http://pandas.pydata.org/)
   - [rpy2](http://rpy.sourceforge.net/rpy2/doc-2.4/html/overview.html)
   - [networkx](http://networkx.github.io/)
- [R](http://cran.r-project.org/)
   - [siggenes](http://www.bioconductor.org/packages/release/bioc/html/siggenes.html)
   - [edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)
   - [limma](http://www.bioconductor.org/packages/release/bioc/html/limma.html)
   - [samr](http://cran.r-project.org/web/packages/samr/index.html)

Installation
-------

- Install dependencies
- Download the pathmark-scripts repository to the desired location
- Run "make" in pathmark-scripts/ to generate source files; may also require downloading the paradigm-scripts repository (see https://github.com/ucscCancer/paradigm-scripts)
- Source init files for pathmark-scripts and paradigm-scripts (init.sh for bash and init.csh for csh)
- Run code on example data in examples/ with "make"

Command-Line
------
```
galaxyPATHMARK.py [options] data_matrix phenotype_matrix pathway_file

data_matrix - sample by feature data file
phenotype_matrix - phenotype by sample dichotomy file
pathway_file - PARADIGM pathway interactions file

-b bootstrap_size - number of bootstrap sample to be generated for robustness estimation
-n null_size - number of null samples to be generated for significance estimation
-p null_permute - permutation method for null model
-m signature_method - differential method for computing signatures
-f filter_parameters - filter threshold coefficients
-t diffusion_time - heat diffusion time for signature scores across the network
-u - apply hub filter that includes hubs with high representation of its children
-z random_seed - fix random seed
```

Folders
------
- bin : executables
- cytoscape : CytoscapeWeb and Cytoscape visualization files
- examples : BCCL inputs for demonstration purposes

Contact
------
Feature requests, comments and requests for clarification should all be sent to the author at <sng@soe.ucsc.edu>. 
I will try to respond quickly to all requests, so feel free to email me!
