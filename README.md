# CAFE

<div>
<h3>
Software for <bold>C</bold>omputational <bold>A</bold>nalysis of gene <bold>F</bold>amily <bold>E</bold>volution
</h3>
</div>

The purpose of CAFE is to analyze changes in gene family size in a way that 
accounts for phylogenetic history and provides a statistical foundation for 
evolutionary inferences. The program uses a birth and death process to model gene 
gain and loss across a user-specified phylogenetic tree. The distribution of family 
sizes generated under this model can provide a basis for assessing the significance 
of the observed family size differences among taxa.

This repository contains code for _CAFE 5_, an updated version of the CAFE code base.
CAFE 5 showcases new functionalities while keeping or improving several of the features 
available in prior releases.

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

[![Build Status](https://travis-ci.org/hahnlab/CAFExp.svg?branch=master)](https://travis-ci.org/hahnlab/CAFExp)

History
-------

The original development of the statistical framework and algorithms are 
described by Hahn, _et al._ (2005) and later implemented in the software package
CAFE by De Bie, _et al._ (2006).

CAFE v2.0 (Hahn, Demuth, _et al._ 2007; Hahn, Han, _et al._ 2007) Included software 
updates and functionality that allowed users to specify different λ values for different 
branches on their input tree.  

CAFE v3.0 (Han, _et al._ 2013) was a major update to CAFE 2.0, with added functionality 
that included: 1) the ability to correct for genome assembly and annotation
error.  2) The ability to estimate separate birth (λ) and death (μ)
rates using the lambdamu command. 3) The ability to estimate error in an
input data set with iterative use of the errormodel command using the
accompanying python script caferror.py. This version also included the
addition of the rootdist command to give the user more control over
simulations.

CAFE v4.0 was primarily a maintenance update. At that time formal issue
tracking and a user forum was introduced.

CAFE v5.0. (Current Release) Another major update, CAFE5 showcases new 
functionalities while keeping or improving several of the features available 
in prior releases.

How to Cite
-----------
___

If you use CAFE5 in your work, please cite the application as

> (xxx - a Zenodo DOI to be determined when released)

Original development of the statistical framework and algorithms implemented in CAFE are 
published in:

- Hahn, M. W., T. De Bie, J. E. Stajich, C. Nguyen, and N. Cristianini. 2005. Estimating the tempo and mode of gene family evolution from comparative genomic data. _Genome Research_ 15:1153–1160.

- De Bie, T., N. Cristianini, J. P. Demuth, and M. W. Hahn. 2006. CAFE: a computational tool for the study of gene family evolution. _Bioinformatics_ 22:1269–1271.

The citation for CAFE v2.0 is:

- Hahn, M. W., J. P. Demuth, and S.-G. Han. 2007. Accelerated rate of gene gain and loss in primates. _Genetics_ 177:1941–1949. Genetics.

The citation for CAFE v3.1 and v4.0 is:

- Han, M. V., G. W. C. Thomas, J. Lugo-Martinez, and M. W. Hahn. 2013. Estimating Gene Gain and Loss Rates in the Presence of Error in Genome Assembly and Annotation Using CAFE 3. _Mol. Biol. Evol._ 30:1987–1997.

New Functionality
-----------------
______
-   Among Family Rate Variation (AFRV) using a discrete gamma model with 
    a jointly optimized alpha shape parameter.  The birth-death model estimates 
    the posterior probabilities of each gene family belonging to different 
    evolutionary rate categories. The rates of each of the K categories are 
    determined in similar fashion to what is done in nucleotide sequence analyses:
    through a discrete gamma distribution that has its alpha parameter
    estimated by maximum likelihood.

-   Ancestral state reconstruction is now done jointly (all ancestral
    states are inferred simultaneously) with Pupko’s algorithm (Pupko et
    al. 2000).

-   The error model is now optimized numerically by CAFE5 directly,
    rather than the Python script required in previous versions.

-   Outputs are now internally parsed by CAFE5 (no external
    dependencies or scripts necessary) into summary tables.

What CAFE5 does
--------
______

CAFE5 implements a birth-death model for evolutionary inferences
about gene family evolution. Its main task is the maximum-likelihood
estimation of a global or local gene family evolutionary rates (lambda
parameter) for a given data set. Briefly, CAFE5 can:

-   Compare scenarios in which the whole phylogeny shares the same
    (global) lambda vs. scenarios in which different parts of the
    phylogeny share different (local) lambdas.

-   Classify specific gene families as "fastly evolving” in at least two
    different ways (see documentation below).

-   Infer gene family counts at all internal nodes (ancestral
    populations) of the phylogeny provided as input. By comparing
    different nodes in the phylogeny, the user will be able to tell
    along which branches gene families have contracted or expanded.

-   Account for non-biological factors (e.g., genome sequencing and
    coverage differences, gene family clustering errors, etc.) leading
    to incorrect gene family counts in input files. This is done with an
    error model.

What CAFE5 does NOT do
-----------

-   Estimate a phylogeny from gene families or gene sequence
    alignments. CAFE5 also does not convert a phylogeny with
    branches in expected substitutions per site into a time tree (an
    ultrametric tree with branch lengths in time units). This task
    should be conducted by the user prior to CAFE5 analyses.

-   Implement clustering algorithms that identify (or verify the
    legitimacy of) gene families. It is entirely up to the user to pick
    a gene family identification method and carry out this task prior to
    CAFE5 analyses.

-   Predict gene family function or infer enrichment of functional
    classes.


Installation
============

Download
--------

The Github page for CAFE5 is https://github.com/hahnlab/CAFExp 

Navigate to a directory that you typically keep source code in and do one of the following:

To copy the master directory as a .zip file:

    $ wget https://github.com/hahnlab/CAFExp/archive/master.zip


Compile
-------

Move into the CAFE5 folder, and type

    $ autoconf
    $ ./configure
    $ make

The CAFE5 executable will be put inside the “bin” folder. Then
just copy the binary file to somewhere on your path, such as
/usr/local/bin. Alternatively, add /path/to/CAFE5/bin/ to your
\$PATH variable (in your .bashrc or .bash profile).


### OSX users

 If you encounter an error during the build that looks like: 
<pre>
    src/matrix_cache.cpp:2:10: <b>fatal error:</b> 'omp.h' file not found
</pre>

You may need to install gcc, find the omp.h file and add it to your path. If you already have gcc installed you may have to find this file and soft link it to a directory that the compiler can find.  Try finding the missing file by:

    $ find / -name omp.h

This found the omp.h file in the [Homebrew](https://brew.sh/) installation of gcc.  I can now soft link (or copy) it to a directory where it can be found.

    ln -sv /usr/local/Cellar/gcc/7.3.0/lib/gcc/7/gcc/x86_64-apple-darwin17.3.0/7.3.0/include/omp.h   /usr/local/include/

Now run 
    
    $ make
 


# Running CAFE5 
----

## Quick Start 


For a typical CAFE analysis, users are most interested in determining two things:
1) Which gene families are rapidly evolving 
2) The branches of the tree on which these families are rapidly evolving

This type of analysis requires a minimum of two input files:
1) A tab-delimited family "counts" file that contains a column for a description of the gene family,
       the unique ID for each family, and a column for each taxon that has count data for each family.
       This file is acquired by first peforming a clustering analysis, often using software such as 
       OrthoMCL, SwiftOrtho, FastOrtho, OrthAgogue, or OrthoFinder and then parsing the output into a table
       like the one below (Note: if a functional description is not desired, include this column anyway with a place holder as below (null)).

Example: mammal_gene_families.txt
```
Desc	Family ID	human	chimp	orang	baboon	gibbon	macaque	marmoset rat	mouse	cat	horse	cow
ATPase	ORTHOMCL1	 52	 55	 54	 57	 54	  56	  56	 53	 52	57	55	 54
(null)	ORTHOMCL2	 76	 51	 41	 39	 45	  36	  37	 67	 79	37	41	 49
HMG box	ORTHOMCL3	 50	 49	 48	 48	 46	  49	  48	 55	 52	51	47	 55
(null)	ORTHOMCL4	 43	 43	 47	 53	 44	  47	  46	 59	 58	51	50	 55
Dynamin	ORTHOMCL5	 43	 40	 43	 44	 31	  46	  33	 79	 70	43	49	 50
......
......
DnaJ	ORTHOMCL10016	 45	 46	 50	 46	 46 	  47	  46	 48	 49	45	44	 48
``` 
2) The other required input file should contain a binary, rooted, ultrametric, tree in Newick format.  Typically
one obtains this tree using one of several molecular dating methods. If you are unsure if your tree is binary,
rooted, or ultrametric CAFE will report this when you try to use it for an analysis. Alternatively, you can use the R package,
Ape with its included functions: is.ultrametric, is.rooted, and is.binary.  

Example: mammals_tree.txt
```
((((cat:68.710507,horse:68.710507):4.566782,cow:73.277289):20.722711,(((((chimp:4.444172,human:4.444172):6.682678,orang:11.126850):2.285855,gibbon:13.412706):7.211527,(macaque:4.567240,baboon:4.567240):16.056992):16.060702,marmoset:36.684935):57.315065):38.738021,(rat:36.302445,mouse:36.302445):96.435575);
```
To get a list of commands just call CAFE with no arguments:

    $ cafe

To estimate lambda with no among family rate variation issue the command:

    $ cafe -i mammal_gene_families.txt -t mammal_tree.txt

To incorporate among family rate variation with both lambda and alpha estimated and three discrete gamma rate categories:

    $ cafe -i mammal_gene_families.txt -t mammal_tree.txt -k 3


To estimate separate lambda values for different lineages in the tree, first identify the branches to which each lambda will apply.
This can be done by making a copy of your tree, and substituting the lambda identifier (1,2,3, etc.) for the branch length values.
For example, to apply a different lambda to the branches leading to human, chimp, and their ancestor, modify the branches as below.

Example chimphuman_separate_lambda.txt:
<pre>
((((cat:1,horse:1):1,cow:1):1,(((((<b>chimp:2,human:2):2</b>,orang:1):1,gibbon:1):1,(macaque:1,baboon:1):1):1,marmoset:1):1):1,(rat:1,mouse:1):1);
</pre>
For this tree, lambda #2 will be applied to branches leading to human, chimp, and their ancestor while lambda #1 will be applied to all other branches of the tree. 


To run this analysis with both lambdas estimated:


    $ cafe -i mammal_gene_families.txt -t mammal_tree.txt -y chimphuman_separate_lambda.txt 

***Caveats***

- **Always** perform multiple runs to ensure convergence, especially if multiple gamma rate categories or lambdas are used.
- More gamma rate categories (-k) does not always mean a better fit to the data. While -k=2 nearly always fits the data better than -k=1, it may be the case that -k=5 has a _worse_ likelihood than -k=3, and convergences between runs is more difficult with more categories. Try several and see what works.
- We recommend using the -o flag to assign a unique name to the output directory for each run so that results from previous runs are not overwritten.

- For all but the simplest of data sets, searching for multiple lambdas with multiple rate categories will result in a failure of convergence to a single optimum between runs. 

---


## Slow Start


CAFE5 performs three different operations on either one or two
models. The operations are

-   Estimate Lambda - the traditional function of CAFE. Takes a tree and
    a file of gene family counts, and performs a maximum likelihood
    calculation to estimate the most likely rate of change across the
    entire tree.

-   Simulate - Given specified values, generate an artificial list of
    gene families that matches the values. To generate a simulation,
    pass the --simulate or -s parameter. Either pass a count of families
    to be simulated with the parameter (--simulate=1000) or pass a
    --rootdist (-f) parameter with a file containing the distribution to
    match (see \[rootdist\] for the file format).

The models are
-   Base - Perform computations as if no gamma function is available

-   Gamma - Perform computations as if each gene family can belong to a
    different evolutionary rate category. To use Gamma modelling, pass
    the -k parameter specifying the number of categories to use.

Unlike earlier versions, CAFE5 does not require a script. All
options are given at once on the command line. Here is an example:

    cafexp -t examples/tree.txt -i genefamilies.txt -p

In this example, the -t parameter specifies a file containing the tree
that CAFE uses; and the -i parameter specifies a list of gene families.
The -p, in this instance given without a parameter, indicates that the
root equilibrium frequency will not be a uniform distribution.

Parameters
----------

-   **--fixed\_alpha, -a**

    Alpha value of the discrete gamma distribution to use in category
    calculations. If not specified, the alpha parameter will be
    estimated by maximum likelihood.

-   **--error_model, -e**

    Error model file path

-   **--rootdist, -f**

    Root distribution file path

-   **--infile, -i**

    character or gene family file path

-	**--n\_gamma\_cats, -k**

    Number of gamma categories to use. If specified, the Gamma model
    will be used to run calculations; otherwise the Base model will be
    used.

-   **--fixed\_lambda, -l**

    Single Lambda value

-   **--fixed\_multiple\_lambdas**

    Multiple lambda values, comma separated

-   **--output\_prefix, -o**

    Directory to place output files. Defaults to "results"

-   **--poisson, -p**

    Use a Poisson distribution for the root frequency distribution.
    Without specifying this, a normal distribution will be used. A value
    can be specified (-p10, or --poisson=10); otherwise the distribution
    will be estimated from the gene families.

-   **--chisquare\_compare, -r**

    Chi square compare

-   **--simulate, -s**

    Simulate families. Either provide an argument of the number of families
	to simulate (-s100, or --simulate=100) or provide a rootdist file giving a set
	of root family sizes to match. Without such a file, the families will be generated
	with root sizes selected randomly between 0 and 100.

-   **--tree, -t**

    Tree file path - Required for estimation

-   **--lambda\_tree, -y**

    Lambda tree file path

Input files
-----------

- 
- Tree files

    A tree file is specified in Newick format.

        (A:1,B:1);

    An example may be found in the examples/test\_tree.txt file.

    When nodes are mentioned in output files and logs, they are named by
    the leaf nodes descending from them, in alphabetical order. For
    example, the parent of A and B in this example will be named “AB”. A
    node containing descendant species of cow, whale, giraffe, and
    manatee would be named “cowgiraffemanateewhale”.

-   Family files

    Family files can be specified in the CAFE input format:

        Desc    Family ID       A       B
        (null)  1       1       2
        (null)  2       2       1
        (null)  3       3       6
        (null)  4       6       3

    The file is tab-separated with a header line giving the order of the
    species. Each line thereafter consists of a description, a family
    ID, and counts for each species in the tree.

    Alternatively, the family file can be specified with a set of lines
    beginning with hashtags containing the species order:

        #A
        #B
        1       2     1
        2       1     2
        3       6     3
        6       3     4

    In this case, the family ID will be in the final column.

-   Root distributions

    A root distribution file takes the format “family\_size
    \[whitespace\] family\_count”, e.g.

        1 1
        2 5
        3 10
        4 15
        5 42

    An example may be found in the
    “examples/poisson\_root\_dist\_1000.txt” file.

-   Error models

    An error model consists of modifications of probabilities of moving
    from one family size to another through the tree. The file is
    structured as a series of lines containing the family size, the
    probability of moving to less than that size, the probability of
    that size staying the same, and the probabilities of the size
    becoming larger. Two header lines must be included: the maximum
    family size to process, and the differential of the probabilities.

        maxcnt:90
        cntdiff -1 0 1
        0 0.00 0.95 0.05
        1 0.05 0.9 0.05
        2 0.05 0.9 0.05
        3 0.05 0.9 0.05
        4 0.05 0.9 0.05

Output
------

All output will be stored to the "results" directory, unless another directory is specified with the "-o"
parameter.


-   _model_\_asr.txt

    The file will be named Base\_asr.txt or Gamma\_asr.txt, based on
    which model is in play. It contains the reconstructed states of the
    families, in the Nexus file format
    (https://en.wikipedia.org/wiki/Nexus\_file). A tree is provided for
    each family,with the expected family size set off with an underscore
    from the node ID. 

    In the case of the Gamma reconstruction, the Lambda multipliers for
    each category are given their own section in this file. In this case,
	ony the fastest families are printed.

	
-   _model_\_family\_results.txt

    The file will be named Base\_family\_results.txt or
    Gamma\_family\_results.txt, based on which model is in play. It
    consists of a header line giving the name of each node in the tree,
    followed by a line consisting of the family ID, an estimate of
    whether the change is significant(’y’ or ’n’) followed by a series
    of ’c’ (constant), ’i’ (increasing), or ’d’ (decreasing) showing
    whether the node is larger, smaller, or consistent with its parent
    node. The characters are separated by tabs.

        #FamilyID   pvalue  *   orang   gibbon  chimp   human     
        0           0.436   n   d       d       c       c 
        1           0.209   n   c       i       c       c
        2           0.002   y   c       c       i       i

    In the Gamma model, an additional set of probabilities are appended,
    representing the likelihood of the family belonging to each gamma
    category.

-   _model_\_clade\_results.txt

    The file will be named Base\_clade\_results.txt or
    Gamma\_clade\_results.txt, based on which model is in play. It
    consists of a header line, “Taxon Increase Decrease”, and for each
    node in the tree, a tab-separated count of the number of families
    which have increased and decreased for that node.

-   _model_\_branch\_probabilities.txt

	Contains a tab-separated list of the probabilities calculated for each clade
	and significant family. Probabilities are displayed as "N/A" if the parent
	and child have the same value. In the case of the Gamma model, only
	contains significant families that are calculated to be rapidly changing.
	
-   _model_\_family\_likelihoods.txt

    Using the Base model, a tab-separated file consisting of the header
    line “\#FamilyID Likelihood of Family”, and additional tab-separated
    lines consisting of the family ID and the posterior probability of
    that family.

    Using the Gamma model, the file contains the following values:

    -   [\#FamilyID]{}

    -   [Gamma Cat Median]{}

    -   [Likelihood of Category]{}

    -   [Likelihood of Family]{}

    -   [Posterior Probability]{}

    -   [Significant]{} The values for each family are listed on each
        tab-separated line.

-   _model_\_results.txt

    A file giving the name of the model that was selected (“Base” or
    “Gamma”), the final likelihood of that model, the final value of the
    rate of change of families (Lambda) that was calculated, and, if an
    error model was specified, the final value of that value (Epsilon).

-	[_model_.txt.change] - A tab-separated file listing, for each family 
        and clade, the difference between it and its parent clade in the 
		reconstruction that was performed.
		
-	[_model_.txt.count] - A tab-separated file listing, for each family 
        and clade, the reconstructed value in that clade.
		
-   simulation\_.txt

    In the case of simulation, a family file is
    generated with simulated data based on the given input parameters.

-   simulation\_truth\_.txt

    When simulating, an additional file is generated with simulated data
    for internal nodes included. The format is otherwise identical to
    the family file. Rather than node names, each internal node is
    assigned an integer ID.

Examples
========

Lambda Search
-------------

Search for a single lambda value using the mammal phylogeny and gene family set in the Examples directory:

    ../bin/cafexp -t mammals_tree.txt -i mammal_gene_families.txt -p -o singlelambda

In earlier versions, the following script would have returned the same
values:

    tree ((((cat:68.710507,horse:68.710507):4.566782,cow:73.277289):20.722711,(((((chimp:4.444172,human:4.444172):6.682678,orang:11.126850):2.285855,gibbon:13.412706):7.211527,(macaque:4.567240,baboon:4.567240):16.056992):16.060702,marmoset:36.684935):57.315065):38.738021,(rat:36.302445,mouse:36.302445):96.435575)
    load -i filtered_cafe_input.txt -t 10 -l reports/run6_caferror_files/cafe_log.txt
    lambda -s -t ((((1,1)1,1)1,(((((1,1)1,1)1,1)1,(1,1)1)1,1)1)1,(1,1)1)

Lambda Search with Multiple Lambdas
-----------------------------------

Search for separate lambda values for the chimp/human clade and the rest of the tree separately, using the mammal phylogeny and gene family set in the Examples directory:

    ../bin/cafexp -t mammals_tree.txt -i mammal_gene_families.txt -p -y chimphuman_separate_lambda.txt -o doublelambda

In earlier versions, the following script would have returned the same
values:

    ((((cat:68,horse:68):4,cow:73):20,(((((chimp:4,human:4):6,orang:11):2,gibbon:13):7,(macaque:4,baboon:4):16):16,marmoset:36):57):38,(rat:36,mouse:36):96) 
    load -i integral_test_families.txt -t 10
    lambda -s -t ((((1,1)1,1)1,(((((2,2)2,2)2,2)2,(1,1)1)1,1)1)1,(1,1)1)

Error Models
------------

Search for a single lambda value using the Newick tree of mammals in the
Samples folder, and the family files from the CAFE tutorial, applying an
error model:

    cafexp -t data/mammals_integral_tree.txt -i data/filtered_cafe_input.txt -p -l 0.01 -e data/cafe_errormodel_0.0548828125.txt

In earlier versions, the following script would have returned the same
values:

    tree ((((cat:68,horse:68):4,cow:73):20,(((((chimp:4,human:4):6,orang:11):2,gibbon:13):7,(macaque:4,baboon:4):16):16,marmoset:36):57):38,(rat:36,mouse:36):96) load -i integral_test_families.txt -t 10
    errormodel -all -model cafe_errormodel_0.0548828125.txt
    lambda -l 0.01 -t ((((1,1)1,1)1,(((((1,1)1,1)1,1)1,(1,1)1)1,1)1)1,(1,1)1) -score

Error Model Estimation
----------------------

Estimate an error model:

    cafexp -t mammals_integral_tree.txt -i filtered_cafe_input.txt -p -e errormodel.txt

In earlier versions, the following script would have returned the same
values:

    load -i filtered_cafe_input.txt -t 4 -l reports/log_run6.txt
    tree ((((cat:68.710507,horse:68.710507):4.566782,cow:73.277289):20.722711,(((((chimp:4.444172,human:4.444172):6.682678,orang:11.126850):2.285855,gibbon:13.412706):7.211527,(macaque:4.567240,baboon:4.567240):16.056992):16.060702,marmoset:36.684935):57.315065):38.738021,(rat:36.302445,mouse:36.302445):96.435575)
    lambda -s -t ((((1,1)1,1)1,(((((1,1)1,1)1,1)1,(1,1)1)1,1)1)1,(1,1)1)
    report reports/report_run6

Troubleshooting
===============

Known Limitations
-----------------

Because the random birth and death process assumes that each family has
at least one gene at the root of the tree, CAFE5 will not provide
accurate results if included gene families were not present in the most
recent common ancestor (MRCA) of all taxa in the tree. For example, even
if all taxa have a gene family size of 0, CAFE will assign the MRCA a
gene family of size 1, and include the family in estimation of the birth
and death rate. This difficulty does not affect analyses containing
families that go extinct subsequent to the root node.

If a change in gene family size is very large on a single branch, CAFE5 
may fail to provide accurate λ estimation and/or die during
computation. To see if this is a problem, look at the likelihood scores
computed during the λ search (reported in the log file if the
job finishes). If ALL scores are “-inf” then there is a problem with
large size changes and CAFE5 has calculated a probability of 0. Removing the
family with the largest difference in size among species and rerunning
CAFE5 should allow λ to be estimated on the remaining data.
If the problem persists, remove the family with the next largest
difference and proceed in a like manner until CAFE5 no longer finds
families with zero probability. However, if rapidly evolving families
are removed, care should be taken in interpretation of the estimated
average rate of evolution for the remaining data.

In very large phylogenetic trees there can be many independent lambda
parameters (2n - 2 in a rooted tree, where n is the number of taxa).
CAFE5 does not always converge to a single global maximum with large
numbers of λ parameters, and therefore can give misleading
results. To check for this you should always run the λ search
multiple times to ensure that the same estimated values are found. Also,
the likelihood of models with more parameters should always be lower
than models with fewer parameters, which may not be true if [CAFE5
]{}has failed to find a global maximum. If CAFE does not converge over
multiple runs, then one should reduce the number of parameters estimated
and try again.

Technical
=========

How does the optimizer work?
----------------------------

The Nelder-Mead optimization algorithm is used. It runs until it can
find a difference of less than 1e-6 in either the calculated score or
the calculated value, or for 10,000 iterations. The parameters that are
used for the optimizer are as follows:

-   rho: 1 (reflection)

-   chi: 2 (expansion)

-   psi: 0.5 (contraction)

-   sigma: 0.5 (shrink)

In some cases, the optimizer suggests values that cannot be calculated
(due to saturation, negative values, or other reasons) In this case, an
infinite score is returned and the optimizer continues.

When optimizing for an alpha value with a set number of clusters, if the
largest multiplier in the longest branch is saturated, the scorer will
return an infinite value. This will be noted at the end of the run with
text like:

    90 values were attempted (10% rejected)

showing that 10

    The following families had failure rates >20% of the time:
    Family6 had 22 failures
    Family9 had 19 failures

Certain options are available at compile-time for the optimizer. If
OPTIMIZER\_STRATEGY\_INITIAL\_VARIANTS is defined, the optimizer will
take several shorter attempts at various values before settling on one
value to continue on with. This may cause the optimizer to take more
iterations to finish but may have greater accuracy. If
OPTIMIZER\_STRATEGY\_PERTURB\_WHEN\_CLOSE is defined, the optimizer will
begin searching a wider range of values when it is getting close to a
solution. This attempts to get the optimizer out of a local optima it
may have found.

How does the simulator choose what lambda to use?
-------------------------------------------------

Although the user specifies the lambda, in order to give more family
variety a multiplier is selected every 50 simulated families. So
if 10,000 families are being simulated, 200 different lambdas will be
used.

When simulating without the gamma model, the multiplier is a random
value based on a normal distribution with a mean of 1 and a standard
deviation of 0.3.

When simulating WITH the gamma model, the multiplier is drawn directly
from a gamma distribution based on the selected alpha and a mean of 1.
If clustering is requested via the -k parameter, the selected cluster
multiplier is modified by a normal distribution with the mean at the
value of the multiplier, and a standard deviation intended to reflect
the number of clusters requested.

Initial Guesses
---------------

One of the most important concerns when searching a parameter space is
what initial values to choose. Each of the three values that CAFE can
search for has a particular initial guess strategy. For lambda values,
the formula is 1 / (longest tree branch times a random number between 0 and 1).
For epsilon values, the initial guesses are taken directly from the
provided error model. For gamma values, a random value taken from an
exponential distribution is used.

In some situations the values may fail. In this case, the scorer will
return an infinite value and the optimizer will retry initialization,
up to a number of attempts determined at compile time. If, after this number
of attempts, the optimizer continues to fail, it will abort. The most
likely cause of failure is too wide a variety of species sizes inside
certain families, and a message will be shown giving the most likely
families to remove from the analysis for success. 	

Acknowledgements
================

Many people have contributed to the CAFE project, either by code or
ideas. Thanks to:

-   Ben Fulton

-   Matthew Hahn

-   Mira Han

-	Fabio Mendes

-   Gregg Thomas

-	Dan Vanderpool


