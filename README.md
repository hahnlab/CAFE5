# CAFE

Software for **C**omputational **A**nalysis of gene **F**amily **E**volution

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


# History

CAFE v4.0 was the first release in a regular series of releases in order to make
CAFE easier and more user-friendly, in addition to adding features and fixing bugs.

CAFE v3.0 was a major update to CAFE v2.1. Major updates in 3.0 included: 1) the ability to correct 
for genome assembly and annotation error when analyzing gene family evolution 
using the **errormodel** command. 2) The ability to estimate separate birth (λ) and 
death (μ) rates using the **lambdamu** command. 3) The ability to estimate error in 
an input data set with iterative use of the errormodel command using the 
accompanying python script **caferror.py**. This version also included the addition of the 
**rootdist** command to give the user more control over simulations.

