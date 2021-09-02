# Mr.MtRobin: Multi-tissue transcriptome-wide Mendelian Randomization method ROBust to INvalid instrumental variables

The goal of `Mr.MtRobin` is to provide tools to conduct transcriptome-wide Mendelian randomization analysis
using the MR-MtRobin algorithm. MR-MtRobin is a Multi-tissue transcriptome-wide Mendelian Randomization method ROBust 
to INvalid instrumental variables that takes summary statistics from complex trait GWAS
and multi-tissue eQTL analyses as input and
uses a reverse regression random slope mixed model to infer whether a gene is
associated with a complex trait.

## Setup

To install and load functions from `Mr.MtRobin`, run the following:

  ```R
  devtools::install_github("kjgleason/Mr.MtRobin")
  library("Mr.MtRobin")
  ```

## Citation

To cite `Mr.MtRobin` in publications, please use:

Kevin J. Gleason, Fan Yang and Lin S. Chen. A robust two-sample Mendelian Randomization method integrating GWAS with multi-tissue eQTL summary statistics. Genet Epidemiol (2021); 45(4):353-371. doi:10.1002/gepi.22380.
