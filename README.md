# NGS-QC-Panel
This is the code repository for the paper "Authentication, characterization and contamination detection of cell lines, xenografts and organoids by barcode deep NGS sequencing" published at NAR Genomics and Bioinformatics. Interested readers can consult the paper for algorithm details. 

## Requirements
The programs were developed and tested in perl(v5.28.1) and R(v3.6), they should run in other versions as well under both Linux and Windows. R can be downloaded from 
<https://cran.r-project.org>. Perl can be downloaded from <https://www.perl.org/get.html>. 

In addition, the following package is required for perl
* Cwd

The following packages are required for R 
* mclust
* MASS
* EnvStats


## Contents
### Programs
There are two perl programs, one R script and an accessory file

**MajorComponentEstimation.pl**:  perl script for inferring the major component of a sample


**ContaminationRatioEstimationS1.pl**: perl script for estimating the heterogeneity/contamination ratio of a sample--step 1


**ContaminationRatioEstimationS2.R**: R script for estimating the heterogeneity/contamination ratio of a sample--step 2


**binomial_threshold.txt**: this file is used by ContaminationRatioEstimationS1.pl 


### Datasets
There are three datasets used for illustrating the inference/estimation.

**mutfre**:  SNP fingerprints for 1055 cell lines

**testedcelllines_Table2**: .depth files for 78 cell line samples used in Table 2 of the paper

**testedcelllines_Table3**: .depth files for 22 cell line samples used in Table 3 of the paper

## Running the programs
1. Inference of the major component of a sample

```perl 
perl MajorComponentEstimation.pl 
```

The command line output can be saved and edited as Sample_MajorComponent.txt in the folder testedcelllines_Table2 and testedcelllines_Table3.

2. Coarse estimation of the heterogeneity/contamination ratio of a sample

```perl 
perl ContaminationRatioEstimationS1.pl 
```

There is one .SNPratio output file for each .depth file, saved in the folder testedcelllines_Table2 and testedcelllines_Table3. There is also a SNPratio.pdf file giving the distribution of heterogeneity/contamination ratios.

3. Precise estimate of the heterogeneity/contamination ratio of a sample
```R
Rscript ContaminationRatioEstimationS2.R
```

This program should be run with pwd set to testedcelllines_Table2 or testedcelllines_Table3, it generates a file called ContaminationRatio.csv.

