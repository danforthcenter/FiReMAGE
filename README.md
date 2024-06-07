# FiReMAGE: Filtering Results of Multispecies Analogous GWAS Experiments

FiReMAGE is a tool written in R that analyzes GWAS datasets of the same traits in multiple species to search for linked orthologous genes. 

R version 4.1.0 (2021-05-18)

	1. R Packages
	2. Directories & Accompanying files
	3. Function
	4. Options
	5. SNP input file format
	6. Customizing ortholog files
	7. Test data

## 1. R Packages
FiReMAGE will use and install the following packages:

	1. plyr and dplyr (Wickham 2011 & Wickham et al. 2018)
	2. tidyr (Wickham 2019)
	3. data.table (Dowle & Srinivasan 2018)
	4. doParallel (Microsoft & Weston 2017)
	5. reader (Cooper 2017)
	6. readr (Wickham & Hester 2020)
	7. iterators (Revolution Analytics & Weston 2017) 
	8. scales (Wickham 2017)
	9. ggplot2 (Wickham 2016)
	10. doRNG (Gaujoux 2023)
	11. stringr (Wickham 2022)
	11. docopt - optional for RStudio use, necessary for use in command line R (Jonge 2016)

## 2. Directories and Accompanying files

	1. Pipeline Rscript: FiReMAGE_[local|server]_version.R - examples for running FiReMAGE locally in RStudio or in remotely in command line. 
	2. functions Rscript: all functions used in the pipeline
	3. org_chromosome_coords.csv - located in the data subdirectory, contains all of the base pair lengths for each chromosome of each species in the comparison
	4. metaTable csv - template located in data/metaTable directory
	5. ortholog tables - data/OrthoFinder_orthologs/Orthogroups_5s.tsv is the file used in scripts. N0.tsv contains extra gene tree info from OrthoFinder
	6. gene files - in the data/phytozome directory, contains gene locations for all genes in each species
	7. GWAS SNPs - data/seedIonome_snps/ has the ionomic datasets for each species used to test FiReMAGE

## 3. Function
The metaTable file contains most of the information that FiReMAGE needs to filter GWAS datasets. Species, linkage ranges, and SNP/loci input information will all need to be updated for each new comparison. GWAS SNPs are read in from each species and collapsed into loci, to prevent redundant overlaps with genes.  GWAS/gene overlaps are identified in each species and then filtered to retain hits to orthologous genes. Only ortholog groups with hits in 3 or more species are kept.

![FiReMAGE pipeline flowchart](https://github.com/danforthcenter/FiReMAGE/blob/main/FiReMAGE_pipeline_flowchart.png?raw=true)

To assess whether the frequency of returned genes is due to noise in the genome, random permutations are used to measure the chance we can generate the same results by selecting SNPs throughout the genome at random. Using the leafcollapsedSNP.csv, the pipeline will generate the permutation datasets, randomly selecting the same number of loci from each chromosome and organism as the actual dataset. Because the actual dataset represents collapsed SNPs with variable ranges due to the amount of nearby SNPs, we replicate this in the permutation dataset by randomly assigning loci ranges from the real dataset to their permutation counterparts. The final permutation dataset with the SNPs and random ranges is saved in the output directory permutation_files/snps/. Note, if using the same SNP data, genome assemblies and linkage ranges, previously generated snp files can be used with the "-r" option. I find this useful for comparing results across different ortholog tables. Permutations are put through the same pipeline as the actual dataset, and false discovery rates and candidate gene likelihoods are calculated. Averages and quantiles from the actual and permutation datasets are saved in the results directory as AllSummaries.csv.

Note: 

	. If the SNP data is different from the most up-to-date versions, or other deviations like custom coordinates are to be used, they should be altered in the organism coordinates file. 

	. Make sure the format matches the test SNP files in the data/seedIonome_snps/ directory, and that the name includes the species name as spelled in the metaTable. This allows for multiple SNP files to be kept in one directory, and only the ones included in the metaTable will be read into the pipeline.
	
	. To work with the FiReMAGE pipeline, SNP files should include columns for the organism (org), trait, base pair position on chromosome, and chromosome (chr). They do not necessarily have to be in order, as the pipeline will order the file with the function snpFilePrep. Chromosome and base pair columns should be numerical (ie. 7, instead of Chr 7, Chr_7). The pipeline expects a different SNP file for each organism in the comparison; no need to concatenate all the files before had.

## 4. Options

	1. -h, --help : Shows a docopt help screen displaying all options
	2. -s <path>, --snps <path> : The path to the directory containing the SNP files for the comparison, required for the operation of FiReMAGE 
	3. -m <path>, --metaTable <path>: Path to read in the metadata table of organisms in the comparison, required for the operation of FiReMAGE. A template and sample metatable input is in ./data/metaTables/
	4. -f <path> --orthologFiles: The path to the directory containing ortholog files
	5. -o <path>, --output <path>: Files will be written to this path, the default is './FM_output/'
	6. -p <int>, --permutations <int>: Number of permutations for random dataset, default is 10. Recommendation is 1000 for full analysis, 100 for initial testing
	7. -r <path>, Path to previous loci permutations, useful for rerunning comparisons with different ortholog tables, default is output directory "-o" 
	8. -c <int>, --cores <int>: Number of cores to use for parallel loops, the default is 1 so it will naturally run on any system 

## 5. Customizing ortholog files
Users should downloaded the appropriate fasta protien files and run according to OrthoFinder v2.0 instructions (https://github.com/davidemms/OrthoFinder, Emms & Kelly, 2019). We obtained fasta protein files (and gene coordinates) from Phytozome v13 via their biomart interface (https://phytozome.jgi.doe.gov/biomart/martview).  The file we used for the ortholog table was in the N0.tsv file in the phylogenetic Hierarchical Orthogroups directory. If using the previous version of OrthoFinder, use the /Orthogroups/Orthogroups.tsv file (this has been deprecated in OrthoFinder v2.0). 

## 6. Test data

A shortened version of our ionomic data can be found in the example_test/ dir to quickly test your FiReMAGE set up before running a full job. An example script for running FiReMAGE in RStudio with the test data is also provided.The Orthogroup_hits.csv should be the same, and the permutation files may be slightly different due to the random distributions. If an exact comparison is desired, users can use the "-r" option and path to example_output. An example of the verbose output from the code can be seen in './test_output/test_log.txt'. This file is what a normal output (including all the warnings) looks like for a run of FiReMAGE. 

## References

Cooper, N. 2017. reader: Suite of Functions to Flexibly Read Data from Files. R package version 1.0.6. https://CRAN.R-project.org/package=reader

De Jonge, E. 2016. docopt: Command-Line Interface Specification Language. R package version 0.4.5.  https://CRAN.R-project.org/package=docopt

Dowle, M. & Srinivasan, A. 2018. data.table: Extension of `data.frame`. R package version 1.11.4. https://CRAN.R-project.org/package=data.table

Emms, D. M., & Kelly, S. 2015. OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome biology, 16(1), 1-14.

Gaujoux, R. 2023. doRNG: Generic Reproducible Parallel Backend for 'foreach' Loops. R package version 1.8.6. https://CRAN.R-project.org/package=doRNG 

Microsoft Corporation & Weston, S. 2017. doParallel: Foreach Parallel Adaptor for the 'parallel' Package. Rpackage version 1.0.11. https://CRAN.R-project.org/package=doParallel

Revolution Analytics and Weston, S. 2017. iterators: Provides Iterator Construct for R. R package version 1.0.9. https://CRAN.R-project.org/package=iterators

Wickham, H. 2011. The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.

Wickham, H. 2016. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York

Wickham, H. 2017. scales: Scale Functions for Visualization. R package version 0.5.0.
 https://CRAN.R-project.org/package=scales

Wickham, H. 2022. stringr: Simple, Consistent Wrappers for Common String Operations. R package version 1.5.0. https://CRAN.R-project.org/package=stringr

Wickham, H. & Henry, L. 2019. tidyr: Easily Tidy Data with 'spread()' and 'gather()' Functions. R package version 0.8.3 https://CRAN.R-project.org/package=tidyr

Wickham, H., & Hester, J. 2020. readr: Read Rectangular Text Data. R package version 1.4.0. https://CRAN.R-project.org/package=readr

Wickham, H., François, R., Henry, L. and Müller, K. 2018. dplyr: A Grammar of Data Manipulation. R package version 0.7.5. https://CRAN.R-project.org/package=dplyr

## Funding Acknowledgment

The creation of this study was supported by the National Science Foundation grant IOS-2309932 to Ivan Baxter and Brian Dilkes and Donald Danforth Plant Science Center (https://www.danforthcenter.org) funds supporting Lauren Whitt.  The datasets used in the work were collected under the funding of: National Science Foundation grant IOS-1126950 to Ivan Baxter, the United States Department of Energy Office of Science (BER; https://science.osti.gov/ber) grant DE-SC0023305 (Brian Dilkes), the Biotechnology and the Biological Sciences Research Council (https://www.ukri.org/councils/bbsrc/) grants BB/J003336/1 (Adam Price, Gareth Norton) and David Salt) and BB/L000113/1 (David Salt)
