# FiReMAGE: Filtering Results of Multispecies Analogous GWAS Experiments

FiReMAGE is a tool written in R that analyzes GWAS datasets of the same traits in multiple species to search for orthologous genes associated with SNPs from their respective species' dataset. 

R version 3.6.0 (2019-04-26)

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
	10. rbenchmark (citation)
	11. docopt - optional for RStudio use, necessary for use in terminal R (Jonge 2016)

## 2. Directories and Accompanying files
FiReMAGE has two separate versions for OrthoFinder and InParanoid orthologs, each directory has the respective files required to run each. Both share the same SNP data. In either version you will need:

	1. Pipeline Rscript: FiReMAGE_[OrthoFinder|InParanoid]_[local|server]_version.R - both directories have examples for running FiReMAGE locally in RStudio or in terminal remotely. 
	2. functions Rscript: functions_[OrthoFinder|InParanoid]_version.R - Rscript containing all the functions used in the pipeline, function version should match pipeline version.
	3. org_chromosome_coords.csv - located in the data subdirectory, contains all of the base pair lengths for each chromosome of each species in the comparison
	4. metaTable csv - template located in data/metaTable directory
	5. ortholog files - for InParanoid, these will be in the data/phytozome directory, for OrthoFinder data/OrthoFinder_orthologs
	6. gene coordinate files - in the data/phytozome directory, contains gene locations for all genes in each species
	7. GWAS SNPs for the species in each comparison

## 3. Function
FiReMAGE v.1.0 comes with the RScript files for both OrthoFinder and InParanoid versions of the pipeline, and the data directory containing organism chromosome coordinates, Phytozome files, ortholog calls, and SNPs for Arabidopsis thaliana, Glycine max, Oryza sativa, and Sorghum bicolor and Zea mays. The versions of these included organisms are listed in currentVersions.txt, included in the data directory. If the SNP data is different from the most up-to-date versions, or other deviations like custom coordinates are to be used, they should be altered in the organism coordinates file. When using real data SNP files, make sure the format matches the test SNP files in the data/seedIonome_snps/ directory, and that the name includes the species name as spelled in the metaTable. This allows for multiple SNP files to be kept in one directory, and only the ones included in the metaTable will be read into the pipeline. Additional Phytozome files for species not included in FiReMAGE v.1.0 can be added by rerunning the OrthoFinder comparison or by downloading InParanoid files from Phytozome and adding the new files to their respective directory. For InParanoid, you will need a phytozome file for each unique comparison of every species in your analysis. For example, if you have an A.thaliana_to_G.max phytozome file, you do not need a G.max_to_A.thaliana phytozome file.

The FiReMAGE pipeline starts out by creating the metaTable, which will indicate the species in the comparison, the linkage disequilibrium (LD) ranges for each species, and the number of chromosomes for each species. It's very important that the metaTable be in the same format as the template because the following functions all refer to the metaTable for information on how many times to iterate, which species is next, etc. 

The next step is to read in the snp files. The function snpFilePrep will read the files in, sort them by trait, organism if applicable, and chromosome and base pair, and combine all the files into one table. From here, the pipeline makes note of all the unique traits in the SNP file, and then collapses SNPS into loci. Two or more SNPs will be condensed into one loci if their LD ranges overlap. Unique names are generated for each loci by their organism name, trait, chromosome and a number. The number of SNPs per trait, per organism is then recorded in the table leafcollapsedSNPs.csv, saved to the output directory. If your SNPs are already collapsed into loci, make sure to mark input_as_loci==TRUE in the metaTable to keep your custom loci ranges.

If using InParanoid, Phytozome comparison files are read in for each unique, "choose-two" comparison and merged into an object called "genefile". This includes all ortholog groups representing 3 or more species. Ortholog groups only containing orthologs between two species will be thrown out. The "present" column shows the number of species from the analysis that are in the ortholog group. If using OrthoFinder, the ortholog file was made in OrthoFinder v2.0 and will be imported from the data/OrthoFinder_orthologs/ directory.

The ortholog genefile is then compared to the SNP file. If the coordinates of a gene and loci overlap all organisms, the genes are saved as an "n-way hit" along with the loci in each species in a file named  Orthogrup_hits.csv written to the output directory. Trait files can include traits not represented in all species, as these will get thrown out later, but it will slow down time in the process as it will iterate through these SNP rows. They will only be thrown out if the number of species representing the trait is below the threshold or is zero in the graph outputs later on. Overall, it's recommended to use SNP files all using the same traits, but it won't prevent the pipeline from working (unless you include a file with SNPs from 20 traits not included in other species - that will just create a lot of extra work to be thrown out later).

To assess whether the frequency of returned genes is due to noise in the genome, random permutations are used to measure the chance we can generate the same results by selecting SNPs throughout the genome at random. Using the leafcollapsedSNP.csv, the pipeline will generate the permutation datasets, randomly selecting the same number of loci from each chromosome and organism as the actual dataset. Random ranges are generated for the permutation SNPs to simulate collapsing. Because the actual dataset represents collapsed SNPs with variable ranges due to the amount of nearby SNPs, we need to replicate this in the permutation dataset. The final permutation dataset with the SNPs and random ranges is saved in the output directory permutation_files/snps/. The permutations are then put into the pipeline again, this time aggregating the set by permutation as well as species and trait. It's corresponding files are in permutation_files/gene_hits/. The random datasets would be too large to display like the real SNP datasets, so these are already summarized in the function and stored in permutation_files/. The AllSummaries.csv gives the summary for the real dataset and permutations (ie. Number of genes and loci in each species combo for each trait) where the average and quantiles of the number of genes and permutations for each trait and species calculated.

A list of the genes found in the actual dataset are listed and ranked in a file called priorityList.csv. It has the genes from all the species and will rank them using pmax() to sort, in parallel, the genes by their false discovery rate, SNP p-value and species combo that returned the gene.

Data from the actual and permutation datasets are displayed as graphs in the OUTDIRPATH/graphs/ directory.  Horizontal barplots display the number of orthologous genes per trait (95 and 99th percentiles) - one frame per species. Preferences for graph output and design can be altered in the main FiReMAGE pipeline Rscript. 

## 4. Options

	1. -h, --help : Shows a docopt help screen displaying all options
	2. -s <path>, --snps <path> : The path to the directory containing the SNP files for the comparison, required for the operation of FiReMAGE. The pathway for the test/tutorial dataset is “./data/el_snps/current/”
	3. -m <path>, --metaTable <path>: Path to read in the metadata table of organisms in the comparison, required for the operation of FiReMAGE. A template and sample metatable input is in ./data/metaTables/
	4. -o <path>, --output <path>: Files will be written to this path, the default is “./FM_output/”
	5. -p <int>, --permutations <int>: Number of permutations for random dataset, default and recommendation is 1000 for full analysis, 100 for initial testing
	6. -c <int>, --cores <int>: Number of cores to use for parallel loops, the default is 1 so it will naturally run on any system
	7. -f <path> --orthologFiles: The path to the directory containing ortholog files.

## 5. SNP input file format

To work with the FiReMAGE pipeline, SNP files should have columns for the organism, trait, base pair position on chromosome, and chromosome number of the SNPs. They do not necessarily have to be in order, as the pipeline with order the file with the function snpFilePrep. Chromosome and base pair columns should be numerical (ie. 7, instead of Chr 7, Chr_7). The pipeline expects a different SNP file for each organism in the comparison; no need to concatenate all the files before had.

## 6. Customizing ortholog files
OrthoFinder: Users should downloaded the appropriate fasta protien files and run according to OrthoFinder v2.0 instructions (https://github.com/davidemms/OrthoFinder, Emms & Kelly, 2019).
InParanoid:New ortholog comparison files downloaded from Phytozome should have no problem being read into FiReMAGE. The function geneFilePrep searches for, and renames, all the column headers generated by Phytozome when exporting tables from their website (when exporting, select the option for adding “human-readable” column headers). Make sure to select a file for all homologs genes of two specific species when generating this file, and add columns for the chromosome, gene start, and gene end of both the ortholog and homolog. We've had trouble downloading even the smallest comparison file from PhytoMine's main web browser, if you are having the same issues try downloading it through their biomart interface (https://phytozome.jgi.doe.gov/biomart/martview). 

## 7. Test data

FiReMAGE comes with a shortened version of our four-way comparison between A. thaliana, G. max, S. bicolor, and Z. mays SNP files as a way to test the set up of FiReMAGE on your local computer. The SNP files are labeled for their species with the suffix _test_data.csv. Unpack, and install the necessary libraries, and it's ready to run (RStudio and command line available). An example of the command line version: 

Rscript FiReMAGE.R -s "./data/seedIonome_snps/" -m "./data/metaTables/CGASmetaTableInput.csv"  -f  "./data/OrthoFinder_orthologs/" -p 1000 -c 8 

The resulting files should be output by default into ./FM_output/. This can be compared to the directory ‘./test_output/’ which has the example files that should be produced from this test. The Orthogroup_hits.csv should be the same, and the permutation files may be slightly different due to the random distributions. If an exact comparison is desired, users can look at the set.seed() function in R. Nevertheless, all of the actual dataset results should be the same. An example of the verbose output from the code can be seen in ‘./test_output/test_log.txt’. This file is what a normal output (including all the warnings) looks like for a run of FiReMAGE. 

## References

Cooper, N. 2017. reader: Suite of Functions to Flexibly Read Data from Files. R package version 1.0.6. https://CRAN.R-project.org/package=reader

De Jonge, E. 2016. docopt: Command-Line Interface Specification Language. R package version 0.4.5.  https://CRAN.R-project.org/package=docopt

Dowle, M. & Srinivasan, A. 2018. data.table: Extension of `data.frame`. R package version 1.11.4. https://CRAN.R-project.org/package=data.table

Emms, D. M., & Kelly, S. 2015. OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome biology, 16(1), 1-14.

Kusnierczyk, W. 2012. rbenchmark: Benchmarking routine for R. R package version 1.0.0. https://CRAN.R-project.org/package=rbenchmark

Microsoft Corporation & Weston, S. 2017. doParallel: Foreach Parallel Adaptor for the 'parallel' Package. Rpackage version 1.0.11. https://CRAN.R-project.org/package=doParallel

Revolution Analytics and Weston, S. 2017. iterators: Provides Iterator Construct for R. R package version 1.0.9. https://CRAN.R-project.org/package=iterators

Wickham, H. 2011. The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.

Wickham, H. 2016. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York

Wickham, H. 2017. scales: Scale Functions for Visualization. R package version 0.5.0.
 https://CRAN.R-project.org/package=scales 

Wickham, H., François, R., Henry, L. and Müller, K. 2018. dplyr: A Grammar of Data Manipulation. R package version 0.7.5. https://CRAN.R-project.org/package=dplyr

Wickham, H., & Hester, J. 2020. readr: Read Rectangular Text Data. R package version 1.4.0. https://CRAN.R-project.org/package=readr

Wickham, H., Hester, J. & Francois, R. 2017. readr: Read Rectangular Text Data. R package version 1.1.1. https://CRAN.R-project.org/package=readr

Wickham, H. & Henry, L. 2019. tidyr: Easily Tidy Data with 'spread()' and 'gather()' Functions. R package version 0.8.3 https://CRAN.R-project.org/package=tidyr

Ziegler, G. 2017. ionomicsUtils: Utility functions for ionomics data analysis. R package version 1.0.
