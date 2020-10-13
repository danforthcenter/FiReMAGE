# FiReMAGE: Filtering Results of Multispecies Analogous GWAS Experiments

FiReMAGE is a tool written in R that analyzes GWAS datasets of the same traits in multiple species to search for orthologous genes associated with SNPs from their respective species' dataset. 

R version 3.5.0 (2018-04-23)

	1. R Packages
	2. Accompanying files
	3. Function
	4. Options
	5. SNP input file format
	6. Customizing ortholog files
	7. Test data

## 1. R Packages
To run CGAS, users will have to pre-install the following packages
	1. plyr and dplyr (Wickham 2011 & Wickham et al. 2018)
	2. tidyr (Wickham 2019)
	3. data.table (Dowle & Srinivasan 2018)
	4. doParallel (Microsoft & Weston 2017)
	5. reader (Cooper 2017)
	6. iterators (Revolution Analytics & Weston 2017) 
	7. scales (Wickham 2017)
	8. ggplot2 (Wickham 2016)
	9. ionomicsUtils - optional; for ionome data (Ziegler 2017)
	10. docopt - optional for RStudio use, necessary for use in terminal R (Jonge 2016)

## 2. Accompanying files
In its most “barebones” state, FiReMAGE needs the following files to operate:
	1. FiReMAGE.R - pipeline
	2. functions.R - Rscript containing all the functions used in the pipeline
	3. org_chromosome_coords.csv - contains all of the base pair lengths for each chromosome of each species in the comparison
	4. metaTable data - template located in data/metaTable directory
	5. Phytozome InParanoid files - files for test dataset located in the data/phytozome/ directory
	6. SNPs for the species in each comparison

## 3. Function
FiReMAGE v.0.5 comes with the RScript files FiReMAGE.R and functions.R, and the data directory containing organism chromosome coordinates, Phytozome InParanoid files for Arabidopsis thaliana, Glycine max, Sorghum bicolor and Zea mays. The versions of these included organisms are listed in currentVersions.txt, included in the data directory. If the SNP data is different from the most up-to-date versions, or other deviations like custom coordinates are to be used, they should be altered in the organism coordinates file. When using real data SNP files, make sure the format matches the test SNP files in the data/el_snps/ directory, and that the name includes the species name as spelled in the metaTable. This allows for multiple SNP files to be kept in one directory, and only the ones included in the metaTable will be read into the pipeline. Additional Phytozome files for species not included in FiReMAGE v.0.5 should be downloaded from Phytozome with header orders matching the test files in the data/phytozome/ directory, and should be added to the data/phytozome/ directory. You will need a phytozome file for each unique comparison of every species in your analysis. For example, if you have an A.thaliana_to_G.max phytozome file, you do not need a G.max_to_A.thaliana phytozome file.

The beginning lines of the pipeline changes slightly depending on if it is run in terminal or in GUI RStudio. When running the pipeline on R in a terminal setting, it is recommended to use the docopt formatting of options (lines 1-40) rather than the regular way of declaring variables in RStudio (lines 42-51). By default, CGASconcise.R comes with the terminal docopt code commented out. To run in terminal, one should uncomment lines 1-40 and comment lines 42-50 to prevent the RStudio code from overwriting the terminal options. 

The FiReMAGE pipeline starts out by creating the metaTable, which will indicate the species in the comparison, the linkage disequilibrium (LD) ranges for each species, and the number of chromosomes for each species. It's very important that the metaTable be in the same format as the template because the following functions all refer to the metaTable for information on how many times to iterate, which species is next, etc. 

The next step is to read in the snp files. The function snpFilePrep will read the files in, sort them by trait, organism if applicable, and chromosome and base pair, and combine all the files into one table. From here, the pipeline makes note of all the unique traits in the SNP file, and then collapses SNPS into loci. Two or more SNPs will be condensed into one loci if their LD ranges overlap. The number of SNPs per trait, per organism is then recorded in the table leafcollapsedSNPs.csv, saved to the output directory. 

Next, Phytozome comparison files are read in for each unique, "choose-two" comparison and merged into an object called "genefile". This includes all ortholog groups representing 3 or more species. Ortholog groups only containing orthologs between two species will be thrown out. The "present" column shows the number of species from the analysis that are in the ortholog group.

The ortholog genefile is then compared to the SNP file. If the coordinates of a gene and loci overlap all organisms, the genes are saved as an "n-way hit" along with the loci in each species in a file named  FinalMerge.csv written to the output directory. Unique names are generated for each loci by their organism name, trait, chromosome and a number. Trait files can include traits not represented in all species, as these will get thrown out later, but it will slow down time in the process as it will iterate through these SNP rows. They will only be thrown out if the number of species representing the trait is below the threshold or is zero in the graph outputs later on. Overall, it's recommended to use SNP files all using the same traits, but it won't prevent the pipeline from working (unless you include a file with SNPs from 20 traits not included in other species - that will just create a lot of extra work to be thrown out later).

It will also save to the output directory csv files containing the loci hitting orthologous genes and how many genes the loci returned, the orthologous genes in AllLoci.csv. To assess whether the frequency of returned genes is due to noise in the genome, random permutations are used to measure the chance we can generate the same results by selecting SNPs throughout the genome at random. Using the leafcollapsedSNP.csv, the pipeline will generate the permutation datasets, randomly selecting the same number of loci from each chromosome and organism as the actual dataset. Random ranges are generated for the permutation SNPs to simulate collapsing. Because the actual dataset represents collapsed SNPs with variable ranges due to the amount of nearby SNPs, we need to replicate this in the permutation dataset. The final permutation dataset with the SNPs and random ranges is saved in the output directory as AllPermuts.csv. The permutations are then put into the pipeline again, this time aggregating the set by permutation as well as species and trait. It's corresponding files are AllPermutationfiles.csv and PermGeneCount.csv. The random datasets would be too large to display like the FinalMerge.csv, so these are already summarized in the function. The AllPermutationfiles.csv gives the summary for each permutation (ie. Number of genes and loci in each species combo for each trait), setting them up to be put into the graphingDF function to have the average and quantiles of the number of genes and permutations for each trait and species calculated. The PermGeneCount.csv is setting up calculations for the multigene/loci table produced later. This is the number of loci that hit 0, 1, 2, etc. number of genes in each random permutation dataset. The permutation column must stay here because the multigene/loci graph will take the frequency of permutations that return these GeneCounts.

A list of the genes found in the actual dataset are listed and ranked in a file called priorityList.csv. It has the genes from all the species and will rank them using pmax() to sort, in parallel, the genes by their false discovery rate, SNP p-value and species combo that returned the gene.

Data from the actual and permutation datasets are displayed as graphs in the OUTDIRPATH/graphs directory.  Horizontal barplots display the number of orthologous genes per trait (95 and 99th percentiles) - one frame per species and the number of loci returning orthologous genes per trait (95 and 99th percentiles) - one frame per species. Vertical barplots display the number of loci associated with multiple orthologous genes - one graph per trait, per species. Preferences for graph output and design can be altered  in FiReMAGE.R. 

## 4. Options
	1. -h, --help : Shows a docopt help screen displaying all options
	2. -s <path>, --snps <path> : The path to the directory containing the SNP files for the comparison, required for the operation of CGAS. The pathway for the test/tutorial dataset is “./data/el_snps/current/”
	3. -m <path>, --metaTable <path>: Path to read in the metadata table of organisms in the comparison, required for the operation of CGAS. A template and sample metatable input for the test/tutorial dataset is “./data/metaTables/CGASmetaTableInput.csv”
	4. -o <path>, --output <path>: Files will be written to this path, the default is “./CGAS_output/”
	5. -p <int>, --permutations <int>: Number of permutations for random dataset, default and recommendation is 1000
	6. -c <int>, --cores <int>: Number of cores to use for parallel loops, the default is 1 so it will naturally run on any system
	7. -f <path> : The path to the directory containing the Phytozome ortholog files. The pathway for the test/tutorial dataset is "./data/phytozome/"

## 5. SNP input file format

To work with the FiReMAGE pipeline, SNP files should have columns for the organism, trait, base pair position on chromosome, and chromosome number of the SNPs. They do not necessarily have to be in order, as the pipeline with order the file with the function snpFilePrep. Chromosome and base pair columns should be numerical (ie. 7, instead of Chr 7, Chr_7). The pipeline expects a different SNP file for each organism in the comparison; no need to concatenate all the files before had.

## 6. Customizing ortholog files
New ortholog comparison files downloaded from Phytozome should have no problem being read into FiReMAGE. The function geneFilePrep searches for, and renames, all the column headers generated by Phytozome when exporting tables from their website (when exporting, select the option for adding “human-readable” column headers). Make sure to select a file for all homologs genes of two specific species when generating this file, and add columns for the chromosome, gene start, and gene end of both the ortholog and homolog. We've had trouble downloading even the smallest comparison file from PhytoMine's main web browser, if you are having the same issues try downloading it through their biomart interface (https://phytozome.jgi.doe.gov/biomart/martview). 

## 7. Test data

FiReMAGE comes with a shortened version of our four-way comparison between A. thaliana, G. max, S. bicolor, and Z. mays SNP files as a way to test the set up of FiReMAGE on your local computer. The SNP files are labeled for their species with the suffix _test_data.csv. Unpack, and install the necessary libraries, and it's ready to run in RStudio console. To run in the terminal, you first need to uncomment the docopt options for use in terminal (lines 31-38) and comment out the editable options for use in Rstudio GUI (41-46). Then in the command line run: 

Rscript FiReMAGE.R -s "./data/el_snps/current/" -m "./data/metaTables/CGASmetaTableInput.csv"  -f  "./data/phytozome/current/" -p 10 -c 1 

The resulting files should be output by default into ./FM_output/. This can be compared to the directory ‘./test_output/’ which has the example files that should be produced from this test. The FinalMerge.csv should be the same, and RandomMerge.csv may be slightly different due to it’s random nature. If an exact comparison is desired, the set.seed lines can be uncommented in the Random_searchNEW.R section on lines 163 and 165 to get the exact results in the ‘./test_output/’ directory. Nevertheless, all of the actual dataset results should be the same. An example of the verbose output from the code can be seen in ‘./test_output/test_log.txt’. This file is what a normal output (including all the warnings) looks like for a run of FiReMAGE. 

References:
Cooper, N. 2017. reader: Suite of Functions to Flexibly Read Data from Files. R package version 1.0.6. https://CRAN.R-project.org/package=reader

De Jonge, E. 2016. docopt: Command-Line Interface Specification Language. R package version 0.4.5.  https://CRAN.R-project.org/package=docopt

Dowle, M. & Srinivasan, A. 2018. data.table: Extension of `data.frame`. R package version 1.11.4. https://CRAN.R-project.org/package=data.table

Microsoft Corporation & Weston, S. 2017. doParallel: Foreach Parallel Adaptor for the 'parallel' Package. Rpackage version 1.0.11. https://CRAN.R-project.org/package=doParallel

Revolution Analytics and Weston, S. 2017. iterators: Provides Iterator Construct for R. R package version 1.0.9. https://CRAN.R-project.org/package=iterators

Wickham, H. 2011. The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.

Wickham, H. 2016. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York

Wickham, H. 2017. scales: Scale Functions for Visualization. R package version 0.5.0.
 https://CRAN.R-project.org/package=scales 

Wickham, H., François, R., Henry, L. and Müller, K. 2018. dplyr: A Grammar of Data Manipulation. R package version 0.7.5. https://CRAN.R-project.org/package=dplyr

Wickham, H., Hester, J. & Francois, R. 2017. readr: Read Rectangular Text Data. R package version 1.1.1. https://CRAN.R-project.org/package=readr

Wickham, H. & Henry, L. 2019. tidyr: Easily Tidy Data with 'spread()' and 'gather()' Functions. R package version 0.8.3 https://CRAN.R-project.org/package=tidyr

Ziegler, G. 2017. ionomicsUtils: Utility functions for ionomics data analysis. R package version 1.0.
