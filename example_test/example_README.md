#Example test README

The example folder contains a small subset of the ionomic datasets to quickly check if FiReMAGE files were downloaded properly. Example_test/example_output/ is provided to compare with the test output. A script for checking on a local device in RStudio is in example_test/scripts/, for testing on a remote server run the FiReMAGE_OrthoFinder_server_version.R with these arguments:

  -s ./data/snps/ 
  -m ./data/metaTables/V1_3s_metaTable.csv 
  -f ./data/OrthoFinder_orthologs/Orthogroups_tutorial_subset.csv 
  -o ./test_output/ 
  -c [number of cpus] 
  -p 10
  
Permutations are expected to vary slightly as they are random and only n=10. If an exact comparison is desired, users can look at the set.seed() function in R. Nevertheless, all of the actual dataset results should be the same. An example of the verbose output from the code can be seen in ‘./test_output/test_log.txt’. This file is what a normal output (including all the warnings) looks like for a run of FiReMAGE.
