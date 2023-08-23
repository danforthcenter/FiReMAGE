# Example test README

The example folder contains a small subset of the ionomic datasets to quickly check if FiReMAGE files were downloaded properly. Example_test/example_output/ is provided to compare with the test output. An example for running FiReMAGE in RStudio is available with these arguments:

  -s ./data/snps/ 
  -m ./data/metaTables/V1_3s_metaTable.csv 
  -f ./data/OrthoFinder_orthologs/Orthogroups_tutorial_subset.csv 
  -o ./test_output/ 
  -c [number of cpus] 
  -p 10
  
Permutations are expected to vary slightly as they are random and only n=10. If an exact comparison is desired, users can set perm_files on line 110 to ./example_test/example_output/. Nevertheless, all of the actual dataset results should be the same. An example of the verbose output from the code is provided: example_verbose_output.txt. This file is what a normal output (including all the warnings) looks like for a run of FiReMAGE.
