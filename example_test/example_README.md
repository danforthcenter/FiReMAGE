#Example test README

The example folder contains a small subset of the ionomic datasets to quickly check if FiReMAGE files were downloaded properly. Example_test/example_output/ is provided to compare with the test output. Permutations are expected to vary slightly as they are random and only n=10. A script for checking on a local device in RStudio is in example_test/scripts/, for testing on a remote server run the FiReMAGE_OrthoFinder_server_version.R with these arguments:

  -s ./data/snps/ 
  -m ./data/metaTables/V1_3s_metaTable.csv 
  -f ./data/OrthoFinder_orthologs/Orthogroups_tutorial_subset.csv 
  -o ./test_output/ 
  -c [number of cpus] 
  -p 10