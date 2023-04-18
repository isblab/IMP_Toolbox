### extract_sampcon
Get the models in cluster 0 in one RMF for further analysis e.g. fit to data 
Usage: 
~/imp-custom/build/setup_environment.sh python extract_sampcon.py sampcon_0_extracted.rmf3 ../A_gsm_clust3.rmf3 cluster.0.sample_A.txt ../B_gsm_clust3.rmf3 cluster.0.sample_B.txt 


### cm
Get boolean contact maps for all protein pairs. Parallelized. By Satwik Pasani. 
Usage:
~/imp-custom/build/setup_environment.sh python cm.py sampcon_0_extracted.rmf3


