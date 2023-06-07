### contact_map_all_pairs

Parallelized version of contact map for all protein pairs. Default code for contact maps.

**NOTE** It assumes that the protein sequence is numbered 1....n. Fix the code if you are modeling the sequence starting from the middle (i.e. without N-terminus)


### cm
Get boolean contact maps for all protein pairs. Parallelized. By Satwik Pasani.
Get the models in cluster 0 in one RMF for further analysis e.g. fit to data

Usage: \
~/imp-custom/build/setup_environment.sh python cm.py sampcon_0_extracted.rmf3 30
