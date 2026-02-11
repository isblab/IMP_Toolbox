This is one way to get distance maps for compact complexes. It is currently slower (and not parallelized) compared to the alternative `contact_map.py`  

To identify the protein-protein interfaces, we compute regions of contact between each protein copy pair. These regions comprise bead pairs in contact; a bead pair is in contact if the average distance between their surfaces across models in the cluster is less than or equal to 10 Å. 
