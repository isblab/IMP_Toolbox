### Using Alphafold information

#### Extract high-confidence regions from Alphafold for use as rigid bodies in IMP 

Use one of these two methods
1. `get_high_confidence_region_from_af2.py` : uses plDDT and PAE to extract high-confidence regions. All residues with plDDT>70.0 are taken, and PAE domains from this subset are extracted to be modeled as rigid bodies in IMP. 
PAE domains correspond to residues with PAE values < 5 (all vs all residues).
This method is more stringent than the method below. 

2. `predict_domains script` : uses Tristan Crolls' clustering based on only PAE to form domains. Less stringent because all PAE values in a domain need not be <5.
Getting domains based on PAE first and then filtering for plDDDT is slightly better, because #TODO

#### Extract contacts or interface residues predicted confidently for use as restraints in IMP
See `https://github.com/isblab/af_pipeline/`
