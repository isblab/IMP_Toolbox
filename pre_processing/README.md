### Using Alphafold information

#### Extract high-confidence regions from Alphafold for use as rigid bodies in IMP 

Use one of these two methods
1. `get_high_confidence_region_from_af2.py` : uses plDDT and PAE to extract high-confidence regions. All residues with plDDT>70.0 are taken, and PAE domains from this subset are extracted to be modeled as rigid bodies in IMP. 
PAE domains correspond to residues with PAE values < 5 (all vs all residues).
This method is more stringent than the method below. 

2. `predict_rigid_bodies.py` : uses Tristan Crolls' clustering based on only PAE to form domains (less stringent because all PAE values in a domain need not be <5) and filters the domains based on per-residue pLDDT.
=======

Currently the latter method is preferred. That is, getting domains based on PAE first and then filtering for plDDDT is slightly better because,
the former script only considers the split of the domains at the borders of confident residue stretch (pLDDT>=70). So, a domain boundary within such a stretch will be missed.
e.g. as shown below, for alpha-actinin monomer, ABD-SR domains is within a confident residue-stretch.

- Output from `predict_rigid_bodies.py`: correct identification of domain boundary of ABD:SR
![image](https://github.com/user-attachments/assets/97cfe31a-e4af-4307-a033-b536c74b846f)
- Output from `get_high_confidence_region_from_af2.py`: missed domain boundary of ABD-SR
![image](https://github.com/user-attachments/assets/0902b16a-5683-46ec-9e92-e9379f28647b)

ABD: Actin-binding domain
SR: Spectrin repeat

#### Extract contacts or interface residues predicted confidently for use as restraints in IMP
See `https://github.com/isblab/af_pipeline/`
