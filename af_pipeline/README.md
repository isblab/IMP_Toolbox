# AF-pipeline

This directory contains scripts to aid in alphafold related workflows

## AFinput
### Create job files for AF server
#### Description
- A general script to create job files for AF server.

**Input:** `.yaml` file in the following format.

```yaml
RNA_DNA_complex_8I54: # job cycle (required)
  - name: "Lb2Cas12a_RNA_DNA_complex" # job name  (not required)
    modelSeeds: [1,2] # 2 models with seeds 1 and 2 (not required)
    entities:
    # protein entity
      - name: "Lb2Cas12a" # (required)
        type: "proteinChain" # (required)
        count: 1
        glycans:
        - - "BMA"
          - 5
        modifications:
        - - "CCD_HY3"
          - 11
    # RNA entity
      - name: "RNA_33"
        type: "rnaSequence"
    # DNA entities
      - name: "DNA_25"
        type: "dnaSequence"
      - name: "DNA_mod"
        type: "dnaSequence"
        modifications: [["CCD_6OG", 2], ["CCD_6MA", 1]]
      - name: "MG"
        type: "ion"
        count: 1
```
- The only required keys are:
  1. job cycle
  2. name and type in entities
- For most use cases, the input will look like this:
```yaml
job_cycle:
# job 1
  - entities:
    - name: "protein_1"
      type: "proteinChain"
    - name: "protein_2"
      type: "proteinChain"
# job 2
  - entities:
    - name: "dna_1"
      type: "dnaSequence"
    - name: "protein_2"
      type: "proteinChain"
```
**Usage:**
- For allowed entity types as well as PTMs, ligands and ions, refer to `allowed_af_things.json` or [JSON file format for AlphaFold Server jobs](https://github.com/google-deepmind/alphafold/tree/main/server) 
- `modelSeeds` can either be an `int` or `list`.
  1. if `isinstance(modelSeeds, int)` -> `modelSeeds = random.sample(range(1, 10 * num_seeds), num_seeds)`
  2. if `isinstance(modelSeeds, list)` -> list taken as is

  Each seed in the list will be a new job.
- Input `yaml` file can contain multiple cycles, each with multiple jobs

```python
from af_pipeline.AFinput import AFInput

proteins = read_json("./input/proteins.json")
protein_sequences = read_fasta("./input/protein_sequences.fasta")
nucleic_acid_sequences = read_fasta("./input/nucleic_acid_sequences.fasta")
input_yml = yaml.load(open("./input/af_server_targets.yaml"), Loader=yaml.FullLoader)

af_input = AFInput(
    protein_sequences=protein_sequences, # required (output of fetch_sequences.py)
    input_yml=input_yml, # required
    nucleic_acid_sequences=nucleic_acid_sequences, # optional only in case of DNA or RNA sequences
    proteins=proteins, # optional if protein_sequences have protein names as headers and they match with input yaml
)

af_input.output_dir = args.output
job_cycles = af_input.create_job_cycles()
af_input.write_job_files(job_cycles=job_cycles)
```

Check the following examples in the examples directory for usage.
- `create_af_jobs.py`

## AFoutput
### Extract rigid bodies
Using Alphafold information extract high-confidence regions from Alphafold for use as rigid bodies in IMP
#### Description
- uses Tristan Crolls' clustering based on only PAE to form domains (less stringent because all PAE values in a domain need not be <5) and filters the domains based on per-residue pLDDT.

**Input:** `.yaml` file in the following format 
```yaml
- structure_path: "path/to/af_structure.cif"  # This is not required
  data_path: "path/to/af_data.json"  # This is required
  selection:  # This is not required
    - id: "A"
      model_region: [400,745]
      af_region: [1,745]
    - id: "B"
      model_region: [635,1118]
      af_region: [630,1118]
```

**Usage:**
- `srtucture_path` and `selection` is optional
- If only `data_path` is provided, the output will be only based on PAE cutoff.

```python
from af_pipeline.AFoutput import AFOutput

af = AFOutput.RigidBodies(data_path=data_path, output_dir="path/to/dir")
domains = af.predict_domains()
rigid_bodies = af.get_rigid_bodies(domains)
```

- To apply pLDDT filter, `structure_path` is needed (structure is used to get per-residue pLDDT)
```python
domains = af.predict_domains()
domains = af.plddt_filtered_domains(domains)
```

- `selection` can be provided if AF-prediction is not full-length or modeled region is not the whole AF region. Note that this is not required, by default the first residue in the AF-prediction will be assigned '1' in the output list
- Set `selected=True` in `get_rigid_bodies()` to apply selection.
- `num_proteins` is the minimum number of proteins to be present in a rigid body (if set to 2, rigid bodies that have at least 2 chains will be given as output)
```python
rigid_bodies = af.get_rigid_bodies(domains, selected=True, num_proteins=2)
```
- The extracted rigid bodies can be saved in `pdb` format.
```python
af.save_rigid_bodies_pdb(rigid_bodies)
```

#### Note
Currently the this method is preferred over `get_high_confidence_region_from_AF2.py`. That is, getting domains based on PAE first and then filtering for plDDDT is slightly better because,
the former script only considers the split of the domains at the borders of confident residue stretch (pLDDT>=70). So, a domain boundary within such a stretch will be missed.
e.g. as shown below, for alpha-actinin monomer, ABD-SR domains is within a confident residue-stretch.

- Output from `predict_rigid_bodies.py`: correct identification of domain boundary of ABD:SR
![image](https://github.com/user-attachments/assets/97cfe31a-e4af-4307-a033-b536c74b846f)
- Output from `get_high_confidence_region_from_af2.py`: missed domain boundary of ABD-SR
![image](https://github.com/user-attachments/assets/0902b16a-5683-46ec-9e92-e9379f28647b)

ABD: Actin-binding domain
SR: Spectrin repeat

Refer to:
- `extract_af_rigid_bodies.py` in `IMP_Toolbox/examples` for usage

#### Extract contacts or interface residues predicted confidently for use as restraints in IMP
WIP
