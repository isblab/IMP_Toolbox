"""
structure
===========

- Assisting with the structure related tasks such as-
  - fetching best structures for a given uniprot id from PDB,
  - calculating buried residue information (e.g. residue depth)
  - splitting structure by chain and obtaining per chain residues
  - transforming and saving structure objects
  - converting between PDB and mmCIF formats

<hr>

### Examples:

<hr>

#### 1. Find best structures for given uniprot ids

```python
from IMP_Toolbox.structure.best_structure import (
    fetch_best_structures,
    make_best_structures_df
)

uniprot_ids = ["P12345", "Q67890"]
best_structures = fetch_best_structures(
    uniprot_ids=uniprot_ids,
    save_path="./output/best_structures.json",
    overwrite=False
)
best_structures_df = make_best_structures_df(best_structures)
best_structures_df.to_csv("./output/best_structures.csv", index=False)
```

#### 2. Get modeled coverage of a chain in a structure given a PDB id

```python
from IMP_Toolbox.structure.best_structure import calculate_modeled_coverage
from IMP_Toolbox.sequence.sequence import query_uniprot_api_for_lengths

pdb_id = "1J6Z"
chain_id = "A"
uniprot_id = "P68135"
uniprot_length = query_uniprot_api_for_lengths(
    uniprot_ids=[uniprot_id]
)[uniprot_id]

coverage, modeled_residues = calculate_modeled_coverage(
    pdb_id=pdb_id,
    chain_id=chain_id,
    uniprot_length=uniprot_length,
    outdir="./output"
)
print(f"Modeled coverage: {coverage}, Modeled residues: {modeled_residues}")
```

#### 3. Get burial information for a structure given a PDB/mmCIF file

```python
from Bio.PDB import PDBList, PDBParser
from IMP_Toolbox.structure.burial import get_burial_info

pdbl = PDBList()
# Download the PDB file
pdb_file = pdbl.retrieve_pdb_file("1J6Z", pdir="./output", file_format="pdb")

burial_info_df = get_burial_info(structure_path=pdb_file)
print(burial_info_df.head())
```

- If you want to include residue depth information, you can set the `include_residue_depth`
  parameter to `True` and provide the path to the MSMS executable for calculating residue depth.:

```python
burial_info_df = get_burial_info(
    structure_path=pdb_file,
    include_residue_depth=True,
    msms_executable="/path/to/msms.x86_64Linux2.2.6.1"
)
print(burial_info_df.head())
```

- There are additional parameters available in the `get_burial_info` function to customize
  the output, such as `residue_selector`, `entity_chain_map`, and `include_residue_depth`.
  You can refer to the function's docstring for more details. Example, to get information
  only for a subset of residues in a specific chain, you can use the `residue_selector` parameter:

```python
# Only include residues 170, 171, and 172 from chain A
residue_selector = {"A": [170, 171, 172]}
burial_info_df = get_burial_info(
    structure_path=pdb_file,
    residue_selector=residue_selector,
    include_residue_depth=True,
    msms_executable="/path/to/msms.x86_64Linux2.2.6.1"
)
print(burial_info_df.head())
```

#### 4. Split a structure by chain and obtain per chain residues

```python
from Bio.PDB.mmtf import MMTFParser
from IMP_Toolbox.structure.tools import split_structure_by_chain, get_per_chain_residues

structure = MMTFParser.get_structure_from_url("1J6Z")
chain_structures = split_structure_by_chain(structure)
chain_residues = get_per_chain_residues(structure)
```

#### 5. Convert between PDB and mmCIF formats

```python
from IMP_Toolbox.structure.tools import pdb_to_mmcif, mmcif_to_pdb
pdb_to_mmcif(
    input_pdb="./path/to/input.pdb",
    output_mmcif="./path/to/output.cif"
)
mmcif_to_pdb(
    input_mmcif="./path/to/input.cif",
    output_pdb="./path/to/output.pdb"
)
```

#### 6. Transform and save structure objects

```python
import numpy as np
from IMP_Toolbox.structure.tools import transform_pdb

# here, the transform_matrix is a 4x3 transformation matrix that defines
# the rotation and translation to be applied to the structure. The first three rows
# represent the rotation matrix, and the last row represents the translation vector.
# You can obtain such a matrix from US-align output
transform_matrix=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [20, 30, 40]])

transform_pdb(
    pdb_file="./path/to/input.pdb",
    out_path="./path/to/transformed.pdb",
    transform_matrix=transform_matrix,
)
```

#### 7. Renumber and save structure objects

```python
from Bio.PDB.mmtf import MMTFParser
from IMP_Toolbox.structure.tools import RenumberResidues, save_structure_obj

structure = MMTFParser.get_structure_from_url("1J6Z")
renumber = RenumberResidues(offset={"A": 101})
renumbered_structure = renumber.renumber_structure(structure=structure)
save_structure_obj(
    structure=renumbered_structure,
    out_file="./path/to/renumbered.pdb",
    save_type="pdb",
)
```

- `res_select_obj` is an optional parameter that allows you to specify a selection
of residues to save. If provided, only the selected residues will be saved in the
output file. You can use the `Bio.PDB.Select` class to create a custom selection
object based on your criteria.
"""

from .burial import (
    get_burial_info,
    get_residue_depth,
)

from .best_structure import (
    get_best_structures,
    fetch_best_structures,
    make_best_structures_df,
)

from .tools import (
    split_structure_by_chain,
    get_per_chain_residues,
    pdb_to_mmcif,
    mmcif_to_pdb,
    transform_pdb,
    save_structure_obj,
)