from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model

def split_structure_by_chain(
    structure: Structure,
) -> dict[str, Structure]:
    """ Split a Biopython Structure object into separate Structure objects for
    each chain.

    ## Arguments:

    - **structure (Structure)**:<br />
        Biopython Structure object to be split by chain.

    ## Returns:

    - **dict[str, Structure]**:<br />
        A dictionary mapping chain IDs to new Structure objects containing only
        that chain.
    """

    chain_structures = {}
    model = structure[0]  # Get the first model

    for chain in model:
        chain_id = chain.get_id()
        new_structure: Structure = Structure(f"{structure.get_id()}_{chain_id}")
        new_structure.add(Model(0))
        new_structure[0].add(chain.copy())
        chain_structures[chain_id] = new_structure

    return chain_structures

def get_per_chain_residues(structure: Structure):

    chain_ranges = {
        chain.get_id(): [
            res.get_id()[1] for res in chain.get_residues()
        ] for chain in structure[0].get_chains()
    }

    return chain_ranges