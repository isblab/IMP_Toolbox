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