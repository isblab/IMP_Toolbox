import os
import pandas as pd
from tqdm import tqdm
from IMP_Toolbox.sequence.sequence import query_uniprot_api_for_lengths
from IMP_Toolbox.utils.api_helpers import (
    request_session,
    request_result
)
from IMP_Toolbox.utils.file_helpers import (
    read_json,
    write_json
)
from IMP_Toolbox.constants.imp_toolbox_constants import (
    APIurl,
)
from IMP_Toolbox.constants.structure_constants import BestStructureCol
from IMP_Toolbox.constants.sequence_constants import UNIPROT_ISOFORM_SEPARATOR

def make_best_structures_df(
    best_structures: dict,
    protein_uniprot_map: dict = {},
) -> pd.DataFrame:
    """ Create a dataframe of best structures from dictionary

    ## Arguments:

    - **best_structures (dict)**:<br />
        The best structures dictionary returned by `get_best_structures`.

    - **uniprot_protein_map (dict, optional):**:<br />
        Mapping from uniprot ids to protein names. Defaults to {}.

    ## Returns:

    - **pd.DataFrame**:<br />
        A dataframe containing the best structures for the given uniprot ids.
    """

    all_best_chains = []

    for uniprot_id, best_structure in best_structures.items():

        uniprot_base = uniprot_id.split(UNIPROT_ISOFORM_SEPARATOR)[0]

        if not isinstance(best_structure, dict):
            print(f"No best structure found for {uniprot_id}")
            continue

        best_chains = best_structure[uniprot_base]
        best_chains = [
            chain for chain in best_chains
            if chain[BestStructureCol.RESOLUTION] is not None
        ]
        best_chains = sorted(
            best_chains, key=lambda x: x[BestStructureCol.RESOLUTION]
        )

        for _, chain in enumerate(best_chains):

            all_best_chains.append(
            {
                BestStructureCol.UNIPROT_ID: uniprot_id,
                BestStructureCol.PROTEIN: protein_uniprot_map.get(uniprot_id, ""),
                BestStructureCol.CHAIN_ID: chain[BestStructureCol.CHAIN_ID],
                BestStructureCol.PDB_ID: chain[BestStructureCol.PDB_ID],
                BestStructureCol.COVERAGE: chain[BestStructureCol.COVERAGE],
                BestStructureCol.MODELED_COVERAGE: chain.get(BestStructureCol.MODELED_COVERAGE, ""),
                BestStructureCol.MODELED_RESIDUES: chain.get(BestStructureCol.MODELED_RESIDUES, ""),
                BestStructureCol.UNP_START: chain[BestStructureCol.UNP_START],
                BestStructureCol.UNP_END: chain[BestStructureCol.UNP_END],
                BestStructureCol.START: chain[BestStructureCol.START],
                BestStructureCol.END: chain[BestStructureCol.END],
                BestStructureCol.RESOLUTION: chain[BestStructureCol.RESOLUTION],
            }
        )

    df = pd.DataFrame(all_best_chains)
    mask = df.duplicated(subset=[BestStructureCol.UNIPROT_ID], keep='first')
    df[BestStructureCol.UNIPROT_ID] = df[BestStructureCol.UNIPROT_ID].mask(mask, '')
    df[BestStructureCol.PROTEIN] = df[BestStructureCol.PROTEIN].mask(mask, '')
    df = df.drop_duplicates(subset=[BestStructureCol.PDB_ID], keep='first')

    return df

def get_best_structures(uniprot_id: str) -> dict:
    """ Get the best structures for a given uniprot id.
    See the correspondoing fetch function :func:`fetch_best_structures`

    ## Arguments:

    - **uniprot_id (str)**:<br />
        The uniprot id for which to fetch the best structures.

    ## Returns:

    - **dict**:<br />
        A dictionary containing the best structures for the given uniprot id.
    """

    req_sess = request_session(max_retries=3)
    get_request = req_sess.get(
        url=APIurl.pdbe_api_best_structures.substitute(uniprot_id=uniprot_id)
    )

    best_structures = request_result(
        get_request, uniprot_id, ignore_error=True
    )

    return best_structures

def get_polymer_coverage(pdb_id: str, chain_id: str) -> dict:
    """ Get the polymer coverage for a given pdb id and chain id.
    See the correspondoing fetch function :func:`fetch_polymer_coverage`

    ## Arguments:

    - **pdb_id (str)**:<br />
        The pdb id for which to fetch the polymer coverage.

    - **chain_id (str)**:<br />
        The chain id for which to fetch the polymer coverage.

    ## Returns:

    - **dict**:<br />
        A dictionary containing the polymer coverage for the given pdb id and chain id.
    """

    req_sess = request_session(max_retries=3)
    get_request = req_sess.get(
        url=APIurl.pdbe_api_poly_coverage.substitute(
            pdb_id=pdb_id, chain_id=chain_id
        )
    )

    polymer_coverage = request_result(
        get_request, f"{pdb_id}_{chain_id}", ignore_error=True
    )

    return polymer_coverage

def fetch_best_structures(
    uniprot_ids: list,
    save_path: str,
    overwrite: bool=False,
) -> dict:
    """ Fetch best structures for given proteins.

    ## Arguments:

    - **uniprot_ids (list)**:<br />
        A list of uniprot ids for which to fetch the best structures.

    - **save_path (str)**:<br />
        The path where the best structures dictionary will be saved as a json file.

    - **overwrite (bool, optional):**:<br />
        Whether to overwrite the existing json file if it already exists. Defaults to False.

    ## Returns:

    - **dict**:<br />
        A dictionary containing the best structures for the given uniprot ids.
    """

    os.makedirs(os.path.dirname(save_path), exist_ok=True)

    if os.path.exists(save_path) and overwrite == False:
        best_structures = read_json(save_path)

    else:
        best_structures = {}
        for uniprot_id in tqdm(uniprot_ids):
            uniprot_base = uniprot_id.split(UNIPROT_ISOFORM_SEPARATOR)[0]
            best_structures[uniprot_id] = get_best_structures(uniprot_base)
        write_json(save_path, best_structures)

    best_structures = add_modeled_coverage(
        best_structures=best_structures,
        outdir=os.path.join(os.path.dirname(save_path), "polymer_coverage")
    )

    return best_structures

def fetch_polymer_coverage(
    pdb_id: str,
    chain_id: str,
    save_path: str,
    overwrite: bool=False,
) -> dict:
    """ Fetch polymer coverage for given pdb id and chain id.

    ## Arguments:

    - **pdb_id (str)**:<br />
        The pdb id for which to fetch the polymer coverage.

    - **chain_id (str)**:<br />
        The chain id for which to fetch the polymer coverage.

    - **save_path (str)**:<br />
        The path where the polymer coverage dictionary will be saved as a json file.

    - **overwrite (bool, optional):**:<br />
        Whether to overwrite the existing json file if it already exists. Defaults to False.

    ## Returns:

    - **dict**:<br />
        A dictionary containing the polymer coverage for the given pdb id and chain id.
    """

    os.makedirs(os.path.dirname(save_path), exist_ok=True)

    if os.path.exists(save_path) and overwrite == False:
        polymer_coverage = read_json(save_path)

    else:
        polymer_coverage = get_polymer_coverage(pdb_id, chain_id)
        write_json(save_path, polymer_coverage)

    return polymer_coverage

def calculate_modeled_coverage(
    pdb_id: str,
    chain_id: str,
    uniprot_length: int,
    outdir: str
):

    polymer_coverage = fetch_polymer_coverage(
        pdb_id=pdb_id,
        chain_id=chain_id,
        save_path=os.path.join(outdir, f"{pdb_id}_{chain_id}_polymer_coverage.json"),
        overwrite=False
    )
    if not isinstance(polymer_coverage, dict):
        return None, None, None
    molecule = polymer_coverage.get(pdb_id, {}).get("molecules", [{}])[0]
    modeled_residues = molecule.get("chains", [{}])[0].get("observed", [])
    modeled_rescount = sum([
        (
            frag["end"]["author_residue_number"] -
            frag["start"]["author_residue_number"] + 1
        )
        for frag in modeled_residues
    ])
    modeled_residues = ",".join([
        (
            f"{frag['start']['author_residue_number']}-"
            f"{frag['end']['author_residue_number']}"
        )
        for frag in modeled_residues
    ]
    )
    coverage = modeled_rescount / uniprot_length

    return coverage, modeled_residues

def add_modeled_coverage(
    best_structures: dict,
    outdir: str,
) -> dict:
    """ Get modeled coverage for best structures.

    ## Arguments:

    - **best_structures (dict)**:<br />
        The best structures dictionary returned by `get_best_structures`.

    ## Returns:

    - **dict**:<br />
        A dictionary containing the best structures with polymer coverage for the given uniprot ids.
    """


    seq_lengths = query_uniprot_api_for_lengths(
        uniprot_ids=list(best_structures.keys()),
    )

    for uniprot_id, best_structure in best_structures.items():

        uniprot_base = uniprot_id.split(UNIPROT_ISOFORM_SEPARATOR)[0]

        if not isinstance(best_structure, dict):
            continue

        best_chains = best_structure[uniprot_base]

        for idx, chain in enumerate(best_chains):

            pdb_id = chain[BestStructureCol.PDB_ID]
            chain_id = chain[BestStructureCol.CHAIN_ID]

            coverage, modeled_residues = calculate_modeled_coverage(
                pdb_id=pdb_id,
                chain_id=chain_id,
                uniprot_length=seq_lengths.get(uniprot_id, 1),
                outdir=outdir
            )
            best_structures[uniprot_id][uniprot_base][idx].update({
                BestStructureCol.MODELED_COVERAGE: f"{coverage:.2f}",
                BestStructureCol.MODELED_RESIDUES: modeled_residues
            })

    return best_structures
