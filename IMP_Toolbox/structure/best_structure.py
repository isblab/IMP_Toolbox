import os
import pandas as pd
from tqdm import tqdm
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

def get_best_structures(uniprot_id: str) -> dict:
    """ Get the best structures for a given uniprot id.

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

def make_best_structures_df(
    best_structures: dict,
    uniprot_protein_map: dict = {},
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

        if not best_structure:
            print(f"No best structure found for {uniprot_id}")

        best_chains = best_structure[uniprot_id.split("-")[0]]
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
                BestStructureCol.PROTEIN: uniprot_protein_map.get(uniprot_id, ""),
                BestStructureCol.CHAIN_ID: chain[BestStructureCol.CHAIN_ID],
                BestStructureCol.PDB_ID: chain[BestStructureCol.PDB_ID],
                BestStructureCol.COVERAGE: chain[BestStructureCol.COVERAGE],
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
            # uniprot_id = uniprot_id.split("-")[0]
            best_structures[uniprot_id] = get_best_structures(
                uniprot_id.split(UNIPROT_ISOFORM_SEPARATOR)[0]
            )
        write_json(save_path, best_structures)

    return best_structures