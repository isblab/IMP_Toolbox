# Description: Fetch best structures for given proteins

import sys
sys.path.append("../")
from utils import request_session, request_result, read_json, write_json
import pandas as pd
import os
from argparse import ArgumentParser
from tqdm import tqdm


def get_best_structures(uniprot_id):
    """Get best structures for given uniprot id

    Args:
        uniprot_id (str): uniprot id

    Returns:
        dict: best structures for given uniprot id
    """

    pisa_url = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures"

    req_sess = request_session(max_retries=3)
    get_request = req_sess.get(
        url=f"{pisa_url}/{uniprot_id}"
    )

    best_structures = request_result(
        get_request, uniprot_id, ignore_error=True
    )

    return best_structures


def make_best_structures_df(best_structures):
    """Create a dataframe of best structures from dictionary

    Args:
        best_structures (dict): best structures for given uniprot ids

    Returns:
        pd.DataFrame: dataframe of best structures
    """

    all_best_chains = []

    for uniprot_id, best_structure in best_structures.items():

        if best_structure:
            best_chains = best_structure[uniprot_id]

            for _, chain in enumerate(best_chains):

                all_best_chains.append(
                {
                    "uniprot_id": uniprot_id,
                    "best_chain": chain["chain_id"],
                    "best_structure": chain["pdb_id"],
                    "coverage": chain["coverage"],
                    "unp_start": chain["unp_start"],
                    "unp_end": chain["unp_end"],
                    "start": chain["start"],
                    "end": chain["end"],
                    "resolution": chain["resolution"],
                }
            )

        else:
            print(f"No best structure found for {uniprot_id}")

    df = pd.DataFrame(all_best_chains)
    mask = df.duplicated(subset=['uniprot_id'], keep='first')
    df['uniprot_id'] = df['uniprot_id'].mask(mask, '')
    df = df.drop_duplicates(subset=['best_structure'], keep='first')

    return df


def fetch_best_structures(proteins_dict, save_path, overwrite=False):
    """Fetch best structures for given proteins

    Args:
        proteins_dict (dict): dictionary of proteins
        overwrite (bool, optional): Defaults to False.

    Returns:
        dict: best structures for given proteins
    """

    if os.path.exists(save_path) and overwrite == False:
        best_structures = read_json(save_path)

    else:
        best_structures = {}
        for uniprot_id in tqdm(proteins_dict.values()):
            uniprot_id = uniprot_id.split("-")[0]
            best_structures[uniprot_id] = get_best_structures(uniprot_id)
        write_json(save_path, best_structures)

    return best_structures


if __name__ == "__main__":
    args = ArgumentParser()
    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="../inputs/cardiac_desmosome_proteins.json",
    )
    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="../output/best_structures.csv",
    )
    args.add_argument(
        "-ow",
        "--overwrite",
        action="store_true",
        required=False,
        default=False,
    )

    args = args.parse_args()
    proteins_dict = read_json(args.input)
    best_structures = fetch_best_structures(
        proteins_dict,
        save_path="../output/best_structures.json",
        overwrite=args.overwrite
    )
    df = make_best_structures_df(best_structures)
    df.to_csv(args.output, index=False)