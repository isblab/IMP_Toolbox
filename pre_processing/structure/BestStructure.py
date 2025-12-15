from utils import request_session, request_result, read_json, write_json
import pandas as pd
import os
from tqdm import tqdm

class BestStructures:
    """Fetch best structures for given proteins
    """

    def __init__(self, uniprot_ids):
        self.uniprot_ids = uniprot_ids

    def get_best_structures(self, uniprot_id: str) -> dict:
        """Get best structures for given uniprot id

        Args:
            uniprot_id (str): uniprot id

        Returns:
            dict: best structures for given uniprot id
        """

        PDBE_API_URL = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures"

        req_sess = request_session(max_retries=3)
        get_request = req_sess.get(
            url=f"{PDBE_API_URL}/{uniprot_id}"
        )

        best_structures = request_result(
            get_request, uniprot_id, ignore_error=True
        )

        return best_structures


    def make_best_structures_df(self, best_structures: dict) -> pd.DataFrame:
        """Create a dataframe of best structures from dictionary

        Args:
            best_structures (dict): best structures for given uniprot ids

        Returns:
            pd.DataFrame: dataframe of best structures
        """

        all_best_chains = []

        for uniprot_id, best_structure in best_structures.items():

            if not best_structure:
                print(f"No best structure found for {uniprot_id}")

            best_chains = best_structure[uniprot_id]
            best_chains = sorted(best_chains, key=lambda x: x['resolution'])

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

        df = pd.DataFrame(all_best_chains)
        mask = df.duplicated(subset=['uniprot_id'], keep='first')
        df['uniprot_id'] = df['uniprot_id'].mask(mask, '')
        df = df.drop_duplicates(subset=['best_structure'], keep='first')

        return df


    def fetch_best_structures(
        self,
        save_path: str,
        overwrite: bool=False,
    ) -> dict:
        """Fetch best structures for given proteins

        Args:
            overwrite (bool, optional): Defaults to False.

        Returns:
            dict: best structures for given proteins
        """

        os.makedirs(os.path.dirname(save_path), exist_ok=True)

        if os.path.exists(save_path) and overwrite == False:
            best_structures = read_json(save_path)

        else:
            best_structures = {}
            for uniprot_id in tqdm(self.uniprot_ids):
                uniprot_id = uniprot_id.split("-")[0]
                best_structures[uniprot_id] = self.get_best_structures(uniprot_id)
            write_json(save_path, best_structures)

        return best_structures