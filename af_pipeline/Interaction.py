from typing import Dict
import numpy as np
import pandas as pd
from sklearn.cluster import MeanShift
from af_pipeline.Initialize import Initialize
from af_pipeline.af_utils import get_interaction_map
from utils import get_key_from_res_range
from af_pipeline.af_utils import offset_interacting_region


class Interaction(Initialize):
    """ Class to handle interaction data for the predicted structure. \n
    One can obtain:
    1. Interaction map: A binary contact map or distance map.
    2. Restraints: Contacts as pairwise residues in a DataFrame format.
    3. Interacting patches: contiguous regions in the interaction map obtained in (1).
    """

    def __init__(
        self,
        struct_file_path: str,
        data_file_path: str,
        af_offset: dict | None = None,
    ):

        super().__init__(
            struct_file_path=struct_file_path,
            data_file_path=data_file_path,
            af_offset=af_offset,
        )

        # Either contact/distance.
        self.interaction_map_type = "contact"
        # Distance threshold in (Angstorm) to define a contact between residue pairs.
        self.contact_threshold = 8
        # pLDDt cutoff to consider a confident prediction.
        self.plddt_cutoff = 70
        # PAE cutoff to consider a confident prediction.
        self.pae_cutoff = 5


    def get_chains_n_indices(self, interacting_region: Dict):
        """
        Obtain the chain IDs and residue indices for the required interacting region. \n
        residue_index = residue_position - 1

        Returns:
            chains (list): list of chain IDs
            mol1_res (list): start and end residue indices for chain 1
            mol2_res (list): start and end residue indices for chain 2
        """

        interacting_region = offset_interacting_region(interacting_region, self.af_offset)
        chain1, chain2 = interacting_region.keys()
        mol1_res1, mol1_res2 = interacting_region[chain1]
        mol1_res1 -= 1
        mol2_res1, mol2_res2 = interacting_region[chain2]
        mol2_res1 -= 1

        return [chain1, chain2], [mol1_res1, mol1_res2], [mol2_res1, mol2_res2]


    def get_required_coords(self, chains: list, mol1_res: list, mol2_res: list):
        """
        Get the coordinates for the interacting region for which confident interactions are required.

        Returns:
            coords1 (np.array): coordinates for chain 1
            coords2 (np.array): coordinates for chain 2
        """

        chain1, chain2 = chains
        start1, end1 = mol1_res
        start2, end2 = mol2_res
        coords1 = self.coords_dict[chain1][start1:end1, :]
        coords2 = self.coords_dict[chain2][start2:end2, :]

        return coords1, coords2


    def get_required_plddt(self, chains: list, mol1_res: list, mol2_res: list):
        """
        Get the plddt for the interacting region for which confident interactions are required.

        Returns:
            plddt1 (np.array): plddt values for chain 1
            plddt2 (np.array): plddt values for chain 2
        """

        chain1, chain2 = chains
        start1, end1 = mol1_res
        start2, end2 = mol2_res
        plddt1 = self.plddt_dict[chain1][start1:end1]
        plddt2 = self.plddt_dict[chain2][start2:end2]

        return plddt1, plddt2


    def get_required_pae(self, chains: list, mol1_res: list, mol2_res: list):
        """
        Get the PAE matrix for the interacting region. \n
        For this we need the cumulative residue index uptil the required residue position.

        Returns:
            pae (np.array): PAE matrix for the interacting region
        """

        chain1, chain2 = chains
        start1, end1 = mol1_res
        start2, end2 = mol2_res

        # Count total residues till start1 and start2.
        cum_start1, cum_start2 = 0, 0
        for chain in self.res_dict:
            if chain == chain1:
                cum_start1 += start1
                break
            else:
                cum_start1 += len(self.res_dict[chain])

        for chain in self.res_dict:
            if chain == chain2:
                cum_start2 += start2
                break
            else:
                cum_start2 += len(self.res_dict[chain])

        cum_end1 = cum_start1 + (end1 - start1)
        cum_end2 = cum_start2 + (end2 - start2)

        pae = self.pae[cum_start1:cum_end1, cum_start2:cum_end2]

        return pae


    def get_interaction_data(self, interacting_region: Dict):
        """
        Get the interaction amp, pLDDT, and PAE for the interacting region.

        Returns:
            interaction_map (np.array): binary contact map or distance map
            plddt1 (np.array): plddt values for chain 1
            plddt2 (np.array): plddt values for chain 2
            pae (np.array): PAE matrix for the interacting region
        """

        interacting_region = offset_interacting_region(interacting_region, self.af_offset)

        chains, mol1_res, mol2_res = self.get_chains_n_indices(interacting_region)

        coords1, coords2 = self.get_required_coords(chains, mol1_res, mol2_res)

        # Create a contact map or distance map as specified.
        interaction_map = get_interaction_map(
            coords1, coords2, self.contact_threshold, self.interaction_map_type
        )

        plddt1, plddt2 = self.get_required_plddt(chains, mol1_res, mol2_res)

        pae = self.get_required_pae(chains, mol1_res, mol2_res)

        return interaction_map, plddt1, plddt2, pae


    def apply_confidence_cutoffs(self, plddt1: np.array, plddt2: np.array, pae: np.array):
        """
        mask low-confidence interactions.

        Returns:
            plddt_matrix (np.array): binary matrix for plddt values >= plddt_cutoff
            pae (np.array): binary matrix for pae values <= pae_cutoff
        """

        plddt1 = np.where(plddt1 >= self.plddt_cutoff, 1, 0)
        plddt2 = np.where(plddt2 >= self.plddt_cutoff, 1, 0)
        plddt_matrix = plddt1 * plddt2.T

        pae = np.where(pae <= self.pae_cutoff, 1, 0)

        return plddt_matrix, pae


    def get_confident_interactions(self, interacting_region: Dict):
        """
        For the specified regions in the predicted structure, obtain all confident interacting residue pairs.

        Returns:
            confident_interactions (np.array): binary map of confident interacting residues
        """

        interacting_region = offset_interacting_region(
            interacting_region=interacting_region,
            af_offset=self.af_offset
        )

        interaction_map, plddt1, plddt2, pae = self.get_interaction_data(
            interacting_region=interacting_region
        )

        plddt_matrix, pae = self.apply_confidence_cutoffs(
            plddt1=plddt1, plddt2=plddt2, pae=pae
        )

        confident_interactions = interaction_map * plddt_matrix * pae

        return confident_interactions


    def get_contacts_as_restraints(
        self,
        prot1_name: str,
        prot2_name: str,
        contact_map: np.array,
        interface_only=True,
    ):
        """
        Given a contact map, convert it the IMP compatible XL-restraint data format.

        Returns:
            df (pd.DataFrame): DataFrame with interacting residue
        """

        idx = np.where(contact_map != 0)

        res1_idx = idx[0]
        res2_idx = idx[1]

        # Index to residue no.
        res1_idx += 1
        res2_idx += 1

        df = pd.DataFrame()
        # Just write the interface residues.
        if interface_only:
            res1_idx = sorted(pd.unique(res1_idx))
            res2_idx = sorted(pd.unique(res2_idx))

            if len(res1_idx) > len(res2_idx):
                diff = len(res1_idx) - len(res2_idx)
                res2_idx = np.append(res2_idx, ["" for x in range(diff)])
            else:
                diff = len(res2_idx) - len(res1_idx)
                res1_idx = np.append(res1_idx, ["" for x in range(diff)])

            df["Protein1"] = [prot1_name] * len(res1_idx)
            df["Residue1"] = res1_idx
            df["Protein2"] = [prot2_name] * len(res2_idx)
            df["Residue2"] = res2_idx

        # Write all interacting residue pairs.
        else:
            df["Protein1"] = [prot1_name] * len(res1_idx)
            df["Residue1"] = res1_idx
            df["Protein2"] = [prot2_name] * len(res2_idx)
            df["Residue2"] = res2_idx

        return df


    def get_interacting_patches(
        self,
        contact_map: np.array,
        interacting_region: dict,
        bandwidth: float | None = None,
    ):
        """Get the interacting patches and the corrresponding segmented map

        Args:
            contact_map (np.array): a binary contact map
            interacting_region (dict): interacting region in the format (after af_offset) \n
                                        {chain1: (start1, end1), chain2: (start2, end2)}
            bandwidth (float, optional): bandwidth for the MeanShift clustering. Defaults to None. \n
                                        If not given, the bandwidth is estimated using sklearn.cluster.estimate_bandwidth

        Returns:
            segmented_map (np.array): segmented map of interacting patches
            patches (dict): interacting patches in the format \n
                            {patch_id: {chain_id1: "1-50,375", chain_id2: "2-10"}}
        """

        patches = {}
        interacting_res_pairs = np.argwhere(contact_map == 1)

        if len(interacting_res_pairs) == 0:  # no interacting residues
            return contact_map, patches

        chain1, chain2 = interacting_region.keys()

        ms = MeanShift(bandwidth=bandwidth)
        ms.fit(interacting_res_pairs)
        labels = ms.labels_
        # cluster_centers = ms.cluster_centers_

        for res_idcs, label in zip(interacting_res_pairs, labels):
            res1_idx, res2_idx = res_idcs
            res1_pos = res1_idx + interacting_region[chain1][0]
            res2_pos = res2_idx + interacting_region[chain2][0]

            if label + 1 not in patches:
                patches[int(label + 1)] = {
                    chain1: [res1_pos],
                    chain2: [res2_pos],
                }

            else:
                patches[int(label + 1)][chain1].append(res1_pos)
                patches[int(label + 1)][chain2].append(res2_pos)

        for label, patch in patches.items():
            for chain in patch.keys():
                patch_residues = list(set(patch[chain]))
                patch[chain] = get_key_from_res_range(patch_residues)

        # segmented map of interacting patches
        segmented_map = np.zeros_like(contact_map)
        for i, res_idx in enumerate(interacting_res_pairs):
            segmented_map[res_idx[0], res_idx[1]] = labels[i] + 1

        return segmented_map, patches
