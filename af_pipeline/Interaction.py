from copy import deepcopy
import os
from typing import Dict
import numpy as np
import pandas as pd
from sklearn.cluster import MeanShift
from af_pipeline._Initialize import _Initialize
# from af_pipeline.archive.af_utils import get_interaction_map #renumber_chain_res_num
from utils import get_key_from_res_range, get_interaction_map, get_patches_from_matrix
# from af_pipeline.archive.af_utils import offset_interacting_region
from utils import save_map

class Interaction(_Initialize):
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
        output_dir: str = "./output/af_output",
    ):

        super().__init__(
            struct_file_path=struct_file_path,
            data_file_path=data_file_path,
            af_offset=af_offset,
        )

        dir_name = os.path.basename(struct_file_path).split(".")[0]
        output_dir = os.path.join(output_dir, f"{dir_name}_patches")
        os.makedirs(output_dir, exist_ok=True)

        # Either contact/distance.
        self.interaction_map_type = "contact"
        # Distance threshold in (Angstorm) to define a contact between residue pairs.
        self.contact_threshold = 8
        # pLDDt cutoff to consider a confident prediction.
        self.plddt_cutoff = 70
        # PAE cutoff to consider a confident prediction.
        self.pae_cutoff = 5
        # self.interacting_region = interacting_region
        self.output_dir = output_dir

        self.save_plot = False
        self.save_table = False

        # if not interacting_regions:
        #     self.interacting_regions = self.create_interacting_regions()

        # else:
        #     self.interacting_regions = interacting_regions
        # for idx, interacting_region in enumerate(self.interacting_regions):
        #     self.save_interaction_info(
        #         interacting_region=interacting_region,
        #     )
        # self.save_interaction_info()

    def save_interaction_info(self, interacting_region, save_plot: bool = False, plot_type: str = "static", save_table: bool = False, interface_only: bool = True, **kwargs):

        # for idx, interacting_region in enumerate(self.interacting_regions):
        # print(interacting_region)
        contact_map = self.get_confident_interactions(
            interacting_region=interacting_region
        )
        # seg_map, interacting_patches = self.get_interacting_patches(
        #     contact_map=contact_map,
        #     interacting_region=interacting_region,
        #     bandwidth=1,
        # )

        interacting_patches = self.get_interacting_patches2(
            contact_map=contact_map,
            interacting_region=interacting_region,
        )

        ir_str = "_".join([f"{k}:{v[0]}-{v[1]}" for k, v in interacting_region.items()])

        save_map(
            # contact_map=seg_map,
            contact_map=contact_map,
            patches=interacting_patches,
            chain1=list(interacting_region.keys())[0],
            chain2=list(interacting_region.keys())[1],
            p1_region=interacting_region[list(interacting_region.keys())[0]],
            p2_region=interacting_region[list(interacting_region.keys())[1]],
            # interacting_region=interacting_region,
            out_file=os.path.join(self.output_dir, f"patches_{ir_str}.html"),
            save_plot=save_plot,
            plot_type=plot_type,
        )

        if save_table:
            # print(kwargs.get("interface_only", True))
            df = self.get_contacts_as_restraints(
                chain_id1=list(interacting_region.keys())[0],
                chain_id2=list(interacting_region.keys())[1],
                p1_region=interacting_region[list(interacting_region.keys())[0]],
                p2_region=interacting_region[list(interacting_region.keys())[1]],
                contact_map=contact_map,
                # interacting_region=interacting_region,
                interface_only=interface_only,
                residue_range=kwargs.get("residue_range", True),
            )

            if interface_only:
                df.to_csv(
                    os.path.join(self.output_dir, f"interface_{ir_str}.csv"),
                    index=False,
                )
            else:
                df.to_csv(
                    os.path.join(self.output_dir, f"pairwise_interactions_{ir_str}.csv"),
                    index=False,
                )

    def create_interacting_regions(self):
        """
        Create interacting regions for all possible chain pairs.

        Returns:
            interacting_regions (list): list of interacting regions
        """

        interacting_regions = []

        chain_pairs = set()

        for chain1 in self.res_dict:
            for chain2 in self.res_dict:
                if chain1 != chain2:
                    pair = tuple(sorted((chain1, chain2)))
                    chain_pairs.add(pair)

        chain_pairs = list(chain_pairs)

        # print(self.af_offset)
        for chain1, chain2 in chain_pairs:
            if self.af_offset:
                ch1_start = self.af_offset[chain1][0] if chain1 in self.af_offset else 1
                ch2_start = self.af_offset[chain2][0] if chain2 in self.af_offset else 1
            else:
                ch1_start = 1
                ch2_start = 1
            ch1_end = ch1_start + len(self.res_dict[chain1]) - 1
            ch2_end = ch2_start + len(self.res_dict[chain2]) - 1
            interacting_regions.append(
                {
                    chain1: (ch1_start, ch1_end),
                    chain2: (ch2_start, ch2_end),
                }
            )

        # for idx, interacting_region in enumerate(interacting_regions):
        #     interacting_regions[idx] = self.renumber.renumber_interacting_region(
        #         interacting_region=interacting_region,
        #     )

        return interacting_regions

    def get_chains_n_indices(self, interacting_region: Dict):
        """
        Obtain the chain IDs and residue indices for the required interacting region. \n
        residue_index = residue_position - 1

        Returns:
            chains (list): list of chain IDs
            mol1_res (list): start and end residue indices for chain 1
            mol2_res (list): start and end residue indices for chain 2
        """

        interacting_region = self.renumber.renumber_interacting_region(
            interacting_region=interacting_region,
        )
        # print(interacting_region)
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

        chains, mol1_res, mol2_res = self.get_chains_n_indices(interacting_region)

        coords1, coords2 = self.get_required_coords(chains, mol1_res, mol2_res)

        # Create a contact map or distance map as specified.
        interaction_map = get_interaction_map(
            coords1, coords2, self.contact_threshold, self.interaction_map_type
        )
        # print(interaction_map.shape)
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

        interaction_map, plddt1, plddt2, pae = self.get_interaction_data(
            interacting_region=interacting_region
        )

        plddt_matrix, pae = self.apply_confidence_cutoffs(
            plddt1=plddt1, plddt2=plddt2, pae=pae
        )
        # print(plddt_matrix.shape, pae.shape, interaction_map.shape)
        confident_interactions = interaction_map * plddt_matrix * pae

        return confident_interactions


    def get_contacts_as_restraints(
        self,
        chain_id1: str,
        chain_id2: str,
        p1_region: tuple,
        p2_region: tuple,
        contact_map: np.array,
        # interacting_region: Dict | None = None,
        interface_only:bool = True,
        residue_range: bool = False,
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

        if interface_only:  # Just write the interface residues.
            res1_idx = sorted(pd.unique(res1_idx))
            res2_idx = sorted(pd.unique(res2_idx))

            start1, start2 = int(p1_region[0]), int(p2_region[0])
            res1_idx = [start1+res_num-1 for res_num in res1_idx]
            res2_idx = [start2+res_num-1 for res_num in res2_idx]


            # if p1_region:
                # res1_idx = [start1+res_num-1 for res_num in res1_idx]
            if residue_range:
                res1_idx = [get_key_from_res_range(res1_idx)]

            if residue_range:
                res2_idx = [get_key_from_res_range(res2_idx)]


            if len(res1_idx) > len(res2_idx):
                diff = len(res1_idx) - len(res2_idx)
                res2_idx = np.append(res2_idx, ["" for x in range(diff)])
            else:
                diff = len(res2_idx) - len(res1_idx)
                res1_idx = np.append(res1_idx, ["" for x in range(diff)])

            print(res1_idx, "res1_idx")
            print(res2_idx, "res2_idx")
            df["Protein1"] = [chain_id1] * len(res1_idx)
            df["Protein2"] = [chain_id2] * len(res2_idx)

            df["Interface_residues1"] = res1_idx
            df["Interface_residues2"] = res2_idx

            # if p1_region:
            #     # res1_idx = [start+res_num-1 for res_num in res1_idx]
            #     df["Residue1"] = get_key_from_res_range(list(res1_idx), as_list=True) if residue_range else res1_idx
            # else:
            #     df["Residue1"] = get_key_from_res_range(list(res1_idx), as_list=True) if residue_range else res1_idx

            # if p2_region:
            #     # start = int(p2_region[0])
            #     # res2_idx = [start+res_num-1 for res_num in res2_idx]
            #     df["Residue2"] = get_key_from_res_range(list(res2_idx), as_list=True) if residue_range else res2_idx
            # else:
            #     df["Residue2"] = get_key_from_res_range(list(res2_idx), as_list=True) if residue_range else res2_idx

        # Write all interacting residue pairs.
        else:
            df["Protein1"] = [chain_id1] * len(res1_idx)
            df["Protein2"] = [chain_id2] * len(res2_idx)

            if p1_region:
                start = int(p1_region[0])
                df["Residue1"] = [start+res_num-1 for res_num in res1_idx]
            else:
                df["Residue1"] = res1_idx

            if p2_region:
                start = int(p2_region[0])
                df["Residue2"] = [start+res_num-1 for res_num in res2_idx]
            else:
                df["Residue2"] = res2_idx

        return df


    def get_interacting_patches(
        self,
        contact_map: np.array,
        interacting_region: dict,
        bandwidth: float = 0, # mod this to 1
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

    def get_interacting_patches2(
        self,
        contact_map: np.array,
        interacting_region: dict,
        # bandwidth: float = 0, # mod this to 1
    ):
        """This is a dirty implementation to get the interacting patches. \n
        This is a temporary solution until we find a better way to get interacting
        patches for the given contact map.

        Args:
            contact_map (np.array): _description_
            interacting_region (dict): _description_

        Returns:
            _type_: _description_
        """

        patches = {}

        chain1, chain2 = interacting_region.keys()

        # interacting_res_pairs = np.argwhere(contact_map == 1)

        # if len(interacting_res_pairs) == 0:  # no interacting residues
        #     return patches

        # for idx, res_idcs in enumerate(interacting_res_pairs):
        #     res1_idx, res2_idx = res_idcs
        #     res1_pos = res1_idx + interacting_region[chain1][0]
        #     res2_pos = res2_idx + interacting_region[chain2][0]
        #     interacting_res_pairs[idx] = np.array([res1_pos, res2_pos])

        # col1_list = []
        # col2_list = []

        # for pairs in interacting_res_pairs:
        #     col1_list.append(str(pairs[0]))
        #     col2_list.append(str(pairs[1]))

        # df0 = pd.DataFrame(columns=[chain1, chain2])
        # df0[chain1] = col1_list
        # df0[chain2] = col2_list

        patches_df = get_patches_from_matrix(contact_map, chain1=chain1, chain2=chain2)

        for patch_idx, patch in patches_df.iterrows():

            ch1_patch = patch[chain1].split("-")
            ch2_patch = patch[chain2].split("-")
            # print(ch1_patch, ch2_patch)
            ch1_patch = [int(x) for x in ch1_patch]
            ch2_patch = [int(x) for x in ch2_patch]
            ch1_patch = np.array(ch1_patch) + interacting_region[chain1][0]
            ch2_patch = np.array(ch2_patch) + interacting_region[chain2][0]
            
            patches[patch_idx] = {
                chain1: f"{ch1_patch[0]}-{ch1_patch[-1]}" if len(ch1_patch) > 1 else str(ch1_patch[0]),
                chain2: f"{ch2_patch[0]}-{ch2_patch[-1]}" if len(ch2_patch) > 1 else str(ch2_patch[0]),
            }

        # for _, patch in patches.items():
        #     for chain in patch.keys():
        #         patch_residues = list(set(patch[chain]))
        #         patch[chain] = get_key_from_res_range(patch_residues)

        return patches