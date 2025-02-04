from copy import deepcopy
import os
from typing import Dict
import numpy as np
from af_pipeline._Initialize import _Initialize
from utils import get_interaction_map, get_patches_from_matrix
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

        self.interaction_map_type = "contact"  # Either contact/distance.
        self.contact_threshold = 8  # Distance threshold in (Angstorm) to define a contact between residue pairs.
        self.plddt_cutoff = 70  # pLDDT cutoff to consider a confident prediction.
        self.pae_cutoff = 5 # PAE cutoff to consider a confident prediction.
        self.output_dir = output_dir

        self.save_plot = False
        self.save_table = False
        # print(self.num_to_idx)

    def save_ppair_interaction(
        self,
        region_of_interest: Dict,
        save_plot: bool = False,
        plot_type: str = "static",
    ):
        """Save the interacting patches for the given region of interest of the protein pair.

        Args:
            region_of_interest (Dict): Dictionary containing the chain IDs and the residue indices for the region of interest.
            save_plot (bool, optional): Outputs the plot if True. Defaults to False.
            plot_type (str, optional): Type of plot to be saved. Defaults to "static"; options: ["static", "interactive", "both"].
        """

        chain1, chain2 = list(region_of_interest.keys())
        p1_region, p2_region = region_of_interest[chain1], region_of_interest[chain2]

        contact_map = self.get_confident_interaction_map(
            region_of_interest=region_of_interest
        )

        if np.unique(contact_map).size == 1:
            print(f"No confident interactions found for {chain1}:{p1_region} and {chain2}:{p2_region}.")

        else:
            interacting_patches = self.get_interacting_patches(
                contact_map=contact_map,
                region_of_interest=region_of_interest,
            )

            file_name = "_".join([
                f"{k}:{v[0]}-{v[1]}" for k, v in region_of_interest.items()
            ])

            save_map(
                contact_map=contact_map,
                patches=interacting_patches,
                chain1=chain1,
                chain2=chain2,
                p1_region=p1_region,
                p2_region=p2_region,
                out_file=os.path.join(self.output_dir, f"patches_{file_name}.html"),
                save_plot=save_plot,
                plot_type=plot_type,
            )


    def create_regions_of_interest(self):
        """
        Create regions of interest for all possible chain pairs.

        Returns:
            regions_of_interest (list): list of regions of interest
        """

        regions_of_interest = []

        chain_pairs = set()

        for chain1 in self.res_dict:
            for chain2 in self.res_dict:
                if chain1 != chain2:
                    pair = tuple(sorted((chain1, chain2)))
                    chain_pairs.add(pair)

        chain_pairs = list(chain_pairs)

        for chain1, chain2 in chain_pairs:

            if self.af_offset:
                ch1_start = self.af_offset[chain1][0] if chain1 in self.af_offset else 1
                ch2_start = self.af_offset[chain2][0] if chain2 in self.af_offset else 1
            else:
                ch1_start = 1
                ch2_start = 1

            ch1_end = ch1_start + len(self.res_dict[chain1]) - 1
            ch2_end = ch2_start + len(self.res_dict[chain2]) - 1

            regions_of_interest.append(
                {
                    chain1: (ch1_start, ch1_end),
                    chain2: (ch2_start, ch2_end),
                }
            )

        return regions_of_interest


    def get_interaction_data(self, region_of_interest: Dict):
        """
        Get the interaction amp, pLDDT, and PAE for the region of interest.

        Returns:
            interaction_map (np.array): binary contact map or distance map
            plddt1 (np.array): plddt values for chain 1
            plddt2 (np.array): plddt values for chain 2
            pae (np.array): PAE matrix for the region of interest
        """

        chain1, chain2 = list(region_of_interest.keys())

        # chain1 start and end indices.
        start_idx1 = self.num_to_idx[chain1][region_of_interest[chain1][0]]
        end_idx1 = self.num_to_idx[chain1][region_of_interest[chain1][1]]

        # chain2 start and end indices.
        start_idx2 = self.num_to_idx[chain2][region_of_interest[chain2][0]]
        end_idx2 = self.num_to_idx[chain2][region_of_interest[chain2][1]]

        pae = self.avg_pae[start_idx1:end_idx1+1, start_idx2:end_idx2+1]

        region_of_interest = self.renumber.renumber_region_of_interest(
            region_of_interest=region_of_interest,
        )

        mol1_res1, mol1_res2 = region_of_interest[chain1]
        mol2_res1, mol2_res2 = region_of_interest[chain2]

        coords1 = self.coords_dict[chain1][mol1_res1-1:mol1_res2, :]
        coords2 = self.coords_dict[chain2][mol2_res1-1:mol2_res2, :]

        plddt1 = self.plddt_dict[chain1][mol1_res1-1:mol1_res2]
        plddt2 = self.plddt_dict[chain2][mol2_res1-1:mol2_res2]

        # Create a contact map or distance map as specified.
        interaction_map = get_interaction_map(
            coords1=coords1,
            coords2=coords2,
            contact_threshold=self.contact_threshold,
            map_type=self.interaction_map_type
        )

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


    def get_confident_interaction_map(self, region_of_interest: Dict):
        """
        For the specified regions in the predicted structure, obtain all confident interacting residue pairs.

        Returns:
            confident_interactions (np.array): binary map of confident interacting residues
        """

        interaction_map, plddt1, plddt2, pae = self.get_interaction_data(
            region_of_interest=region_of_interest
        )

        plddt_matrix, pae = self.apply_confidence_cutoffs(
            plddt1=plddt1, plddt2=plddt2, pae=pae
        )

        confident_interactions = interaction_map * plddt_matrix * pae

        return confident_interactions


    def get_interacting_patches(
        self,
        contact_map: np.array,
        region_of_interest: dict,
    ):
        """This is a dirty implementation to get the interacting patches. \n
        This is a temporary solution until we find a better way to get interacting
        patches for the given contact map.

        Args:
            contact_map (np.array): _description_
            region_of_interest (dict): _description_

        Returns:
            _type_: _description_
        """

        patches = {}

        chain1, chain2 = region_of_interest.keys()

        patches_df = get_patches_from_matrix(contact_map, chain1=chain1, chain2=chain2)

        for patch_idx, patch in patches_df.iterrows():

            ch1_patch = patch[chain1]
            ch2_patch = patch[chain2]

            ch1_patch = sorted([int(x) for x in ch1_patch])
            ch2_patch = sorted([int(x) for x in ch2_patch])

            ch1_patch = np.array(ch1_patch) + region_of_interest[chain1][0]
            ch2_patch = np.array(ch2_patch) + region_of_interest[chain2][0]

            patches[patch_idx] = {
                chain1: f"{ch1_patch[0]}-{ch1_patch[-1]}" if len(ch1_patch) > 1 else str(ch1_patch[0]),
                chain2: f"{ch2_patch[0]}-{ch2_patch[-1]}" if len(ch2_patch) > 1 else str(ch2_patch[0]),
            }

        return patches