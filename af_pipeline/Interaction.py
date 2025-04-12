import os
from typing import Dict
import warnings
import numpy as np
import pandas as pd
from af_pipeline._Initialize import _Initialize
from utils import get_interaction_map, get_patches_from_matrix, save_map


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
        idr_chains: list = [],
    ):

        super().__init__(
            struct_file_path=struct_file_path,
            data_file_path=data_file_path,
            af_offset=af_offset,
        )

        dir_name = os.path.basename(struct_file_path).split(".")[0]
        output_dir = os.path.join(output_dir, f"{dir_name}_patches")

        self.interaction_map_type = "contact"  # Either contact/distance.
        self.contact_threshold = 8  # Distance threshold in (Angstorm) to define a contact between residue pairs.
        self.plddt_cutoff = 70  # pLDDT cutoff to consider a confident prediction.
        self.idr_plddt_cutoff = 50  # pLDDT cutoff for IDR chains.
        self.pae_cutoff = 5 # PAE cutoff to consider a confident prediction.
        self.output_dir = output_dir
        self.idr_chains = idr_chains # List of chains that are disordered

        self.save_plot = False
        self.save_table = False


    def create_regions_of_interest(self):
        """
        Create regions of interest for all possible chain pairs.

        Returns:
            regions_of_interest (list): list of regions of interest
        """

        regions_of_interest = []
        token_chain_ids = self.token_chain_ids
        chain_pairs = set()

        for chain1 in set(token_chain_ids):
            for chain2 in set(token_chain_ids):
                if chain1 != chain2:
                    pair = tuple(sorted([chain1, chain2]))
                    chain_pairs.add(pair)

        chain_pairs = list(chain_pairs)

        for chain1, chain2 in chain_pairs:

            ch1_start = self.renumber.renumber_chain_res_num(
                chain_res_num=1,
                chain_id=chain1
            )
            ch1_end = self.renumber.renumber_chain_res_num(
                chain_res_num=self.lengths_dict[chain1],
                chain_id=chain1
            )
            ch2_start = self.renumber.renumber_chain_res_num(
                chain_res_num=1,
                chain_id=chain2
            )
            ch2_end = self.renumber.renumber_chain_res_num(
                chain_res_num=self.lengths_dict[chain2],
                chain_id=chain2
            )

            region_of_interest = {
                chain1: [ch1_start, ch1_end],
                chain2: [ch2_start, ch2_end],
            }

            # region_of_interest = self.renumber.renumber_region_of_interest(
            #     region_of_interest=region_of_interest,
            # )

            regions_of_interest.append(region_of_interest)

        return regions_of_interest


    def get_interaction_data(self, region_of_interest: Dict):
        """
        Get the interaction amp, pLDDT, and PAE for the region of interest.

        Args:
            region_of_interest (Dict): Dictionary containing the chain IDs and the residue numbers for the region of interest.

        Returns:
            interaction_map (np.array): binary contact map or distance map
            plddt1 (np.array): plddt values for chain 1
            plddt2 (np.array): plddt values for chain 2
            pae (np.array): PAE matrix for the region of interest
        """

        chain1, chain2 = list(region_of_interest.keys())
        p1_region, p2_region = region_of_interest[chain1], region_of_interest[chain2]

        # chain1 start and end indices.
        start_idx1 = self.num_to_idx[chain1][p1_region[0]]
        end_idx1 = self.num_to_idx[chain1][p1_region[1]]

        # chain2 start and end indices.
        start_idx2 = self.num_to_idx[chain2][p2_region[0]]
        end_idx2 = self.num_to_idx[chain2][p2_region[1]]

        avg_pae = self.avg_pae[start_idx1:end_idx1+1, start_idx2:end_idx2+1]

        coords1 = np.array(self.coords_list[start_idx1:end_idx1+1])
        coords2 = np.array(self.coords_list[start_idx2:end_idx2+1])

        plddt1 = {
            chain1: np.array(self.plddt_list[start_idx1:end_idx1+1])
        }
        plddt2 = {
            chain2: np.array(self.plddt_list[start_idx2:end_idx2+1])
        }

        # Create a contact map or distance map as specified.
        interaction_map = get_interaction_map(
            coords1=coords1.reshape(-1, 3),
            coords2=coords2.reshape(-1, 3),
            contact_threshold=self.contact_threshold,
            map_type=self.interaction_map_type
        )

        return interaction_map, plddt1, plddt2, avg_pae


    def apply_confidence_cutoffs(
        self,
        plddt1: dict,
        plddt2: dict,
        pae: np.array
    ):
        """
        mask low-confidence interactions.

        Args:
            plddt1 (dict): pLDDT values for chain 1
            plddt2 (dict): pLDDT values for chain 2

        Returns:
            plddt_matrix (np.array): binary matrix for plddt values >= plddt_cutoff
            pae (np.array): binary matrix for pae values <= pae_cutoff
        """

        chain1, chain2 = next(iter(plddt1)), next(iter(plddt2))
        plddt1, plddt2 = plddt1[chain1], plddt2[chain2]
        plddt1, plddt2 = plddt1.reshape(-1, 1), plddt2.reshape(-1,1)

        ch1_cutoff = ch2_cutoff = self.plddt_cutoff
        if chain1 in self.idr_chains:
            ch1_cutoff = self.idr_plddt_cutoff
        if chain2 in self.idr_chains:
            ch2_cutoff = self.idr_plddt_cutoff

        plddt1 = np.where(plddt1 >= ch1_cutoff, 1, 0)
        plddt2 = np.where(plddt2 >= ch2_cutoff, 1, 0)
        plddt_matrix = plddt1 * plddt2.T

        pae = np.where(pae <= self.pae_cutoff, 1, 0)

        return plddt_matrix, pae


    def get_confident_interaction_map(self, region_of_interest: Dict):
        """
        For the specified regions in the predicted structure, obtain all confident interacting residue pairs.

        Returns:
            confident_interactions (np.array): binary map of confident interacting residues
        """

        interaction_map, plddt1, plddt2, avg_pae = self.get_interaction_data(
            region_of_interest=region_of_interest
        )

        plddt_matrix, pae_matrix = self.apply_confidence_cutoffs(
            plddt1=plddt1, plddt2=plddt2, pae=avg_pae
        )

        confident_interactions = interaction_map * plddt_matrix * pae_matrix

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
            contact_map (np.array): binary contact map.
            region_of_interest (dict): region of interest for the protein pair.

        Returns:
            patches (dict): interacting patches for the given region of interest of the protein pair.
        """

        patches = {}

        chain1, chain2 = region_of_interest.keys()
        p1_region, p2_region = region_of_interest[chain1], region_of_interest[chain2]

        if np.unique(contact_map).tolist() == [0]: # No interactions found.
            warnings.warn(
                f"No interacting patches found for {chain1}:{p1_region} and {chain2}:{p2_region}."
            )
            return patches

        patches_df = get_patches_from_matrix(
            matrix=contact_map,
            chain1=chain1,
            chain2=chain2
        )

        for patch_idx, patch in patches_df.iterrows():

            ch1_patch = patch[chain1]
            ch2_patch = patch[chain2]

            ch1_patch = sorted([int(x) for x in ch1_patch])
            ch2_patch = sorted([int(x) for x in ch2_patch])

            ch1_patch = np.array(ch1_patch) + region_of_interest[chain1][0]
            ch2_patch = np.array(ch2_patch) + region_of_interest[chain2][0]

            # patches[patch_idx] = {
            #     chain1: (
            #         f"{ch1_patch[0]}-{ch1_patch[-1]}"
            #         if len(ch1_patch) > 1
            #         else str(ch1_patch[0])
            #     ),
            #     chain2: (
            #         f"{ch2_patch[0]}-{ch2_patch[-1]}"
            #         if len(ch2_patch) > 1
            #         else str(ch2_patch[0])
            #     ),
            # }

            patches[patch_idx] = {
                chain1: np.array(ch1_patch),
                chain2: np.array(ch2_patch),
            }

        return patches


    def save_ppair_interaction(
        self,
        region_of_interest: Dict,
        save_plot: bool = False,
        plot_type: str = "static",
        p1_name: str | None = None,
        p2_name: str | None = None,
        concat_residues: bool = True,
        contact_probability: bool = True,
    ):
        """Save the interacting patches for the given region of interest of the protein pair.

        Args:
            region_of_interest (Dict): Dictionary containing the chain IDs and the residue indices for the region of interest.
            save_plot (bool, optional): Outputs the plot if True. Defaults to False.
            plot_type (str, optional): Type of plot to be saved. Defaults to "static"; options: ["static", "interactive", "both"].
            p1_name (str, optional): Name of the first protein. Defaults to None.
            p2_name (str, optional): Name of the second protein. Defaults to None.
            concat_residues (bool, optional): Whether to concatenate the residues into residue ranges. Defaults to True.
            contact_probability (bool, optional): Whether to add contact probability column to the output. Defaults to True.
        """

        chain1, chain2 = list(region_of_interest.keys())
        p1_region, p2_region = region_of_interest[chain1], region_of_interest[chain2]

        contact_map = self.get_confident_interaction_map(
            region_of_interest=region_of_interest
        )

        interacting_patches = self.get_interacting_patches(
            contact_map=contact_map,
            region_of_interest=region_of_interest,
        )

        if p1_name and p2_name:
            p_names = {
                chain1: p1_name,
                chain2: p2_name,
            }
            dir_name_to_replace = "_".join([p1_name, p2_name])
        else:
            p_names = {
                chain1: chain1,
                chain2: chain2,
            }
            dir_name_to_replace = None


        if len(interacting_patches) > 0:

            file_name = "_".join([
                f"{p_names[k]}_{k}:{v[0]}-{v[1]}" for k, v in region_of_interest.items()
            ])

            if dir_name_to_replace:
                dir_name = os.path.basename(self.struct_file_path).split(".")[0]
                # self.output_dir = self.output_dir.replace(dir_name, dir_name_to_replace)

            os.makedirs(self.output_dir, exist_ok=True)

            save_map(
                contact_map=contact_map,
                avg_contact_probs_mat=self.avg_contact_probs_mat,
                patches=interacting_patches,
                chain1=chain1,
                chain2=chain2,
                p1_name=p_names[chain1],
                p2_name=p_names[chain2],
                p1_region=p1_region,
                p2_region=p2_region,
                out_file=os.path.join(self.output_dir, f"patches_{file_name}.html"),
                save_plot=save_plot,
                plot_type=plot_type,
                concat_residues=concat_residues,
                contact_probability=contact_probability,
                num_to_idx=self.num_to_idx,
                idx_to_num=self.idx_to_num,
            )