import numpy as np
from af_pipeline._Initialize import _Initialize
from af_pipeline.pae_to_domains.pae_to_domains import (
    parse_pae_file,
    domains_from_pae_matrix_igraph,
    domains_from_pae_matrix_networkx,
)
import os
from collections import defaultdict
from af_pipeline.Parser import ResidueSelect
from utils import get_key_from_res_range, save_pdb, convert_false_to_true, fill_up_the_blanks


class RigidBodies(_Initialize):
    """Class to predict rigid bodies from a PAE file.
    - The rigid bodies (pseudo-domains) are predicted based on the PAE matrix (Graph-based community clustering approach by Tristan Croll).
    - The rigid bodies can be further filtered based on the pLDDT cutoff.
    - The rigid bodies can be saved as PDB files and/or plain txt format specifying the chains and residues in each rigid body.
    """

    def __init__(
        self,
        data_path: str,
        structure_path: str | None = None,
        af_offset: dict | None = None,
    ):

        super().__init__(
            data_file_path=data_path,
            struct_file_path=structure_path,
            af_offset=af_offset,
        )

        self.library = "igraph"
        self.pae_power = 1
        self.pae_cutoff = 5
        self.resolution = 0.5
        self.plddt_cutoff = 70
        self.patch_threshold = 10


    def predict_domains(
        self,
        num_res: int = 5,
        num_proteins: int = 1,
        plddt_filter: bool = True
    ):
        """Predict domains from a PAE file.
        - Two implementations are available:
            1. igraph based
            2. networkx based

        (1) is significantly faster than (2)

        Args:
            num_res (int): Minimum number of residues in a rigid body
            num_proteins (int): Minimum number of proteins in a rigid body
            plddt_filter (bool): Filter the residues based on the pLDDT cutoff

        Raises:
            ValueError: Invalid library specified. Use 'igraph' or 'networkx'

        Returns:
            domains (list): List of domains in which each domain is a rigid body dictionary

        A rigid body dictionary is of the form:
        - {
            chain_id1: [res_num, ...],
            chain_id2: [res_num, ...],
            ...
        }
        """

        pae_path = self.data_file_path
        pae_matrix = parse_pae_file(pae_file=pae_path)

        if self.library == "igraph":
            f = domains_from_pae_matrix_igraph

        elif self.library == "networkx":
            f = domains_from_pae_matrix_networkx

        else:
            raise ValueError("Invalid library specified. Use 'igraph' or 'networkx")

        domains = f(
            pae_matrix,
            pae_power=self.pae_power,
            pae_cutoff=self.pae_cutoff,
            graph_resolution=self.resolution,
        )

        for idx, domain in enumerate(domains):

            if isinstance(domain, frozenset):
                domain = list(domain)

            rb_dict = self.domain_to_rb_dict(domain=domain)

            if plddt_filter:
                rb_dict = self.filter_plddt(
                    rb_dict=rb_dict,
                    patch_threshold=self.patch_threshold,
                )

            domains[idx] = rb_dict

        # Filter out domains with less than num_proteins proteins
        domains = [
            rb_dict
            for rb_dict in domains
            if len(rb_dict) >= num_proteins
        ]

        # Filter out domains with less than num_res residues
        domains = [
            rb_dict
            for rb_dict in domains
            if sum([len(res_list) for res_list in rb_dict.values()]) >= num_res
        ]

        return domains


    def domain_to_rb_dict(self, domain: list):
        """Convert the domain list to a dictionary of rigid bodies.
        - The rigid bodies are represented as a dictionary with chain_id as the key and
            a list of residue numbers as the value.

        Args:
            domain (list): list of residue indices in the domain

        Returns:
            rb_dict (dict): pseudo-rigid body in the form of a dictionary

        Example:
            if predicted structure has chains: A (20 aa), B (30 aa), C (50 aa) \n
            and detected domain is [0, 1, 2, 3, 4, 5, 20, 21, 22, 23, 54, 55, 56, 57, 58] \n
            rb_dict = {
                'A': [1, 2, 3, 4, 5],
                'B': [1, 2, 3, 4],
                'C': [5, 6, 7, 8, 9]
            }
        """

        rb_dict = defaultdict(list)

        for res_idx in domain:

            res_num = self.idx_to_num[res_idx].get("res_num")
            chain_id = self.idx_to_num[res_idx].get("chain_id")

            rb_dict[chain_id].append(res_num)

        return rb_dict


    def filter_plddt(
        self,
        rb_dict: dict,
        patch_threshold: int = 5,
    ):
        """Filter the residues in the rigid bodies based on the pLDDT cutoff.
        - If the pLDDT score of a residue is less than the cutoff, it is removed from the rigid body.
        Args:
            rb_dict (dict): dictionary of rigid bodies
            patch_threshold (int): minimum number of contiguous residues for which the pLDDT score is above the cutoff

        Returns:
            rb_dict (dict): dictionary of rigid bodies with residues filtered based on the pLDDT cutoff
        """

        # Filter the residues in each chain in the rigid body based on the pLDDT cutoff
        for chain_id, rb_res_num_list in rb_dict.items():

            confident_residues = []

            # sorted list of residue numbers in the rigid body
            rb_res_num_arr = np.array(sorted(rb_res_num_list))

            # sorted list of residue indices in the rigid body
            plddt_res_num_arr = np.array([self.num_to_idx[chain_id][res_num] for res_num in rb_res_num_list])

            # True/False array based on the pLDDT cutoff
            # for e.g. plddt_arr = [70, 78, 90, 65, 65, 80, 90]
            # tf_plddt_filtered = [True, True, True, False, False, True, True] for cutoff = 70
            tf_plddt_filtered = np.array(self.plddt_list)[plddt_res_num_arr] >= self.plddt_cutoff

            # Convert the pLDDT scores to True/False based on the threshold
            # for e.g. if arr = [True, False, False, False, True, True] and threshold = 3
            # the output will be [True, True, True, True, True, True]
            tf_plddt_filtered = convert_false_to_true(
                arr=tf_plddt_filtered,
                threshold=patch_threshold,
            )

            # Get the residue numbers of the confident residues
            confident_residues = rb_res_num_arr[tf_plddt_filtered]
            rb_dict[chain_id] = confident_residues.tolist()

        # Remove chains which have no confident residues
        empty_chains = []

        for chain_id, confident_residues in rb_dict.items():
            if len(confident_residues) == 0:
                empty_chains.append(chain_id)

        for chain_id in empty_chains:
            del rb_dict[chain_id]

        return rb_dict


    def save_rigid_bodies(
        self,
        domains: list,
        output_dir: str,
        output_format: str = "txt",
        save_structure: bool = True,
        no_plddt_filter_for_structure: bool = False,
    ):
        """Save the rigid bodies to a text file."""

        output_dir = os.path.join(output_dir)
        dir_name = os.path.basename(self.struct_file_path).split(".")[0]
        output_dir = os.path.join(output_dir, dir_name)

        os.makedirs(output_dir, exist_ok=True)

        file_name = (
            os.path.basename(self.struct_file_path).split(".")[0] + "_rigid_bodies"
        )

        if output_format == "txt":
            file_name += ".txt"
            output_path = os.path.join(output_dir, file_name)

            with open(output_path, "w") as f:

                for idx, rb_dict in enumerate(domains):
                    f.write(f"Rigid Body {idx}\n")

                    for chain_id, res_list in rb_dict.items():

                        if len(res_list) > 0:
                            f.write(
                                f"{chain_id}:{get_key_from_res_range(res_range=res_list)}\n"
                            )

                    f.write("\n")

        if save_structure:

            structure = self.renumber.renumber_structure(
                structure=self.structureparser.structure,
            )

            for idx, rb_dict in enumerate(domains):

                if no_plddt_filter_for_structure:
                    for chain_id, res_list in rb_dict.items():
                        if len(res_list) > 0:
                            res_list = fill_up_the_blanks(res_list)
                            rb_dict[chain_id] = res_list

                output_path = os.path.join(output_dir, f"rigid_body_{idx}.cif")

                save_pdb(
                    structure=structure,
                    out_file=output_path,
                    res_select_obj=ResidueSelect(rb_dict),
                    save_type="cif",
                    preserve_header_footer=False,
                )