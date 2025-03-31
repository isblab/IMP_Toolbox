from pprint import pprint
import numpy as np
from af_pipeline._Initialize import _Initialize
from af_pipeline.Interaction import Interaction
from af_pipeline.pae_to_domains.pae_to_domains import (
    parse_pae_file,
    domains_from_pae_matrix_igraph,
    domains_from_pae_matrix_networkx,
    domains_from_pae_matrix_label_propagation
)
import os
from collections import defaultdict
from af_pipeline.Parser import ResidueSelect
from utils import get_key_from_res_range, save_structure_obj, convert_false_to_true, fill_up_the_blanks, get_interaction_map
from itertools import combinations


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
        idr_chains: list = [],
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
        self.plddt_cutoff_idr = 50
        self.patch_threshold = 0
        self.random_seed = 99
        self.idr_chains = idr_chains


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
            ch1: [res_num, ...],
            ch2: [res_num, ...],
            ...
        }
        """

        pae_path = self.data_file_path
        pae_matrix = parse_pae_file(pae_file=pae_path)

        if self.library == "igraph":
            f = domains_from_pae_matrix_igraph

        elif self.library == "networkx":
            f = domains_from_pae_matrix_networkx

        elif self.library == "label_propagation":
            f = domains_from_pae_matrix_label_propagation

        else:
            raise ValueError("Invalid library specified. Use 'igraph' or 'networkx")

        if f == domains_from_pae_matrix_igraph or f == domains_from_pae_matrix_networkx:
            domains = f(
                pae_matrix,
                pae_power=self.pae_power,
                pae_cutoff=self.pae_cutoff,
                graph_resolution=self.resolution,
            )
        elif f == domains_from_pae_matrix_label_propagation:
            domains = f(
                pae_matrix,
                pae_power=self.pae_power,
                pae_cutoff=self.pae_cutoff,
                # random_seed=99,
                # random_seed=9,
                random_seed=self.random_seed,
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
            # TODO shouldnt we keep only residues that pass plddt ? i.e. patch_threshold = 0 not 10, allowing IMP to fill in missing regions.
            # Yes, it is  set to 0, but I was setting it to 0 while using the RigidBodies class. I've changed default to 0 in the function.

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
            such that, actual residue numbers are
            A: [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39] \n
            B: [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49] \n
            C: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49] \n
            and detected domain is [0, 1, 2, 3, 4, 5, 20, 21, 22, 23, 54, 55, 56, 57, 58] \n
            rb_dict = {
                'A': [20, 21, 22, 23, 24, 25],
                'B': [20, 21, 22, 23, 24],
                'C': [5, 6, 7, 8, 9]
            }
        """ #TODO indices 0, 20 are not represented in output in example.
        # TODO shouldnt it be 'A': [1,2,3,4,5,6] for instance?
        # The function is correct I think. you can take a look at the RenumberResidues.residue_map function in Parser.py. Copilot generated incorrect docstring for the example I guess. Fixed it now.

        rb_dict = defaultdict(list)

        for res_idx in domain:

            res_num = self.idx_to_num[res_idx].get("res_num")
            chain_id = self.idx_to_num[res_idx].get("chain_id")

            rb_dict[chain_id].append(res_num)

        return rb_dict


    def filter_plddt(
        self,
        rb_dict: dict,
        patch_threshold: int = 0,
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
            if chain_id in self.idr_chains:
                tf_plddt_filtered = np.array(self.plddt_list)[plddt_res_num_arr] >= self.plddt_cutoff_idr
            else:
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

                save_structure_obj(
                    structure=structure,
                    out_file=output_path,
                    res_select_obj=ResidueSelect(rb_dict),
                    save_type="cif",
                    preserve_header_footer=False,
                )
                # TODO rename function name to save_structure or something since we are using CIF not PDB
                # renamed to save_structure_obj since save_structure is used as a flag in save_rigid_bodies function


    def chain_pair_condition(self, chain_pair, interface_type):
        """Check the interface type based on the chain pair.

        Args:
            chain_pair (tuple): Pair of chain IDs
            interface_type (str): IDR-R, R-R, IDR-IDR, any-any, IDR-any, R-any

        Returns:
            bool: True if the chain pair satisfies the interface type condition, False otherwise
        """

        condition_idr_r = (
            (chain_pair[0] in self.idr_chains and chain_pair[1] not in self.idr_chains)
            or (chain_pair[1] in self.idr_chains and chain_pair[0] not in self.idr_chains)
        )

        condition_r_r = (
            chain_pair[0] not in self.idr_chains and chain_pair[1] not in self.idr_chains
        )

        condition_idr_idr = (
            chain_pair[0] in self.idr_chains and chain_pair[1] in self.idr_chains
        )

        condition_idr_any = (
            chain_pair[0] in self.idr_chains or chain_pair[1] in self.idr_chains
        )

        condition_r_any = (
            chain_pair[0] not in self.idr_chains or chain_pair[1] not in self.idr_chains
        )

        condition_any_any = True

        return locals()[f"condition_{interface_type.lower().replace('-', '_')}"]


    def get_interface_residues(
        self,
        domains: list,
        contact_threshold: int = 8,
        as_matrix: bool = False,
    ):
        """Get the interface residues / interface residue pairs for each interface in each rigid body.
        Each interface is defined by a tuple of chain IDs.

        Args:
            domains (list): List of rigid bodies
            contact_threshold (int, optional): Defaults to 8.
            as_matrix (bool, optional): Defaults to False.

        Returns:
            all_interface_residues (dict): Interface residues for each chain in each rigid body.
        """

        # CB coordinates of all residues for each chain.
        # TODO I think this is CA currently not CB
        # it is CB, I had mentioned it incorrectly in the gdoc and the docstring. See StructureParser.extract_perresidue_quantity function in Parser.py

        coords = np.array(self.coords_list)

        # all v all contact map
        contact_map = get_interaction_map(
            coords1=coords,
            coords2=coords,
            contact_threshold=contact_threshold,
            map_type="contact",
        )

        # (chain1, chain2): ([ch1_res_idx], [chain2_res_indices])
        all_interface_residues = defaultdict(dict)

        for rb_idx, rb_dict in enumerate(domains):

            unique_chains_in_rb = [chain_id for chain_id in rb_dict if rb_dict[chain_id]]

            if len(unique_chains_in_rb) < 2:
                continue

            unique_chain_pairs = sorted(list(combinations(unique_chains_in_rb, 2)))

            for chain_pair in unique_chain_pairs:

                ch1, ch2 = chain_pair

                # convert residue numbers to residue indices
                ch1_res_idx = [
                    self.num_to_idx[ch1][res_num] for res_num in rb_dict[ch1]
                ]
                ch2_res_idx = [
                    self.num_to_idx[ch2][res_num] for res_num in rb_dict[ch2]
                ]

                chain_pair_contact_map = np.zeros_like(contact_map)

                chain_pair_contact_map[np.ix_(ch1_res_idx, ch2_res_idx)] = 1
                chain_pair_contact_map[np.ix_(ch2_res_idx, ch1_res_idx)] = 1

                # 1 if both residues are in rigid body and make contact, 0 otherwise
                chain_pair_contact_map = chain_pair_contact_map * contact_map

                # if residue pairs are needed, do this
                if as_matrix:

                    all_interface_residues[rb_idx][chain_pair] = chain_pair_contact_map
                    continue

                contact_res_indices = np.unique(np.where(chain_pair_contact_map == 1)[0])

                ch1_contact_res_idx = list(
                    set(ch1_res_idx).intersection(contact_res_indices)
                )
                ch2_contact_res_idx = list(
                    set(ch2_res_idx).intersection(contact_res_indices)
                )

                if len(ch1_contact_res_idx+ch2_contact_res_idx) > 0:
                    all_interface_residues[rb_idx][chain_pair] = (
                        ch1_contact_res_idx,
                        ch2_contact_res_idx,
                    )
            #         print(
            #             f"Number of interface residues between {ch1},{ch2} ="
            #             f" {ch1}: {len(ch1_contact_res_idx)}"
            #             f" {ch2}: {len(ch2_contact_res_idx)}"
            #             f" Total: {len(ch1_contact_res_idx + ch2_contact_res_idx)}"
            #         )
            # print("-"*50)

        return all_interface_residues

    # TODO lots of duplication.
    # TODO merge functions that calculate interface pLDDT a) get_iplddt b) get_idr_plddt: calculate for all interfaces (chain pairs) by default. Have a flag to turn this off (turned on by default).
    # idr plddt is now merged in get_average_plddt function. chain_type = idr will give the idr plddt values.
    # also, for idr plddt, should the avg idr plddt for a chain be calculated on the residues in rigid body or all residues in the chain (that are in the prediction)?
    # iplddt can be calculated for all interfaces by default in get_iplddt. interface_type flag can be used to filter the interfaces.

    def get_ipLDDT(
        self,
        domains: list,
        interface_type: str = "any-any",
    ):
        """ Get the interface pLDDT values for each interface (of interface type) in each rigid body.

        Args:
            domains (list): List of rigid bodies
            interface_type (str, optional): Type of interface (any-any, idr-r, r-r, r-any, idr-idr, idr-any). Defaults to "any-any".

        Returns:
            all_iplddt_values (dict): Interface pLDDT values for each interface (of interface type) in each rigid body.
        """

        all_iplddt_values = defaultdict(dict)

        assert interface_type in ["any-any", "idr-r", "r-r", "r-any", "idr-idr", "idr-any"]; "Invalid interface type, choose from any-any, idr-r, r-r, r-any, idr-idr, idr-any"

        all_interface_residues = self.get_interface_residues(
            domains=domains,
            contact_threshold=8,
            as_matrix=False
        )

        for rb_idx, interface_res_dict in all_interface_residues.items():

            all_iplddt_values[rb_idx] = defaultdict(list)

            for chain_pair, interface_residues in interface_res_dict.items():

                if self.chain_pair_condition(chain_pair, interface_type):

                    ch1_contact_res_idx, ch2_contact_res_idx = interface_residues

                    ch1_plddt_vals = [
                        self.plddt_list[res_idx] for res_idx in ch1_contact_res_idx
                    ]
                    ch2_plddt_vals = [
                        self.plddt_list[res_idx] for res_idx in ch2_contact_res_idx
                    ]

                    all_iplddt_values[rb_idx][chain_pair] = ch1_plddt_vals + ch2_plddt_vals

                    print(
                        f"Average ipLDDT between {chain_pair[0]},{chain_pair[1]} = {np.mean(ch1_plddt_vals + ch2_plddt_vals):.2f}"
                    )

            _all_iplddt_vals = [
                iplddt
                for chain_pair_iplddt in all_iplddt_values[rb_idx].values()
                for iplddt in chain_pair_iplddt
            ]

            print("-" * 50)
            print(f"Average ipLDDT for interface type {interface_type} in rigid body {rb_idx} = {np.mean(_all_iplddt_vals):.2f}")
            print("-" * 50)

        return all_iplddt_values


    def get_average_pLDDT(
        self,
        domains: list,
        chain_type: str = "any",
    ):
        """Get the average pLDDT value of each chain of given chain_type in each rigid body.

        Args:
            domains (list): List of rigid bodies
            chain_type (str, optional): Type of chain (any, idr, r). Defaults to "any".

        Returns:
            avg_plddt_values (list): Average pLDDT values for each chain of given chain_type in each rigid body.
        """

        all_chain_plddt_dict = defaultdict(dict)

        assert chain_type in ["any", "idr", "r"]; "Invalid chain type, choose from any, idr, r"

        chain_ids = [chain_id for chain_id in self.lengths_dict.keys()]
        allowed_chain_ids = []

        if chain_type == "any":
            allowed_chain_ids = chain_ids
        elif chain_type == "idr":
            allowed_chain_ids = self.idr_chains
        elif chain_type == "r":
            allowed_chain_ids = [chain_id for chain_id in chain_ids if chain_id not in self.idr_chains]

        for rb_idx, rb_dict in enumerate(domains):

            all_chain_plddt_dict[rb_idx] = defaultdict(list)

            for chain_id, res_list in rb_dict.items():

                if chain_id in allowed_chain_ids:

                    all_chain_plddt_dict[rb_idx][chain_id].extend(
                        [self.plddt_list[self.num_to_idx[chain_id][res_num]] for res_num in res_list]
                    )
                    print(f"Average pLDDT of {chain_id} in rigid body {rb_idx} = {np.mean(all_chain_plddt_dict[rb_idx][chain_id]):.2f}")

            _all_chain_plddt_vals = [
                plddt
                for plddt_vals in all_chain_plddt_dict[rb_idx].values()
                for plddt in plddt_vals
            ]

            print("-"*50)
            print(f"Average pLDDT for all chains of chain_type '{chain_type}' in rigid body {rb_idx} = {np.mean(_all_chain_plddt_vals):.2f}")
            print("-"*50)

        return all_chain_plddt_dict


    def get_ipae(self, domains: list):
        """Get the interface PAE values for each interface in each rigid body.

        Args:
            domains (list): List of rigid bodies

        Returns:
            all_ipae_values (dict): Interface PAE values for each interface in each rigid body.
        """

        # we need residue pairs for this
        all_interface_residues = self.get_interface_residues(
            domains=domains,
            contact_threshold=8,
            as_matrix=True
        )

        all_ipae_values = defaultdict(dict)

        for rb_idx, interface_res_dict in all_interface_residues.items():

            for chain_pair, chain_pair_contact_map in interface_res_dict.items():

                chain_pair_pae_map = chain_pair_contact_map * self.pae
                chain_pair_pae_vals = chain_pair_pae_map[chain_pair_pae_map != 0].tolist()

                if len(chain_pair_pae_vals) > 0:
                    all_ipae_values[rb_idx][chain_pair] = chain_pair_pae_vals

                    print(
                        f"Average iPAE between {chain_pair[0]},{chain_pair[1]} = {np.mean(chain_pair_pae_vals):.2f}"
                    )

            _all_ipae_vals = [
                ipae_val
                for chain_pair_pae in all_ipae_values[rb_idx].values()
                for ipae_val in chain_pair_pae
            ]

            print("-"*50)
            print(f"Average iPAE for rigid body {rb_idx} is {np.mean(_all_ipae_vals):.2f}")
            print("-"*50)

        return all_ipae_values