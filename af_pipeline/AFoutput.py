#! WIP
# Description: This file contains the classes to handle the output of the AlphaFold2/3 pipeline.
# RigibBodies -> find rigid bodies in the predicted structures

# TODO: combine af_pipeline functionality within this script
#   - Interaction class to get confident interactions (output: contact map, restraints) --> patches (meanshift)
#   - additional work -> adjust res_num for the output of af_pipeline scripts
# TODO: add one of the following
#   - cxc output selecting low plddt residues
#   - function to handle contiguous patches of low plddt residues
#   - change the bfactor of residues in pdb file so that all the atoms in a residue have the same bfactor (of representative atom)


from collections import defaultdict
import os
from af_pipeline.pae_to_domains.pae_to_domains import (
    parse_pae_file,
    domains_from_pae_matrix_igraph,
    domains_from_pae_matrix_networkx,
)
from utils import read_json
from af_pipeline.parser import AfParser, ResidueSelect


class Initialize:

    def __init__(
        self,
        output_dir: str,
        data_path: str,
        structure_path: str | None = None,
        pred_selection: list = [],
    ):
        self.structure_path = structure_path
        self.data_path = data_path
        self.data = read_json(data_path)
        self.output_dir = output_dir
        self.pred_selection = pred_selection

        self.chain_ids = self.get_chain_ids() # [A, B, C]
        self.res_num_idx = self.get_protein_residues() # {A: [(1, 0), (2, 1)], B: [(1, 2), (2, 3)]}
        self.protein_lenghts, self.total_length = self.get_protein_lengths()  # {A: 2, B: 2}, 4

    def get_chain_ids(self):
        """Get the protein chains from the data file.
        list of unique chain ids.
        - [A, B, C]
        """

        chain_ids = list(set(self.data["token_chain_ids"]))

        return chain_ids

    def get_protein_residues(self):
        """Get the protein residues from the data file.
        chain wise residues are stored in a dictionary.
        - chain_id : [(res_num, res_idx), ...]
        """

        res_chains = self.data["token_chain_ids"]  # A, A, B, ...
        res_nums = self.data["token_res_ids"]  # 1, 2, 1, ...

        # adjust res_num for each chain as per af prediction
        chain_wise_offset = {}
        for p_select in self.pred_selection:
            chain = p_select["id"]
            af_start, _ = p_select["af_region"]
            chain_wise_offset[chain] = af_start - 1

        res_indices = list(range(len(res_nums)))  # 0, 1, 2, ...
        res_num_idx = {}
        for ch_id, res_num, res_idx in zip(res_chains, res_nums, res_indices):
            actual_res_num = res_num + chain_wise_offset[ch_id]
            if ch_id not in res_num_idx:
                res_num_idx[ch_id] = [(actual_res_num, res_idx)]
            else:
                res_num_idx[ch_id].append((actual_res_num, res_idx))

        return res_num_idx

    def get_protein_lengths(self):
        """Get the protein lengths from the data file.
        - chain_id : length
        """

        protein_lengths = {}
        total_length = 0
        for chain, residues in self.res_num_idx.items():
            protein_lengths[chain] = len(residues)
            total_length += len(residues)

        return protein_lengths, total_length

    def chains_of_interest(self):
        """Get the chains of interest as specified in the modeled region.
        if modeled region is not specified, get all chains.
        - list of chain ids. [A, B]
        """

        if len(self.pred_selection) == 0:
            chains_of_i = self.chain_ids
        else:
            chains_of_i = []
            for p_select in self.pred_selection:
                chains_of_i.append(p_select["id"])

        return chains_of_i

    def residues_of_interest(self):
        """Get the residues of interest as specified in the modeled region.
        if modeled region is not specified, get all residues.
        - chain_id : [(res_num, res_idx), ...]
        """

        if len(self.pred_selection) == 0:
            residues_of_i = self.res_num_idx
        else:
            residues_of_i = {}
            for p_select in self.pred_selection:

                chain = p_select["id"]
                start, end = p_select["model_region"]
                af_start, af_end = p_select["af_region"]

                assert af_start <= start <= end <= af_end; "model region not covered by af region"

                sel_res = range(start, end+1)
                sel_res = [
                    res_tuple
                    for res_tuple
                    in self.res_num_idx[chain]
                    if res_tuple[0] in sel_res
                ]
                residues_of_i[chain] = sel_res

        return residues_of_i


class RigidBodies(Initialize):

    def __init__(
        self,
        output_dir: str,
        data_path: str,
        structure_path: str | None = None,
        pred_selection: list = [],
    ):
        super().__init__(
            output_dir=output_dir,
            data_path=data_path,
            structure_path=structure_path,
            pred_selection=pred_selection,
        )
        self.af_parser = None
        self.library = "igraph"

        self.pae_power = 1
        self.pae_cutoff = 5
        self.resolution = 0.5
        self.plddt_cutoff = 70

    def predict_domains(self):
        """Predict domains from a PAE file.
        - Please Check pae_to_domains.py for more details

        Args:
            pae_file (str): Path to the PAE file

        Raises:
            ValueError: Invalid library specified. Use 'igraph' or 'networkx'

        Returns:
            domains (list): List of residues in each domain
            - [[0, 1, 2, 3], [4, 5, 6], ...]
        """

        pae_path = self.data_path
        pae_matrix = parse_pae_file(pae_path)

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
                domains[idx] = list(domain)

        return domains

    def select_roi_domains(self, domains):
        """Select the domains of interest and/or regions of interest within domains.
        domains = [[0, 1, 2, 3], [4, 5, 6], ...]
        residues_of_i = {A: [(1, 0), (2, 1)], B: [(1, 2), (2, 3)]}
        """

        residues_of_i = self.residues_of_interest()
        flat_residues_of_i = [
            res[1]  # res_idx
            for _, res_list in residues_of_i.items()
            for res in res_list
        ]

        selected_domains = []

        for domain in domains:
            sel_domain = []
            for res in domain:
                if res in flat_residues_of_i:
                    sel_domain.append(res)

            selected_domains.append(sel_domain)

        return selected_domains

    def get_plddt(self):
        """Get the plddt values of the residues.
        - residue numbering in the list is as per res_idx.
        """

        assert self.structure_path is not None; "Structure path not provided."

        if self.af_parser is None:
            self.af_parser = AfParser(self.structure_path, self.data_path)
        plddt = self.af_parser.get_ca_plddt()
        flat_plddt = [0] * self.total_length

        for chain, res_list in self.res_num_idx.items():
            for res in res_list:
                res_num, res_idx = res
                flat_plddt[res_idx] = plddt[chain][res_num-1]

        return flat_plddt

    def plddt_filtered_domains(self, domains, selected=False):
        """Filter the residues based on plddt values.
        """

        if selected:
            domains = self.select_roi_domains(domains)

        plddt = self.get_plddt()

        filtered_domains = []

        for domain in domains:
            filtered_domain = []
            for res_idx in domain:
                if plddt[res_idx-1] >= self.plddt_cutoff:
                    filtered_domain.append(res_idx)

            filtered_domains.append(filtered_domain)

        return filtered_domains

    def get_rigid_bodies(self, domains, num_proteins=1, selected=False):
        """Get the rigid bodies from the domains.
        - num_proteins: minimum number of proteins to be present in a rigid body
        - selected: if True, only the residues of interest are considered (model_region)

        Note:
        - The residues are numbered as per res_num which is the actual residue number if the modeled region is specified.
        - Otherwise, the residues are numbered as per AF prediction starting from 1.

        Returns:
        - list of rigid bodies
            - [{chain_id: [res_num, ...]}, ...]
        """

        assert num_proteins > 0; "Number of proteins in a rigid body should be greater than 0."

        if selected:
            domains = self.select_roi_domains(domains)

        rigid_bodies = []

        for idx, domain in enumerate(domains):
            chain_wise_domain_residues = defaultdict(list)

            for chain, res_list in self.res_num_idx.items():
                for res in res_list:
                    res_num, res_idx = res
                    if res_idx in domain:
                        chain_wise_domain_residues[chain].append(res_num)

            if len(domain) > 0:
                rigid_bodies.append(chain_wise_domain_residues)

        # all rigid bodies with at least num_proteins
        rigid_bodies = [rb for rb in rigid_bodies if len(rb) >= num_proteins]

        return rigid_bodies

    def save_rigid_bodies_pdb(self, rigid_bodies):
        """Save the rigid bodies as PDB files.
        """

        assert self.structure_path is not None; "Structure path not provided."


        if self.af_parser is None:
            self.af_parser = AfParser(self.structure_path, self.data_path)

        chain_wise_offset = {}
        for p_select in self.pred_selection:
            chain = p_select["id"]
            af_start, _ = p_select["af_region"]
            chain_wise_offset[chain] = af_start - 1

        for idx, rb in enumerate(rigid_bodies):
            pdb_file = f"rigid_body_{idx}.pdb"
            res_dict = self.af_parser.get_residue_positions()

            rb_res_dict = {}
            for chain_id, res_list in res_dict.items():
                rb_res_dict[chain_id] = [
                    res
                    for res in res_list
                    if (res[0]+chain_wise_offset[chain_id]) in rb[chain_id]
                ]

        ResidueSelect(rb_res_dict)
        self.af_parser.save_pdb(ResidueSelect(rb_res_dict), f"{self.output_dir}/{pdb_file}")

# class ContactMap:

#     def __init__(self, chain1, chain2):
#         self.chain1 = chain1
#         self.chain2 = chain2

# class Restraints:

#     def __init__(self, chain1, chain2):
#         self.chain1 = chain1
#         self.chain2 = chain2