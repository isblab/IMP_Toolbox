from collections import defaultdict
from af_pipeline.pae_to_domains.pae_to_domains import (
    parse_pae_file,
    domains_from_pae_matrix_igraph,
    domains_from_pae_matrix_networkx,
)
from utils import read_json
from af_pipeline.parser import AfParser


class Initialize:

    def __init__(
        self,
        output_dir: str,
        data_path: str,
        structure_path: str | None = None,
        modeled_region: dict = {},
    ):
        self.structure_path = structure_path
        self.data_path = data_path
        self.data = read_json(data_path)
        self.output_dir = output_dir
        self.modeled_region = modeled_region

        self.chain_ids = self.get_chain_ids() # [A, B, C]
        self.res_num_idx = self.get_protein_residues() # {A: [(1, 0), (2, 1)], B: [(1, 2), (2, 3)]}
        self.protein_lenghts, self.total_length = self.get_protein_lengths()

        # self.chains_of_i = self.chains_of_interest()
        # self.residues_of_i = self.residues_of_interest()

    def get_chain_ids(self):
        """Get the protein chains from the data file.
        list of unique chain ids.
        [A, B, C]
        """
        chain_ids = list(set(self.data["token_chain_ids"]))

        return chain_ids

    def get_protein_residues(self):
        """Get the protein residues from the data file.
        chain wise residues are stored in a dictionary.
        chain_id : [(res_num, res_idx), ...]
        """
        res_chains = self.data["token_chain_ids"] # A, A, B, ...
        res_nums = self.data["token_res_ids"] # 1, 2, 1, ...
        res_indices = list(range(len(res_nums))) # 0, 1, 2, ...
        res_num_idx = {}
        for ch_id, res_num, res_idx in zip(res_chains, res_nums, res_indices):
            if ch_id not in res_num_idx:
                res_num_idx[ch_id] = [(res_num, res_idx)]
            else:
                res_num_idx[ch_id].append((res_num, res_idx))
            # print(ch_id, res_num, res_idx)

        return res_num_idx

    def get_protein_lengths(self):
        """Get the protein lengths from the data file.
        chain_id : length
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
        list of chain ids. [A, B]
        """
        if len(self.modeled_region) == 0:
            chains_of_i = self.chain_ids
        else:
            chains_of_i = []
            for p_name, p_info in self.modeled_region.items():
                chains_of_i.append(p_info["chain"])

        return chains_of_i

    def residues_of_interest(self):
        """Get the residues of interest as specified in the modeled region.
        if modeled region is not specified, get all residues.
        chain_id : [(res_num, res_idx), ...]
        """
        if len(self.modeled_region) == 0:
            residues_of_i = self.res_num_idx
            # print(residues_of_i)
        else:
            residues_of_i = {}
            for p_name, p_info in self.modeled_region.items():
                chain = p_info["chain"]
                start, end = p_info["range"]
                offset = start - 1
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
        modeled_region: dict = {},
    ):
        super().__init__(
            output_dir=output_dir,
            data_path=data_path,
            structure_path=structure_path,
            modeled_region=modeled_region,
        )
        self.af_parser = None
        self.library = "igraph"

        self.pae_power = 1
        self.pae_cutoff = 5
        self.resolution = 0.5
        self.plddt_cutoff = 70

    def predict_domains(self):
        """Predict domains from a PAE file.
        Please Check pae_to_domains.py for more details

        Args:
            pae_file (str): Path to the PAE file

        Raises:
            ValueError: Invalid library specified. Use 'igraph' or 'networkx'

        Returns:
            domains (list): List of residues in each domain
            [[0, 1, 2, 3], [4, 5, 6], ...]
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

    def select_domains(self, domains):
        """Select the domains of interest and/or regions of interest within domains.
        domains = [[0, 1, 2, 3], [4, 5, 6], ...]
        residues_of_i = {A: [(1, 0), (2, 1)], B: [(1, 2), (2, 3)]}
        """

        residues_of_i = self.residues_of_interest()
        flat_residues_of_i = [
            res[1] # res_idx
            for chain, res_list
            in residues_of_i.items()
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
        residue numbering in the list is as per res_idx.
        """

        assert self.structure_path is not None; "Structure path not provided."

        self.af_parser = AfParser(self.structure_path, self.data_path)
        plddt = self.af_parser.get_ca_plddt()
        flat_plddt = [0] * self.total_length

        for chain, res_list in self.res_num_idx.items():
            for res in res_list:
                res_num, res_idx = res
                # print(res_idx)
                flat_plddt[res_idx] = plddt[chain][res_num-1]

        return flat_plddt

    def plddt_filtered_domains(self, domains):
        """Filter the residues based on plddt values.
        """

        selected_domains = self.select_domains(domains=domains)
        plddt = self.get_plddt()

        filtered_domains = []
        for domain in selected_domains:
            filtered_domain = []
            for res_idx in domain:
                if plddt[res_idx-1] >= self.plddt_cutoff:
                    filtered_domain.append(res_idx)
            filtered_domains.append(filtered_domain)

        return filtered_domains

    def get_rigid_bodies(self, domains, num_proteins=1):
        """Get the rigid bodies from the domains.
        """

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

        rigid_bodies = [rb for rb in rigid_bodies if len(rb) >= num_proteins]

        return rigid_bodies

# class ContactMap:

#     def __init__(self, chain1, chain2):
#         self.chain1 = chain1
#         self.chain2 = chain2

# class Restraints:

#     def __init__(self, chain1, chain2):
#         self.chain1 = chain1
#         self.chain2 = chain2