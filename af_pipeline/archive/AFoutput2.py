#! WIP
# Description: This file contains the classes to handle the output of the AlphaFold2/3 pipeline.
# RigibBodies -> find rigid bodies in the predicted structures

# TODO: combine af_pipeline functionality within this script
#   - Interaction class to get confident interactions (output: contact map, restraints) --> patches (meanshift)
#   - get interaction patches from the contact map for defining the restraints
#   - additional work -> adjust res_num for the output of af_pipeline scripts
# TODO: add one of the following
#   - cxc output selecting low plddt residues
#   - function to handle contiguous patches of low plddt residues
#   - change the bfactor of residues in pdb file so that all the atoms in a residue have the same bfactor (of representative atom)


from collections import defaultdict
import os
import numpy as np
from sklearn.cluster import MeanShift
from af_pipeline.pae_to_domains.pae_to_domains import (
    parse_pae_file,
    domains_from_pae_matrix_igraph,
    domains_from_pae_matrix_networkx,
)
from utils import get_key_from_res_range, read_json, write_json
from af_pipeline.Parser import AfParser, ResidueSelect
from af_pipeline.main import Interaction

"""
This script will contain modules for performing analysis for the AF2/3 prediction.
Define a class to perform housekeeping jobs like parsing, extracting coords, plddt, pae, etc.
Derive child classes for downstream analysis.
"""

class Initialize(AfParser):
    def __init__(
        self,
        data_file_path: str,
        struct_file_path: str | None = None,
        pred_selection: list = [],
        output_dir: str = "./af_output",
    ):

        # AF2/3 structure file path (optional for rigid bodies).
        self.struct_file_path = struct_file_path
        # AF2/3 structure data path.
        self.data_file_path = data_file_path
        # output directory to save the results.
        self.output_dir = output_dir
        # pred_selection: chain ids and regions to consider (optional).
        self.pred_selection = pred_selection

        self.get_data_attributes()

        if self.struct_file_path:
            self.get_structure_attributes()

    def get_data_attributes(self):

        self.data = self.get_data_dict()
        self.chain_ids = self.get_chain_ids(data=self.data)
        self.avg_pae = self.get_pae(data=self.data)
        self.lengths_dict = self.get_chain_lengths(data=self.data)


    def get_structure_attributes(self):
        """
        Extract the following from the input files:
            1. Ca coordinates.
            2. Ca pLDDT.
            3. Average PAE matrix.
        """

        # Residue positions of all residues for each chain.
        self.res_dict = self.get_residue_positions()
        self.lengths_dict = self.get_chain_lengths(self.res_dict)
        # Ca-coords of all residues for each chain.
        self.coords_dict = self.get_ca_coordinates()
        # Ca-plddt of all residues for each chain.
        self.plddt_dict = self.get_ca_plddt()
        # Average PAE matrix.
        # data = self.get_data_dict()
        # self.avg_pae = self.get_pae(data)
        # Get minPAE for each residue.
        self.min_pae_dict = self.get_min_pae(
            avg_pae = self.avg_pae,
            lengths_dict = self.lengths_dict,
            mask_intrachain = True,
            return_dict = True,
        )


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
            af_start, _ = p_select.get("af_region", (1, self.lengths_dict[chain]))
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
                start, end = p_select.get("model_region", (1, self.lengths_dict[chain]))
                af_start, af_end = p_select.get("af_region", (1, self.lengths_dict[chain]))

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
        - The residues are numbered as per res_num which is the actual residue number if
            the modeled region is specified.
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
            af_start, _ = p_select.get("af_region", (1, self.lengths_dict[chain]))
            chain_wise_offset[chain] = af_start - 1

        # note: rigid_bodies are as per actual residue numbering (UniProt)
        for idx, rb in enumerate(rigid_bodies):

            pdb_file = f"rigid_body_{idx}.pdb"
            res_dict = self.af_parser.get_residue_positions() # af prediction numbering

            rb_res_dict = {}
            for chain_id, res_list in res_dict.items():
                rb_res_dict[chain_id] = [
                    res
                    for res in res_list
                    if (res[0]+chain_wise_offset[chain_id]) in rb[chain_id]
                ]

        ResidueSelect(rb_res_dict)
        self.af_parser.save_pdb(ResidueSelect(rb_res_dict), f"{self.output_dir}/{pdb_file}")


class ContactMap(Interaction):

    def __init__(
        self,
        structure_path: str, # AF3 structure file
        data_path: str, # AF3 json file
        pred_selection: list = [], # offsets and chains to consider
        output_dir: str = "./af_output",
    ):
        super().__init__(structure_path, data_path)

        self.output_dir = output_dir
        self.pred_selection = pred_selection
        self.contact_map = None
        self.restraints = None

    def get_interacting_regions(self):
        """Get the interacting regions for all possible pairs of chains.

        The user input for model_region in pred_selection is as per UniProt numbering. \n
        However, the residue numbering in the data file is as per AF prediction (i.e. starting from 1). \n
        The offset is calculated as the difference between the model_region and af_region,\n
        so that start and end of the interacting region is calculated as per the AF prediction numbering.

        for e.g.
        if model_region is (635, 1118) and af_region is (600, 1118) for chain A, \n
        the offset is 600 - 1 = 599.\n
        and start = 635 - 599 = 36, \n
            end = 1118 - 599 = 519 \n
        i.e. we are concerened about residues 36 to 519 in chain A in AF3 prediction.

        Returns:
            interacting_regions (list): interacting regions in the format \n
                                        [{chain1: (start1, end1), chain2: (start2, end2)}, ...]
        """

        interacting_regions = []

        if len(self.pred_selection) == 0:  # consider all chains
            self.pred_selection = [
                {
                    "id": chain,
                    "model_region": (1, self.lengths_dict[chain]),
                    "af_region": (1, self.lengths_dict[chain]),
                }
                for chain in self.lengths_dict.keys()
            ]

        else:
            for p_select1 in self.pred_selection:
                for p_select2 in self.pred_selection:

                    chain_id1 = p_select1["id"]
                    chain_id2 = p_select2["id"]

                    offset = (
                        p_select1.get("af_region", (1, self.lengths_dict[chain_id1]))[0]
                        - 1
                    )
                    model_start, model_end = p_select1.get(
                        "model_region", (1, self.lengths_dict[chain_id1])
                    )

                    start1, end1 = model_start - offset, model_end - offset

                    offset = (
                        p_select2.get("af_region", (1, self.lengths_dict[chain_id2]))[0]
                        - 1
                    )
                    model_start, model_end = p_select2.get(
                        "model_region", (1, self.lengths_dict[chain_id2])
                    )

                    start2, end2 = model_start - offset, model_end - offset

                    if chain_id1 != chain_id2:
                        if any(
                            chain_id1 in d.keys() and chain_id2 in d.keys()
                            for d in interacting_regions
                        ):
                            continue
                        interacting_regions.append(
                            {
                                chain_id1: (start1, end1),
                                chain_id2: (start2, end2),
                            }
                        )

        return interacting_regions

    def get_contact_map(self, interacting_region):
        """Create a contact map."""

        contact_map = self.get_confident_interactions(
            interacting_region=interacting_region
        )
        return contact_map

    # def get_contacting_residues(self, interacting_region: dict):
    #     """Get the interacting residues in the interacting region.

    #     Args:
    #         interacting_region (dict): interacting region in the format \n
    #                                     {chain1: (start1, end1), chain2: (start2, end2)}

    #     Returns:
    #         residues (list): interacting residues in the format \n
    #                         [(res1, res2), ...]
    #     """

    #     contact_map = self.get_contact_map(interacting_region)
    #     residues = []

    #     contacting_residues = np.argwhere(contact_map == 1)
    #     chain1, chain2 = interacting_region.keys()

    #     ch_dict1 = next(item for item in self.pred_selection if item["id"] == chain1)
    #     ch_dict2 = next(item for item in self.pred_selection if item["id"] == chain2)

    #     chain1_offset = interacting_region[chain1][0] - 1 + ch_dict1["af_region"][0] - 1
    #     chain2_offset = interacting_region[chain2][0] - 1 + ch_dict2["af_region"][0] - 1

    #     # adjust the residue numbers
    #     for res in contacting_residues:
    #         res1, res2 = res
    #         res1 += chain1_offset + 1
    #         res2 += chain2_offset + 1
    #         residues.append((res1, res2))

    #     return residues

    def get_interacting_patches(
        self, contact_map: np.ndarray, interacting_region: dict, bandwidth=None
    ):
        """Get the interacting patches and the corrresponding segmented map

        Args:
            contact_map (np.ndarray): a binary contact map
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
        contacting_residues = np.argwhere(contact_map == 1)

        if len(contacting_residues) == 0:  # no interacting residues
            return contact_map, patches

        chain1, chain2 = interacting_region.keys()

        # ch_dict1 = next(item for item in self.pred_selection if item["id"] == chain1)
        # ch_dict2 = next(item for item in self.pred_selection if item["id"] == chain2)

        # # interacting region: start and end are as per AF prediction numbering
        # # we have to add the offset to get the residue numbers as per UniProt numbering
        # chain1_offset = (
        #     interacting_region[chain1][0]
        #     + ch_dict1.get("af_region", (1, self.lengths_dict[chain1]))[0]
        #     - 1
        # )
        # chain2_offset = (
        #     interacting_region[chain2][0]
        #     + ch_dict2.get("af_region", (1, self.lengths_dict[chain2]))[0]
        #     - 1
        # )

        ms = MeanShift(bandwidth=bandwidth)
        ms.fit(contacting_residues)
        labels = ms.labels_
        # cluster_centers = ms.cluster_centers_

        for res, label in zip(contacting_residues, labels):
            res1, res2 = res
            # res1 += chain1_offset
            # res2 += chain2_offset

            if label + 1 not in patches:
                patches[int(label + 1)] = {
                    chain1: [res1],
                    chain2: [res2],
                }

            else:
                patches[int(label + 1)][chain1].append(res1)
                patches[int(label + 1)][chain2].append(res2)

        for label, patch in patches.items():
            for chain in patch.keys():
                patch_residues = list(set(patch[chain]))
                patch[chain] = get_key_from_res_range(patch_residues)

        # segmented map of interacting patches
        segmented_map = np.zeros_like(contact_map)
        for i, point in enumerate(contacting_residues):
            segmented_map[point[0], point[1]] = labels[i] + 1

        return segmented_map, patches

    def save_patches(self, plots: str = None):
        """Save the interacting patches as plots and json file.

        Args:
            plots (str, optional): Defaults to "static". \n
                                "static" for static plots using matplotlib and \n
                                "interactive" for interactive plots using plotly.

        Returns:
            interacting_patches (dict): interacting patches in the format \n
                                    {chain1_chain2: {patch_id: {chain_id1: "1-50,375", chain_id2: "2-10"}}}
        """

        interacting_regions = self.get_interacting_regions()
        file_name = os.path.basename(self.struct_file_path).split(".")[0]
        os.makedirs(f"{self.output_dir}/{file_name}", exist_ok=True)

        interacting_patches = {}

        for interacting_region in interacting_regions:

            chain1, chain2 = interacting_region.keys()

            contact_map = self.get_contact_map(interacting_region)

            segmented_map, patches = self.get_interacting_patches(
                contact_map, interacting_region, bandwidth=None
            )

            interacting_patches[f"{chain1}_{chain2}"] = patches

            # plots
            # one can use generate_colors(n) from utils.py to get n distinct colors
            if plots in ["static", "interactive"]:
                ch_dict1 = next(item for item in self.pred_selection if item["id"] == chain1)
                ch_dict2 = next(item for item in self.pred_selection if item["id"] == chain2)

                p1_region = ch_dict1.get("model_region", (1, self.lengths_dict[chain1]))
                p2_region = ch_dict2.get("model_region", (1, self.lengths_dict[chain2]))

                ytick_vals = np.arange(
                    0, p1_region[1] - p1_region[0] + 1, p1_region[1] - p1_region[0]
                )
                xtick_vals = np.arange(
                    0, p2_region[1] - p2_region[0] + 1, p2_region[1] - p2_region[0]
                )
                ytick_labels = np.arange(p1_region[0], p1_region[1] + 1, p1_region[1] - p1_region[0])
                xtick_labels = np.arange(p2_region[0], p2_region[1] + 1, p2_region[1] - p2_region[0])

                if plots == "interactive":
                    import plotly.graph_objects as go

                    fig = go.Figure(
                        data=go.Heatmap(
                            z=segmented_map,
                            colorscale="hsv"
                        )
                    )
                    fig.update_layout(title=f"{chain1}_{chain2} patches")
                    fig.update_xaxes(title_text=f"{chain1}")
                    fig.update_yaxes(title_text=f"{chain2}")
                    fig.update_layout(
                        xaxis=dict(
                            tickvals=xtick_vals,
                            ticktext=xtick_labels,
                        ),
                        yaxis=dict(
                            tickvals=ytick_vals,
                            ticktext=ytick_labels,
                        ),
                    )
                    fig.write_html(
                        f"{self.output_dir}/{file_name}/{chain1}_{chain2}.html",
                        full_html=False,
                    )

                elif plots == "static":
                    from matplotlib import pyplot as plt

                    plt.imshow(
                        segmented_map,
                        cmap="hsv",
                        interpolation="nearest"
                    )
                    plt.title(f"{chain1}_{chain2} patches")
                    plt.xlabel(f"{chain1}")
                    plt.ylabel(f"{chain2}")
                    plt.xticks(ticks=xtick_vals, labels=xtick_labels)
                    plt.yticks(ticks=ytick_vals, labels=ytick_labels)
                    plt.savefig(
                        f"{self.output_dir}/{file_name}/{chain1}_{chain2}_patches.png"
                    )
                    plt.close()

                else:
                    raise ValueError("Invalid plot type. Use 'static' or 'interactive'")

        write_json(f"{self.output_dir}/{file_name}/patches.json", interacting_patches)

        return interacting_patches


# class Restraints:

#     def __init__(self, chain1, chain2):
#         self.chain1 = chain1
#         self.chain2 = chain2
