# import sys
import pandas as pd
# from pprint import pprint
import time
import warnings
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import matplotlib.patches
import numpy as np
from tqdm import tqdm
from af_pipeline._Initialize import _Initialize
from af_pipeline.pae_to_domains.pae_to_domains import (
    domains_from_pae_matrix_igraph,
    domains_from_pae_matrix_networkx,
    domains_from_pae_matrix_label_propagation
)
import os
from collections import defaultdict
from af_pipeline.Parser import ResidueSelect
from utils import (
    get_key_from_res_range,
    save_structure_obj,
    convert_false_to_true,
    fill_up_the_blanks,
    get_interaction_map
)
from itertools import combinations, product
# from functools import wraps

# def time_it(func):
#     @wraps(func)
#     def wrapper(*args, **kwargs):
#         start_time = time.perf_counter()
#         result = func(*args, **kwargs)
#         end_time = time.perf_counter()
#         print(f"Function {func.__name__} took: {end_time - start_time:.4f} seconds")
#         return result
#     return wrapper


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
        **kwargs,
    ):

        super().__init__(
            data_file_path=data_path,
            struct_file_path=structure_path,
            af_offset=af_offset,
            **kwargs,
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
        num_res: int = 1,
        num_proteins: int = 1,
        plddt_filter: bool = True
    ):
        """Predict domains from a PAE file.
        - Three implementations are available:
            1. igraph based
            2. networkx based
            3. label_propagation based

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

        print("Predicting domains...")
        start_time = time.time()

        pae_matrix = self.pae

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
                random_seed=self.random_seed,
            )

        # domains is a list of frozensets or lists
        # each frozenset/list contains residue indices for residues in a domain
        for idx, domain in enumerate(domains):

            if isinstance(domain, frozenset):
                domain = list(domain)

            # rb_dict is a dictionary of rigid bodies
            # each rigid body is represented as a dictionary with chain_id as the key and a list of residue numbers as the value
            rb_dict = self.domain_to_rb_dict(domain=domain)

            # removing residues with pLDDT score below the cutoff
            if plddt_filter:
                rb_dict = self.filter_plddt(
                    rb_dict=rb_dict,
                    patch_threshold=self.patch_threshold,
                )

            domains[idx] = rb_dict

        # Remove domains with number of proteins less than `num_proteins`
        domains = [
            rb_dict
            for rb_dict in domains
            if len(rb_dict) >= num_proteins
        ]

        # Remove domains with number of residues less than `num_res`
        domains = [
            rb_dict
            for rb_dict in domains
            if sum([len(res_list) for res_list in rb_dict.values()]) >= num_res
        ]

        end_time = time.time()
        print(f"Done predicting pseudo-rigid domains in {end_time - start_time:.2f} seconds")

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
        structure_file_type: str = "pdb",
        no_plddt_filter_for_structure: bool = False,
        pae_plot: bool = False,
        rb_assessment: dict | None = None,
    ):
        """ Save the rigid bodies to a file and/or save the structure of the rigid bodies and assess the rigid bodies.
        - The rigid bodies are saved in a plain text format with the chain IDs and residue numbers.
        - The structure of the rigid bodies can be saved in PDB or CIF format. For rigid bodies with modifications, it is recommended to use PDB format.
        - The PAE plot can be saved to visualize the rigid bodies in the PAE matrix.
        - The rigid bodies can be assessed based on the interface residues, number of contacts, interface PAE and pLDDT, average PAE and plDDT and minimum PAE.
        - The assessment is saved in an Excel file.

        parameters for rigid body assessment:\n
        - `as_average`: \n
        whether to report only the average of assessment metric to the output file. Defaults to False. \n
        - `symmetric_pae`: \n
        whether to report a single average PAE value or assymetric PAE value for PAE assessment metrics. Defaults to False. \n

        Args:
            domains (list): list of rigid bodies, where each rigid body is a dictionary with chain IDs as keys and residue numbers as values.
            output_dir (str): Directory to save the output files.
            output_format (str, optional): Defaults to "txt". ("txt" or "csv")
            save_structure (bool, optional): Whether to save the structure of the rigid bodies. Defaults to True.
            structure_file_type (str, optional): File type to save the structure. Defaults to "pdb". ("pdb" or "cif")
            no_plddt_filter_for_structure (bool, optional): Whether to save the structure without filtering based on pLDDT. Defaults to False.
            pae_plot (bool, optional): Whether to save the PAE plot for the rigid bodies. Defaults to False.
            rb_assessment (dict | None, optional): Dictionary containing parameters for rigid body assessment.
        """

        dir_name = os.path.basename(self.struct_file_path).split(".")[0]
        output_dir = os.path.join(output_dir, dir_name)

        os.makedirs(output_dir, exist_ok=True)

        file_name = (
            os.path.basename(self.struct_file_path).split(".")[0] + "_rigid_bodies"
        )

        # txt or csv output format
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

        elif output_format == "csv":
            file_name += ".csv"
            output_path = os.path.join(output_dir, file_name)

            rows = []
            for idx, rb_dict in enumerate(domains):
                for chain_id, res_list in rb_dict.items():
                    if len(res_list) > 0:
                        rows.append({
                            "Rigid Body": idx,
                            "Chain ID": chain_id,
                            "Residues": get_key_from_res_range(res_range=res_list),
                        })

            df = pd.DataFrame(rows)
            df.to_csv(output_path, index=False)

        else:
            raise ValueError(
                f"Invalid output format: {output_format}. Use 'txt' or 'csv'."
            )

        ##################################################
        # Save the structure of the rigid bodies
        if save_structure:

            if structure_file_type == "cif":
                warnings.warn(
                    """
                    Protein or nucleotide modifications are stored as HETATM for which sequence connectivity
                    is lost in CIF format. \n
                    Please use PDB format to save the structure with modifications.
                    """
                )

            # Renumber the structure to match the actual sequence numbering if af_offset is provided
            structure = self.renumber.renumber_structure(
                structure=self.structureparser.structure,
            )

            for idx, rb_dict in enumerate(domains):

                # In the following case, the txt or csv ouput will have pLDDT filtered residues
                # but, the structure file will ignore this filter
                # use this flag when you don't want missing residues in the structure file
                if no_plddt_filter_for_structure:
                    for chain_id, res_list in rb_dict.items():
                        if len(res_list) > 0:
                            res_list = fill_up_the_blanks(res_list)
                            rb_dict[chain_id] = res_list

                output_path = os.path.join(output_dir, f"rigid_body_{idx}.{structure_file_type}")

                save_structure_obj(
                    structure=structure,
                    out_file=output_path,
                    res_select_obj=ResidueSelect(rb_dict),
                    save_type=structure_file_type,
                    preserve_header_footer=False,
                )

        ##################################################
        # Save the PAE plot for the rigid bodies
        # the region of the PAE matrix corresponding to the rigid bodies will be highlighted
        if pae_plot:
            for rb_idx, rb_dict in enumerate(domains):

                # patches are the highlighted rectangles in the PAE matrix
                patches = []

                for chain_id1, res_list1 in rb_dict.items():

                    for chain_id2, res_list2 in rb_dict.items():

                        res_idxs_1 = [
                            self.num_to_idx[chain_id1][res_num] for res_num in res_list1
                        ]
                        res_idxs_2 = [
                            self.num_to_idx[chain_id2][res_num] for res_num in res_list2
                        ]
                        res_idx_range_1 = get_key_from_res_range(res_range=res_idxs_1, as_list=True)
                        res_idx_range_2 = get_key_from_res_range(res_range=res_idxs_2, as_list=True)

                        for res_idx_1 in res_idx_range_1:

                            for res_idx_2 in res_idx_range_2:

                                if "-" in res_idx_1 and "-" in res_idx_2:
                                    xy_ = (int(res_idx_2.split("-")[0]), int(res_idx_1.split("-")[0])) # xy (0,0) coordinates for the rectangle
                                    h_ = int(res_idx_1.split("-")[1]) - int(res_idx_1.split("-")[0]) + 1 # patch height
                                    w_ = int(res_idx_2.split("-")[1]) - int(res_idx_2.split("-")[0]) + 1 # patch width

                                    if h_ > 0 and w_ > 0:
                                        patches.append(
                                            [xy_, h_, w_]
                                        )

                fig = plt.figure(figsize=(20, 20))
                plt.rcParams['font.size'] = 16
                plt.rcParams['axes.titlesize'] = 28
                plt.rcParams['axes.labelsize'] = 22
                plt.rcParams['xtick.labelsize'] = 13
                plt.rcParams['ytick.labelsize'] = 13
                plt.imshow(
                    self.pae,
                    # cmap="Greens_r",
                    cmap="Greys_r",
                    vmax=31.75,
                    vmin=0,
                    interpolation="nearest",
                    )

                for xy, h, w in patches:
                    rect = matplotlib.patches.Rectangle(
                        xy,
                        w,
                        h,
                        linewidth=0,
                        # edgecolor="green",
                        facecolor="lime",
                        alpha=0.5,
                    )
                    plt.gca().add_patch(rect)

                cumu_len = 0
                ticks = []
                ticks_labels = []

                for chain_id, p_length in self.lengths_dict.items():
                    if chain_id != "total":
                        cumu_len += p_length

                        if cumu_len != self.pae.shape[1]:
                            plt.axhline(y=cumu_len, color='red', linestyle='--', linewidth=0.75)
                            plt.axvline(x=cumu_len, color='red', linestyle='--', linewidth=0.75)

                        if self.af_offset is not None:

                            ticks_labels.extend(["\n" + f"{self.af_offset[chain_id][0]}" , f"{self.af_offset[chain_id][1]}" + "\n" ])
                            ticks.extend([cumu_len-p_length, cumu_len]) if cumu_len-p_length not in ticks else ticks.extend([cumu_len-p_length+1, cumu_len])

                        else:
                            ticks_labels.extend(["\n" + "1", f"{self.lengths_dict[chain_id]}" + "\n"])
                            ticks.extend([cumu_len-p_length, cumu_len]) if cumu_len-p_length not in ticks else ticks.extend([cumu_len-p_length+1, cumu_len])

                plt.xlim(0, self.pae.shape[0])
                plt.ylim(0, self.pae.shape[1])

                plt.gca().invert_yaxis()
                plt.yticks(ticks, ticks_labels)

                plt.xticks(ticks, ticks_labels, rotation=90, ha='center')
                plt.title(f"Predicted aligned error (PAE)", pad=20)
                plt.xlabel("Scored residue")
                plt.ylabel("Aligned residue")

                ax = plt.gca()

                divider = make_axes_locatable(ax)
                cax = divider.append_axes("bottom", size="5%", pad=1.2)
                plt.colorbar(
                    label="Predicted Alignment Error (PAE)",
                    orientation="horizontal",
                    cax=cax,
                )

                plt.savefig(os.path.join(output_dir, f"rigid_body_{rb_idx}.png"), transparent=True)
                plt.close(fig)

        ##################################################
        # Save the assessment of rigid bodies
        if rb_assessment:

            _start = time.time()
            print("Assessing rigid bodies...")

            assessment_file_name = (
                os.path.basename(self.struct_file_path).split(".")[0] + "_rb_assessment.xlsx"
            )
            save_path = os.path.join(output_dir, assessment_file_name)

            coords = np.array(self.coords_list)

            contact_map = get_interaction_map(
                coords1=coords,
                coords2=coords,
                contact_threshold=8,
                map_type="contact",
            )

            for rb_idx, rb_dict in enumerate(domains):

                rb_save_path = save_path.replace(
                    ".xlsx", f"_rb_{rb_idx}.xlsx"
                )

                rb_assess = RigidBodyAssessment(
                    rb_dict=rb_dict,
                    num_to_idx=self.num_to_idx,
                    idx_to_num=self.idx_to_num,
                    contact_map=contact_map,
                    plddt_list=self.plddt_list,
                    pae=self.pae,
                    lengths_dict=self.lengths_dict,
                    save_path=rb_save_path,
                    symmetric_pae=rb_assessment.get("symmetric_pae", False),
                    as_average=rb_assessment.get("as_average", False),
                    idr_chains=self.idr_chains,
                )

                rb_assess.save_rb_assessment()

            print(f"Time taken to save rigid body assessment: {time.time() - _start:.2f} seconds")


class RigidBodyAssessment:

    def __init__(
        self,
        rb_dict: dict,
        num_to_idx: dict,
        idx_to_num: dict,
        contact_map: np.ndarray,
        plddt_list: np.ndarray,
        pae: np.ndarray,
        lengths_dict: dict,
        save_path: str,
        **kwargs,
    ):
        """Initialize a rigid body object."""

        self.rb_dict = rb_dict
        self.num_to_idx = num_to_idx
        self.idx_to_num = idx_to_num
        self.symmetric_pae = kwargs.get("symmetric_pae", False)
        self.as_average = kwargs.get("as_average", False)
        self.idr_chains = kwargs.get("idr_chains", [])

        self.unique_chains = self.get_unique_chains()
        self.chain_pairs = self.get_chain_pairs()
        self.rb_res_binary_map = self.get_rb_res_binary_map(lengths_dict=lengths_dict)
        self.rb_res_pairs = self.get_rb_res_pairs()

        self.per_chain_plddt = self.get_per_chain_plddt(plddt_list=plddt_list)
        self.per_chain_avg_plddt = self.get_per_chain_avg_plddt()
        self.pairwise_pae = self.get_pairwise_pae(pae=pae)

        self.interface_res_pairs = self.get_interface_res_pairs(contact_map=contact_map)
        self.per_chain_interface_residues = self.get_per_chain_interface_residues()
        self.num_contacts = self.get_num_contacts()
        self.num_interface_residues = self.get_num_interface_residues()

        self.pairwise_ipae = self.get_pairwise_ipae(pae=pae)

        self.per_chain_iplddt = self.get_per_chain_iplddt(plddt_list=plddt_list)
        self.per_chain_avg_iplddt = self.get_per_chain_average_iplddt()

        self.pairwise_avg_iplddt = self.get_pairwise_avg_iplddt()

        self.pairwise_min_pae = self.get_pairwise_min_pae()
        self.pairwise_avg_pae = self.get_pairwise_avg_pae(symmetric_pae=self.symmetric_pae)
        self.pairwise_avg_ipae = self.get_pairwise_avg_ipae(symmetric_pae=self.symmetric_pae)

        self.overall_assessment = self.get_overall_assessment()

        self.save_path = save_path

    # @time_it
    def save_rb_assessment(self):
        """ Save the assessment of the rigid bodies to an Excel file.

        The assessment includes:
        - Per chain assessment: Average pLDDT, Average iLDDT, Number of interface residues, Chain type (IDR or R)
        - Per chain pair assessment: Number of interface residues, Number of contacts, Average PAE, Average iPAE, Minimum PAE, Average iLDDT for each chain, Chain type (IDR or R) for each chain
        - Overall assessment: Average pLDDT, Average iLDDT, Number of interface residues, Chain type (IDR or R)

        The assessment is saved in an Excel file with three sheets:
        - "Chain Wise Assessment": Contains per chain assessment data.
        - "Chain Pairwise Assessment": Contains per chain pair assessment data.
        - "Overall Assessment": Contains overall assessment data.
        """

        chain_wise_assessment_rows = []
        chain_pairwise_assessment_rows = []
        overall_assessment_rows = []

        for chain_id in self.unique_chains:
            chain_wise_assessment_rows.append({
                "Chain ID": chain_id,
                "Average pLDDT": self.per_chain_avg_plddt[chain_id],
                "Average ipLDDT": self.per_chain_avg_iplddt.get(chain_id, np.nan),
                "Number of Interface Residues": len(self.per_chain_interface_residues[chain_id]),
                "Chain Type": "IDR" if chain_id in self.idr_chains else "R",
            })

        for chain_pair in self.chain_pairs:
            chain1, chain2 = chain_pair

            if self.as_average:

                if self.symmetric_pae:
                    chain_pairwise_assessment_rows.append({
                        "Chain Pair": f"{chain1}-{chain2}",
                        "Number of Interface Residues": self.num_interface_residues[chain_pair],
                        "Number of Contacts": self.num_contacts[chain_pair],
                        "Average PAE": self.pairwise_avg_pae[chain_pair],
                        "Average iPAE": self.pairwise_avg_ipae[chain_pair],
                        "Minimum PAE": self.pairwise_min_pae[chain_pair],
                        "Average ipLDDT chain1": self.pairwise_avg_iplddt[chain_pair].get(chain1, np.nan),
                        "Average ipLDDT chain2": self.pairwise_avg_iplddt[chain_pair].get(chain2, np.nan),
                        "Chain Type 1": "IDR" if chain1 in self.idr_chains else "R",
                        "Chain Type 2": "IDR" if chain2 in self.idr_chains else "R",
                    })
                else:
                    chain_pairwise_assessment_rows.append({
                        "Chain Pair": f"{chain1}-{chain2}",
                        "Number of Interface Residues": self.num_interface_residues[chain_pair],
                        "Number of Contacts": self.num_contacts[chain_pair],
                        "Average PAE ij": self.pairwise_avg_pae[chain_pair]["ij"],
                        "Average PAE ji": self.pairwise_avg_pae[chain_pair]["ji"],
                        "Average iPAE ij": self.pairwise_avg_ipae[chain_pair]["ij"] if chain_pair in self.pairwise_avg_ipae else np.nan,
                        "Average iPAE ji": self.pairwise_avg_ipae[chain_pair]["ji"] if chain_pair in self.pairwise_avg_ipae else np.nan,
                        "Minimum PAE ij": self.pairwise_min_pae[chain_pair]["ij"] if chain_pair in self.pairwise_min_pae else np.nan,
                        "Minimum PAE ji": self.pairwise_min_pae[chain_pair]["ji"] if chain_pair in self.pairwise_min_pae else np.nan,
                        "Average ipLDDT chain1": self.pairwise_avg_iplddt[chain_pair].get(chain1, np.nan),
                        "Average ipLDDT chain2": self.pairwise_avg_iplddt[chain_pair].get(chain2, np.nan),
                        "Chain Type 1": "IDR" if chain1 in self.idr_chains else "R",
                        "Chain Type 2": "IDR" if chain2 in self.idr_chains else "R",
                    })

            else:

                if self.symmetric_pae:
                    for res1_idx, res2_idx in self.interface_res_pairs[chain_pair]:
                        ipae_val = (
                            self.pairwise_ipae[chain_pair]["ij"].get((res1_idx, res2_idx), np.nan) +
                            self.pairwise_ipae[chain_pair]["ji"].get((res2_idx, res1_idx), np.nan)
                        ) / 2
                        chain_pairwise_assessment_rows.append({
                            "Chain Pair": f"{chain1}-{chain2}",
                            "Residue Pair": f"{res1_idx}-{res2_idx}",
                            "iPAE": ipae_val,
                            "ipLDDT res1": self.per_chain_iplddt[chain1].get(res1_idx, np.nan),
                            "ipLDDT res2": self.per_chain_iplddt[chain2].get(res2_idx, np.nan),
                            "Chain Type 1": "IDR" if chain1 in self.idr_chains else "R",
                            "Chain Type 2": "IDR" if chain2 in self.idr_chains else "R",
                        })
                else:
                    for res1_idx, res2_idx in self.interface_res_pairs[chain_pair]:
                        ipae_ij = self.pairwise_ipae[chain_pair]["ij"].get((res1_idx, res2_idx), np.nan)
                        ipae_ji = self.pairwise_ipae[chain_pair]["ji"].get((res2_idx, res1_idx), np.nan)
                        chain_pairwise_assessment_rows.append({
                            "Chain Pair": f"{chain1}-{chain2}",
                            "Residue Pair": f"{res1_idx}-{res2_idx}",
                            "iPAE ij": ipae_ij,
                            "iPAE ji": ipae_ji,
                            "ipLDDT res1": self.per_chain_iplddt[chain1].get(res1_idx, np.nan),
                            "ipLDDT res2": self.per_chain_iplddt[chain2].get(res2_idx, np.nan),
                            "Chain Type 1": "IDR" if chain1 in self.idr_chains else "R",
                            "Chain Type 2": "IDR" if chain2 in self.idr_chains else "R",
                        })

        overall_assessment_keys = {
            "Number of Chains": "num_chains",
            "Number of Interacting Chain Pairs": "num_interacting_chain_pairs",
            "Number of Interface Residues": "num_interface_residues",
            "Number of Contacts": "num_contacts",
            "Average ipLDDT": "avg_iplddt",
            "Average IDR ipLDDT": "avg_idr_iplddt",
            "Average iPAE ij": "avg_ipae_ij",
            "Average iPAE ji": "avg_ipae_ji",
        }

        for col_head, key in overall_assessment_keys.items():
            if self.overall_assessment.get(key, np.nan) is not np.nan:
                overall_assessment_rows.append({
                    "Key": col_head,
                    "Value": self.overall_assessment.get(key)
                })

        chain_pairwise_assessment_df = pd.DataFrame(chain_pairwise_assessment_rows)
        chainwise_assessment_df = pd.DataFrame(chain_wise_assessment_rows)
        overall_assessment_df = pd.DataFrame(overall_assessment_rows)

        df_dict = {
            "chain_pairwise_assessment": chain_pairwise_assessment_df,
            "chainwise_assessment": chainwise_assessment_df,
            "overall_assessment": overall_assessment_df,
        }

        for k, df_ in df_dict.items():
            df_dict[k] = df_.fillna(np.nan)
            df_dict[k] = df_.map(lambda x: round(x, 2) if isinstance(x, (int, float)) else x)

        with pd.ExcelWriter(self.save_path, engine='openpyxl', mode='w') as writer:
            for sheet_name, df in df_dict.items():

                if df.empty:
                    warnings.warn(f"Skipping empty DataFrame for sheet: {sheet_name}")
                    continue

                df.to_excel(
                    writer,
                    sheet_name=sheet_name,
                    index=False,
                )

    # @time_it
    def get_unique_chains(self):
        """Get unique chains in the rigid body.

        Returns:
            unique_chains (list): List of unique chain IDs in the rigid body.
        """

        unique_chains = [
            chain_id
            for chain_id in self.rb_dict.keys()
            if len(self.rb_dict[chain_id]) > 0
        ]

        return unique_chains

    # @time_it
    def get_chain_pairs(self):
        """Get all unique chain pairs in the rigid body.

        Returns:
            chain_pairs (list): List of tuples containing unique chain pairs.
            Each tuple contains two chain IDs.
        """

        chain_pairs = list(combinations(self.unique_chains, 2))

        return [tuple(pair) for pair in chain_pairs]

    # @time_it
    def get_rb_res_binary_map(self, lengths_dict):
        """Get a binary map of residues in the rigid body.

        Returns:
            rb_res_binary_map (np.ndarray): A binary map of residues in the rigid body.
            The shape is (total_length, total_length) where total_length is the sum of lengths of all chains.
            The value is 1 if the residue is part of the rigid body, 0 otherwise.
        """

        total_len = lengths_dict.get("total", 0)
        rb_res_binary_map = np.zeros((total_len, total_len), dtype=int)
        all_rb_interface_res_idxs = []

        for chain_id, res_list in self.rb_dict.items():

            res_idxs = [self.num_to_idx[chain_id][res_num] for res_num in res_list]
            all_rb_interface_res_idxs.extend(res_idxs)

        all_rb_interface_res_idxs = np.unique(all_rb_interface_res_idxs)

        rb_res_binary_map[
            np.ix_(all_rb_interface_res_idxs, all_rb_interface_res_idxs)
        ] = 1

        return rb_res_binary_map

    # @time_it
    def get_rb_res_pairs(self):
        """Get all unique residue pairs in the rigid body.

        Returns:
            rb_res_pairs (defaultdict): A dictionary where keys are chain pairs (tuples) and values are lists of residue index pairs.
            Each residue index pair is a tuple of indices from the two chains in the rigid body.
        """

        rb_res_pairs = defaultdict(list)

        for chain_pair in self.chain_pairs:

            chain1, chain2 = chain_pair

            res1_list = self.rb_dict[chain1]
            res2_list = self.rb_dict[chain2]

            res1_idxs = [self.num_to_idx[chain1][res_num] for res_num in res1_list]
            res2_idxs = [self.num_to_idx[chain2][res_num] for res_num in res2_list]

            # Create pairs of residues from the two chains
            pairs = list(product(res1_idxs, res2_idxs))
            rb_res_pairs[chain_pair].extend(pairs)

        return rb_res_pairs

    # @time_it
    def get_interface_res_pairs(
        self,
        contact_map: np.ndarray,
    ):
        """ Get interface residue pairs from the contact map.

        Args:
            contact_map (np.ndarray): A binary contact map where 1 indicates a contact between residues and 0 indicates no contact.

        Returns:
            interface_res_pairs (defaultdict): A dictionary where keys are chain pairs (tuples) and values are lists of residue index pairs.
        """

        interface_res_pairs = defaultdict(list)
        contacting_res_indices = np.argwhere(contact_map == 1)

        for chain1, chain2 in tqdm(self.chain_pairs):

            res1_list = self.rb_dict[chain1]
            res2_list = self.rb_dict[chain2]

            res1_idxs = [self.num_to_idx[chain1][res_num] for res_num in res1_list]
            res2_idxs = [self.num_to_idx[chain2][res_num] for res_num in res2_list]

            mask1 = np.isin(contacting_res_indices[:, 0], res1_idxs)
            mask2 = np.isin(contacting_res_indices[:, 1], res2_idxs)

            mask = mask1 & mask2
            contacting_res_pairs = set(
                map(tuple, contacting_res_indices[mask])
            )

            if len(contacting_res_pairs) > 0:
                interface_res_pairs[(chain1, chain2)].extend(
                    list(contacting_res_pairs)
                )

        return interface_res_pairs

    # @time_it
    def get_per_chain_interface_residues(self):
        """Get interface residues for each chain.

        Returns:
            per_chain_interface_residues (defaultdict): A dictionary where keys are chain IDs and values are lists of residue indices.
            Each list contains the indices of residues that are part of any of the interfaces that the chain is involved in.
        """

        per_chain_interface_residues = defaultdict(list)

        for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

            chain1, chain2 = chain_pair

            for res1_idx, res2_idx in interacting_res_pairs:

                per_chain_interface_residues[chain1].append(
                    res1_idx
                ) if res1_idx not in per_chain_interface_residues[chain1] else None

                per_chain_interface_residues[chain2].append(
                    res2_idx
                ) if res2_idx not in per_chain_interface_residues[chain2] else None

        return per_chain_interface_residues

    # @time_it
    def get_num_interface_residues(self):
        """Get the number of interface residues for each chain pair.

        Returns:
            num_interface_residues (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the number of unique interface residues.
            Each key is a tuple of two chain IDs, and the value is the count of unique residues that interact between those chains.
        """

        num_interface_residues = defaultdict(int)

        for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

            unique_interface_residues = np.unique(np.array(interacting_res_pairs).flatten())
            num_interface_residues[chain_pair] = len(unique_interface_residues)

        return num_interface_residues

    # @time_it
    def get_num_contacts(self):
        """Get the number of contacts for each chain pair.

        Returns:
            num_contacts (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the number of contacts.
            Each key is a tuple of two chain IDs, and the value is the count of contacts between those chains.
        """

        num_contacts = defaultdict(int)

        for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

            num_contacts[chain_pair] = len(interacting_res_pairs)

        return num_contacts

    # @time_it
    def get_per_chain_plddt(self, plddt_list):
        """ Get per-chain pLDDT scores from a list of pLDDT scores.

        Args:
            plddt_list (list): A list of pLDDT scores for all residues in the structure.

        Returns:
            per_chain_plddt (defaultdict): A dictionary where keys are chain IDs and values are numpy arrays of pLDDT scores for residues in that chain.
        """

        per_chain_plddt = defaultdict(np.ndarray)

        for chain_id, res_list in self.rb_dict.items():

            res_idxs = [
                self.num_to_idx[chain_id][res_num] for res_num in res_list
            ]

            plddt_scores = np.array(plddt_list)[res_idxs]
            per_chain_plddt[chain_id] = plddt_scores

        return per_chain_plddt

    # @time_it
    def get_per_chain_avg_plddt(self):
        """ Get the average pLDDT score for each chain.

        Returns:
            per_chain_avg_plddt (dict): A dictionary where keys are chain IDs and values are the average pLDDT scores for that chain.
        """

        return {
            chain_id: np.mean(plddt_scores)
            for chain_id, plddt_scores in self.per_chain_plddt.items()
        }

    # @time_it
    def get_per_chain_iplddt(self, plddt_list):
        """ Get per-chain ipLDDT scores from a list of pLDDT scores.

        Args:
            plddt_list (list): A list of pLDDT scores for all residues in the structure.

        Returns:
            per_chain_iplddt (defaultdict): A dictionary where keys are chain IDs and values are dictionaries mapping residue indices to their pLDDT scores.
        """

        per_chain_iplddt = defaultdict(dict)

        for chain_id, interface_res_idxs in self.per_chain_interface_residues.items():

            for res_idx in interface_res_idxs:

                per_chain_iplddt[chain_id][res_idx] = plddt_list[res_idx]

        return per_chain_iplddt

    # @time_it
    def get_per_chain_average_iplddt(self):
        """ Get the average ipLDDT score for each chain.

        Returns:
            per_chain_avg_iplddt (dict): A dictionary where keys are chain IDs and values are the average ipLDDT scores for that chain.
        """

        return {
            chain_id: np.mean(list(iplddt_scores.values()))
            for chain_id, iplddt_scores in self.per_chain_iplddt.items()
        }

    # @time_it
    def get_pairwise_pae(self, pae):
        """ Get pairwise PAE values for each chain pair.

        Args:
            pae (np.ndarray): A 2D numpy array representing the predicted aligned error (PAE) matrix.

        Returns:
            pairwise_pae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are dictionaries containing PAE values for residue pairs.
        """

        pairwise_pae = defaultdict(np.ndarray)

        for chain_pair in self.chain_pairs:

            rb_chain_pair_res = self.rb_res_pairs[chain_pair]

            rb_pae_vals_ij = [
                pae[res1_idx, res2_idx] for res1_idx, res2_idx in rb_chain_pair_res
            ]

            rb_pae_vals_ji = [
                pae[res2_idx, res1_idx] for res1_idx, res2_idx in rb_chain_pair_res
            ]

            if len(rb_pae_vals_ij) > 0:
                pairwise_pae[chain_pair] = {
                    "ij": rb_pae_vals_ij,
                    "ji": rb_pae_vals_ji,
                }

        return pairwise_pae

    # @time_it
    def get_pairwise_avg_pae(self, symmetric_pae: bool = False):
        """ Get the average PAE for each chain pair.

        Args:
            symmetric_pae (bool, optional): If True, calculates the average PAE symmetrically for both directions (ij and ji).

        Returns:
            pairwise_avg_pae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the average PAE values.
        """

        if symmetric_pae:
            pairwise_avg_pae = defaultdict(float)
        else:
            pairwise_avg_pae = defaultdict(dict)

        for chain_pair in self.chain_pairs:
            if symmetric_pae:
                pairwise_avg_pae[chain_pair] = (
                    np.mean(
                        self.pairwise_pae[chain_pair]["ij"] +
                        self.pairwise_pae[chain_pair]["ji"]
                    ) / 2
                )
            else:
                pairwise_avg_pae[chain_pair]["ij"] = np.mean(self.pairwise_pae[chain_pair]["ij"])
                pairwise_avg_pae[chain_pair]["ji"] = np.mean(self.pairwise_pae[chain_pair]["ji"])

        return pairwise_avg_pae

    # @time_it
    def get_pairwise_min_pae(self, symmetric_pae: bool = False):
        """Get the minimum PAE for each chain pair.

        Args:
            symmetric_pae (bool, optional): If True, calculates the minimum PAE symmetrically for both directions (ij and ji).

        Returns:
            pairwise_min_pae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the minimum PAE values.
            If symmetric_pae is True, the minimum PAE is calculated as the minimum of both directions (ij and ji).
            If symmetric_pae is False, the minimum PAE is calculated separately for each direction.
        """

        if symmetric_pae:
            pairwise_min_pae = defaultdict(float)
        else:
            pairwise_min_pae = defaultdict(dict)

        for chain_pair, pae_dict in self.pairwise_pae.items():
            if symmetric_pae:
                pairwise_min_pae[chain_pair] = np.min(
                    np.min(pae_dict["ij"]), np.min(pae_dict["ji"])
                )
            else:
                pairwise_min_pae[chain_pair]["ij"] = np.min(pae_dict["ij"])
                pairwise_min_pae[chain_pair]["ji"] = np.min(pae_dict["ji"])

        return pairwise_min_pae

    # @time_it
    def get_pairwise_ipae(self, pae):
        """ Get pairwise iPAE values for each chain pair.

        Args:
            pae (np.ndarray): A 2D numpy array representing the predicted aligned error (PAE) matrix.

        Returns:
            pairwise_ipae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are dictionaries containing iPAE values for residue pairs.
        """

        pairwise_ipae = defaultdict(dict)

        for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

            pairwise_ipae[chain_pair] = {
                "ij" : {
                    (res1_idx, res2_idx): pae[res1_idx, res2_idx]
                    for res1_idx, res2_idx in interacting_res_pairs
                },
                "ji" : {
                    (res2_idx, res1_idx): pae[res2_idx, res1_idx]
                    for res1_idx, res2_idx in interacting_res_pairs
                }
            }

        return pairwise_ipae

    # @time_it
    def get_pairwise_avg_ipae(self, symmetric_pae: bool = False):
        """ Get the average iPAE for each chain pair.

        Args:
            symmetric_pae (bool, optional): If True, calculates the average iPAE symmetrically for both directions (ij and ji).

        Returns:
            pairwise_avg_ipae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the average iPAE values.
        """

        if symmetric_pae:
            pairwise_avg_ipae = defaultdict(float)
        else:
            pairwise_avg_ipae = defaultdict(dict)

        for chain_pair, ipae_dict in self.pairwise_ipae.items():

            if symmetric_pae:
                pairwise_avg_ipae[chain_pair] = (
                    np.mean(list(ipae_dict["ij"].values()) + list(ipae_dict["ji"].values())) / 2
                )
            else:
                pairwise_avg_ipae[chain_pair]["ij"] = np.mean(list(ipae_dict["ij"].values()))
                pairwise_avg_ipae[chain_pair]["ji"] = np.mean(list(ipae_dict["ji"].values()))

        return pairwise_avg_ipae

    # @time_it
    def get_pairwise_avg_iplddt(self):
        """ Get the average ipLDDT for each chain pair.

        Returns:
            pairwise_avg_iplddt (defaultdict): A dictionary where keys are chain pairs (tuples) and values are dictionaries containing average ipLDDT values for each chain in the pair.
        """

        pairwise_avg_iplddt = defaultdict(dict)

        for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

            chain1, chain2 = chain_pair

            iplddt1_values = [
                self.per_chain_iplddt[chain1].get(res1_idx, np.nan)
                for res1_idx, res2_idx in interacting_res_pairs
            ]

            iplddt2_values = [
                self.per_chain_iplddt[chain2].get(res2_idx, np.nan)
                for res1_idx, res2_idx in interacting_res_pairs
            ]

            pairwise_avg_iplddt[chain_pair][chain1] = np.mean(iplddt1_values)
            pairwise_avg_iplddt[chain_pair][chain2] = np.mean(iplddt2_values)

        return pairwise_avg_iplddt


    # @time_it
    def get_overall_assessment(self):
        """ Get overall assessment of the rigid body.

        Returns:
            overall_assessment (dict): A dictionary containing overall statistics about the rigid body.
            It includes the number of chains, number of interacting chain pairs, number of interface residues,
            number of contacts, average ipLDDT, average IDR ipLDDT, average iPAE ij, and average iPAE ji.
        """

        overall_assessment = {}

        overall_assessment["num_chains"] = len(self.unique_chains)

        overall_assessment["num_interacting_chain_pairs"] = len(self.interface_res_pairs)

        overall_assessment["num_interface_residues"] = sum(
            len(res_list)
            for res_list in self.per_chain_interface_residues.values()
        )

        overall_assessment["num_contacts"] = sum(
            len(contact_pairs)
            for contact_pairs in self.interface_res_pairs.values()
        )

        global_iplddt_scores = [
            iplddt
            for iplddt_scores in self.per_chain_iplddt.values()
            for iplddt in iplddt_scores.values()
        ]

        overall_assessment["avg_iplddt"] = (
            np.mean(global_iplddt_scores) if global_iplddt_scores else np.nan
        )

        global_idr_iplddt_scores = [
            iplddt
            for chain_id, iplddt_scores in self.per_chain_iplddt.items()
            for iplddt in iplddt_scores.values()
            if chain_id in self.idr_chains
        ]

        overall_assessment["avg_idr_iplddt"] = (
            np.mean(global_idr_iplddt_scores) if global_idr_iplddt_scores else np.nan
        )

        global_ipae_ij_scores = [
            ipae
            for ipae_dict in self.pairwise_ipae.values()
            for ipae in ipae_dict["ij"].values()
        ]

        global_ipae_ji_scores = [
            ipae
            for ipae_dict in self.pairwise_ipae.values()
            for ipae in ipae_dict["ji"].values()
        ]

        overall_assessment["avg_ipae_ij"] = (
            np.mean(global_ipae_ij_scores) if global_ipae_ij_scores else np.nan
        )

        overall_assessment["avg_ipae_ji"] = (
            np.mean(global_ipae_ji_scores) if global_ipae_ji_scores else np.nan
        )

        return overall_assessment
