import textwrap
import pandas as pd
from collections import defaultdict
from IMP_Toolbox.analysis.rmf_to_xyzr import RMFToXYZRConverter
from IMP_Toolbox.chimerax.rmf_selection import get_rmf_to_residue_map
from IMP_Toolbox.utils.mutation_utils import split_missense_mutation
from IMP_Toolbox.utils.special_helpers import parse_topology_file
from IMP_Toolbox.utils.obj_helpers import get_key_from_res_range
from IMP_Toolbox.constants.structure_constants import RES_COLOR_MAP

class MutationMapper:
    """ Generic methods for mapping mutations to the structure """

    def __init__(
        self,
        mutation_impact_df: pd.DataFrame,
    ):
        self.mutation_impact_df = mutation_impact_df
        self.chimerax_res_sels = {
            "interface": defaultdict(set),
            "core": defaultdict(set),
            "exposed": defaultdict(set),
            "ptm": defaultdict(set),
        }
        self.afm_attr_scores = defaultdict(list)
        self._is_mapped = False

    @property
    def is_mapped(self):
        return self._is_mapped

    def check_if_mapped(self):
        if not self.is_mapped:
            raise RuntimeError(
                "Mutations not mapped. Please call get_mapped_residues()."
            )

    def get_mapped_residues(self):
        """ Map mutations to the structure based on the mutation_impact_df. """

        self.mapped_interface = self.get_residues_mapped_to_interface()
        self.mapped_core = self.get_residues_mapped_to_core()
        self.mapped_exposed = self.get_residues_mapped_to_exposed()
        self.mapped_ptm = self.get_residues_mapped_to_ptm()
        self._is_mapped = True

    def get_residues_mapped_to_core(self):
        """ Obtain residues mapped to the protein core from the mutation_impact_df.

        ## Returns:

        - **dict**:<br />
            A dictionary where keys are gene names and values are sets of residue numbers
            that are mapped to the protein core.
        """

        mutations_impacting_core = defaultdict(set)
        mutation_impact_df_core = self.mutation_impact_df[
            self.mutation_impact_df["Affected Interfaces or Core"] == "core"
        ]
        for idx, row in mutation_impact_df_core.iterrows():
            gene = row["Gene"]
            res_num = row["Residue Number"]
            mutations_impacting_core[gene].add(int(res_num))
        return mutations_impacting_core

    def get_residues_mapped_to_interface(self):
        """ Obtain residues mapped to the interface from the mutation_impact_df

        ## Returns:

        - **dict**:<br />
            A dictionary where keys are gene names and values are sets of residue numbers
            that are mapped to the interface.
        """

        mutations_impacting_interface = defaultdict(set)
        mutation_impact_df_interface = self.mutation_impact_df[
            (self.mutation_impact_df["Affected Interfaces or Core"] != "core") &
            (self.mutation_impact_df["Affected Interfaces or Core"] != "") &
            self.mutation_impact_df["Affected Interfaces or Core"].notnull()
        ]
        for idx, row in mutation_impact_df_interface.iterrows():
            gene = row["Gene"]
            res_num = row["Residue Number"]
            mutations_impacting_interface[gene].add(int(res_num))
        return mutations_impacting_interface

    def get_residues_mapped_to_ptm(self):
        """ Obtain residues mapped to PTM sites from the mutation_impact_df

        ## Returns:

        - **dict**:<br />
            A dictionary where keys are gene names and values are sets of residue numbers
            that are mapped to PTM sites.
        """

        mutations_impacting_ptm = defaultdict(set)
        mutation_impact_df_ptm = self.mutation_impact_df[
            self.mutation_impact_df["PTM"].notnull() &
            (self.mutation_impact_df["PTM"] != "")
        ]
        for idx, row in mutation_impact_df_ptm.iterrows():
            gene = row["Gene"]
            res_num = row["Residue Number"]
            mutations_impacting_ptm[gene].add(int(res_num))
        return mutations_impacting_ptm

    def get_residues_mapped_to_exposed(self):
        """ Obtain residues mapped to exposed residues from the mutation_impact_df

        ## Returns:

        - **dict**:<br />
            A dictionary where keys are gene names and values are sets of residue numbers
            that are mapped to exposed residues.
        """

        mutations_impacting_exposed = defaultdict(set)
        mutation_impact_df_exposed = self.mutation_impact_df[
            (self.mutation_impact_df["Affected Interfaces or Core"] == "") |
            (self.mutation_impact_df["Affected Interfaces or Core"].isnull())
        ]
        mutation_impact_df_exposed = mutation_impact_df_exposed[
            (mutation_impact_df_exposed["PTM"] == "") |
            (mutation_impact_df_exposed["PTM"].isnull())
        ]
        for idx, row in mutation_impact_df_exposed.iterrows():
            gene = row["Gene"]
            res_num = row["Residue Number"]
            mutations_impacting_exposed[gene].add(int(res_num))
        return mutations_impacting_exposed

    def generate_chimerax_cmds_for_mapped_mutations(
        self,
        outpath: str,
        mutated_regions: list = ["interface", "core", "exposed", "ptm"],
        color: bool = True,
    ):
        """ Generate a .cxc file with ChimeraX commands to select and optionally
        color residues based on their mapped mutation regions.

        ## Arguments:

        - **outpath (str)**:<br />
            Path to save the generated .cxc file.

        - **mutated_regions (list, optional):**:<br />
            List of mutated regions to include in the ChimeraX commands.

        - **color (bool, optional):**:<br />
            Whether to color the selected residues based on their mutated regions.
        """

        self.check_if_mapped()
        chimerax_cmds = []
        for mutated_region, selections in self.chimerax_res_sels.items():
            if mutated_region not in mutated_regions:
                continue
            for chain, res_nums in selections.items():
                res_range = get_key_from_res_range(res_nums)
                chimerax_cmds.append(f"sel add $1/{chain}:{res_range}")
                if color:
                    chimerax_cmds.append(
                        f"col $1/{chain}:{res_range} {RES_COLOR_MAP[mutated_region]}"
                    )

        if not outpath.endswith(".cxc"):
            outpath += ".cxc"

        with open(outpath, "w") as f:
            for cmd in chimerax_cmds:
                f.write(cmd + "\n")

        print("""
        Use the following commands in ChimeraX to execute the generated commands:
        open <your_ccm_file>. OR open <your_pdb_file>.pdb
        runscript <your_generated_commands>.cxc #<model_id>
        """)

    def generate_chimerax_attributes_for_mapped_mutations(
        self,
        outpath: str,
        operation: str = "average",
        match_mode: str = "any",
    ):

        self.check_if_mapped()
        chimerax_attrs = []
        for sel_key, scores in self.afm_attr_scores.items():
            if operation == "average":
                avg_score = sum(scores) / len(scores)
            elif operation == "max":
                avg_score = max(scores)
            elif operation == "min":
                avg_score = min(scores)
            afm_attr = f"\t{sel_key}\t{avg_score}"
            chimerax_attrs.append(afm_attr)

        attribute_template = textwrap.dedent(f"""
        attribute: afm_score
        match mode: {match_mode}
        recipient: residues
        """)
        chimerax_attrs = attribute_template + "\n".join(chimerax_attrs)

        if not outpath.endswith(".defattr"):
            outpath += ".defattr"

        with open(outpath, "w") as f:
            f.write(chimerax_attrs)

        print("""
        Use the following commands in ChimeraX to execute the generated attribute definitions:
        open <your_ccm_file>.rmf3 OR open <your_pdb_file>.pdb
        open <your_generated_attributes>.defattr format defattr models #<model_id>
        """)

class MutationMapperPDB(MutationMapper):
    """ Class for mapping mutations to PDB structures. """

    def __init__(
        self,
        entity_chain_map: dict,
        modeled_ranges: dict,
        mutation_impact_df: pd.DataFrame,
        af_missense_cutoff: float = 0.95,
    ):

        super().__init__(mutation_impact_df=mutation_impact_df)

        self.entity_chain_map = entity_chain_map
        self.modeled_ranges = modeled_ranges
        self.af_missense_cutoff = af_missense_cutoff
        self.molwise_chain_map = self.get_molwise_chain_map()

    def get_molwise_chain_map(self):
        """ Obtaine molecule-wise chains.

        ## Returns:

        - **dict**:<br />
            Dictionary where keys are molecule names and values are lists of
            chains corresponding to that molecule.
        """

        molwise_chain_map = defaultdict(list)
        for chain, mol in self.entity_chain_map.items():
            molwise_chain_map[mol].append(chain)
        return molwise_chain_map

    def get_mutations_mapped_to_pdb(self):
        """ Obtain mutations mapped to the given pdb.

        ## Returns:

        - **tuple**:<br />
            A tuple containing two dictionaries:
            1. A dictionary where keys are mutated regions and values are dictionaries
                with chains as keys and sets of residue numbers as values.
            2. A dictionary where keys are residue selections and values are lists of
                AlphaMissense scores corresponding to those selections.
        """

        self.get_mapped_residues()

        for idx, row in self.mutation_impact_df.iterrows():
            protein = row["Protein"]
            if len(self.modeled_ranges) == 2:
                start, end = self.modeled_ranges[protein]
            else:
                start, end = 0, 0
            gene = row["Gene"]
            mutations = row["Mutation"].split("\n")
            afm_scores = row["AlphaMissense Score"].split("\n")

            for mutation, afm_score in zip(mutations, afm_scores):

                wt, res_num, mut = split_missense_mutation(mutation)
                chains = self.molwise_chain_map.get(protein, [])

                if (
                    int(res_num) not in range(start, end+1) or
                    float(afm_score) < float(self.af_missense_cutoff)
                ):
                    continue

                for chain in chains:

                    if int(res_num) in self.mapped_interface.get(gene, set()):
                        self.chimerax_res_sels["interface"][chain].add(res_num)
                    elif int(res_num) in self.mapped_core.get(gene, set()):
                        self.chimerax_res_sels["core"][chain].add(res_num)
                    elif int(res_num) in self.mapped_exposed.get(gene, set()):
                        self.chimerax_res_sels["exposed"][chain].add(res_num)
                    if int(res_num) in self.mapped_ptm.get(gene, set()):
                        self.chimerax_res_sels["ptm"][chain].add(res_num)

                    self.afm_attr_scores[f"/{chain}:{res_num}"].append(afm_score)

        return self.chimerax_res_sels, self.afm_attr_scores

class MutationMapperRMF(MutationMapper):
    """ Class for mapping mutations to RMF cluster-center model. """

    def __init__(
        self,
        ccm_rmf_path: str,
        topology_file: str,
        mutation_impact_df: pd.DataFrame,
        modeled_ranges: dict = {},
        af_missense_cutoff: float = 0.95,
    ):

        super().__init__(mutation_impact_df=mutation_impact_df)

        self.ccm_rmf_path = ccm_rmf_path
        self.modeled_ranges = modeled_ranges
        self.af_missense_cutoff = af_missense_cutoff

        rmf_to_xyzr_converter = RMFToXYZRConverter(
            rmf_file=self.ccm_rmf_path,
            frame_subset="0",
            num_cores=1,
        )
        molwise_xyzr = rmf_to_xyzr_converter.convert_rmf_to_xyzr()
        xyzr_keys = list(molwise_xyzr.keys())
        chain_map, resolution_map = parse_topology_file(topology_file)
        rmf_to_residue_map = get_rmf_to_residue_map(
            xyzr_keys, chain_map, resolution_map
        )
        self.chain_map = chain_map
        self.rmf_to_residue_map = rmf_to_residue_map
        self.molwise_chain_map = self.get_molwise_chain_map()

    def get_molwise_chain_map(self):
        """ Obtain molecule-wise chains

        ## Returns:

        - **dict**:<br />
            Dictionary where keys are molecule names and values are lists of
            chains corresponding to that molecule.
        """

        molwise_chain_map = defaultdict(list)
        for molcopy, chain in self.chain_map.items():
            molwise_chain_map[molcopy.split("_")[0]].append(chain)
        return molwise_chain_map

    def get_mutations_mapped_to_ccm(self):
        """ Obtain mutations mapped to the give rmf file

        ## Returns:

        - **tuple**:<br />
            A tuple containing two dictionaries:
            1. A dictionary where keys are mutated regions and values are dictionaries
                with chains as keys and sets of residue numbers as values.
            2. A dictionary where keys are residue selections and values are lists of
                AlphaMissense scores corresponding to those selections.
        """

        self.get_mapped_residues()

        for idx, row in self.mutation_impact_df.iterrows():
            protein = row["Protein"]
            if len(self.modeled_ranges) == 2:
                start, end = self.modeled_ranges[protein]
            else:
                start, end = 0, 0
            gene = row["Gene"]
            mutations = row["Mutation"].split("\n")
            afm_scores = row["AlphaMissense Score"].split("\n")

            for mutation, afm_score in zip(mutations, afm_scores):

                wt, res_num, mut = split_missense_mutation(mutation)
                chains = self.molwise_chain_map.get(protein, [])

                if (
                    int(res_num) not in range(start, end+1) or
                    float(afm_score) < float(self.af_missense_cutoff)
                ):
                    continue

                for copy_idx, chain in enumerate(chains):

                    imp_particle = f"{protein}_{copy_idx}_{res_num}"
                    imp_res_num = self.rmf_to_residue_map[imp_particle][1]

                    if int(res_num) in self.mapped_interface.get(gene, set()):
                        self.chimerax_res_sels["interface"][chain].add(imp_res_num)
                    elif int(res_num) in self.mapped_core.get(gene, set()):
                        self.chimerax_res_sels["core"][chain].add(imp_res_num)
                    elif int(res_num) in self.mapped_exposed.get(gene, set()):
                        self.chimerax_res_sels["exposed"][chain].add(imp_res_num)
                    if int(res_num) in self.mapped_ptm.get(gene, set()):
                        self.chimerax_res_sels["ptm"][chain].add(imp_res_num)

                    self.afm_attr_scores[f"/{chain}:{imp_res_num}"].append(afm_score)

        return self.chimerax_res_sels, self.afm_attr_scores