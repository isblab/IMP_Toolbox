import os
import tempfile
import pandas as pd
import Bio.PDB.Structure
from pathlib import Path
from Bio.PDB import FastMMCIFParser, PDBParser
from collections import defaultdict
from IMP_Toolbox.structure import (
    get_burial_info,
    get_per_chain_residues,
    split_structure_by_chain,
)
from IMP_Toolbox.utils.mutation_utils import split_missense_mutation
from IMP_Toolbox.utils.obj_helpers import get_res_range_from_key
from IMP_Toolbox.mutations.af_missense import (
    PairwiseSequenceAlignment,
    fetch_fasta_dict_for_af_missense
)
from IMP_Toolbox.constants.structure_constants import ProteinCoreCriteria
from IMP_Toolbox.ptms.ptmd import get_ptmd_data
from IMP_Toolbox.constants.mutation_constants import AA_MAP_ONE_TO_THREE
from IMP_Toolbox.analysis.interaction import Interaction
from IMP_Toolbox.structure.tools import save_structure_obj

class VariantsEnricher:
    """ Enrich variants with structural information including:
    - Affected interfaces
    - Disease association
    - Residue depth
    - Relative solvent accessibility
    """

    def __init__(
        self,
        chain_mappings: list,
        clinvar_mutations_df: pd.DataFrame,
        protein_sequences: dict,
        protein_gene_map: dict,
        protein_uniprot_map: dict,
        interaction_map_dir: str,
        alpha_missense_dir: str,
        pairwise_alignments_dir: str,
        include_ptm: bool = False,
        include_vus: bool = False,
        af_missense_cutoff: float = 0.95,
        ptmd_dir: str = "./ptms",
        msms_executable: str = os.environ.get("MSMS_BIN", "/home/$USER/Software/msms_i86_64Linux2_2.6.1/msms.x86_64Linux2.2.6.1"),
    ):

        self.chain_mappings = chain_mappings
        self.clinvar_mutations_df = clinvar_mutations_df
        self.protein_sequences = protein_sequences
        self.protein_gene_map = protein_gene_map
        self.gene_protein_map = {v: k for k, v in protein_gene_map.items()}
        self.protein_uniprot_map = protein_uniprot_map
        self.interaction_map_dir = interaction_map_dir
        self.msms_executable = msms_executable
        self.include_ptm = include_ptm
        self.include_vus = include_vus
        self.af_missense_cutoff = af_missense_cutoff

        assert os.path.exists(self.msms_executable), f"MSMS executable not found at {self.msms_executable}. Please provide correct path."

        ############################################################################
        # Get interface residues for each protein from interaction patches
        ############################################################################
        interacting_patches_csvs = [
            os.path.join(interaction_map_dir, f) for f in os.listdir(interaction_map_dir)
            if f.startswith("patches_") and f.endswith(".csv")
        ]
        self.interface_dict = Interaction.extract_interface_residues_from_patches(
            interacting_patches_csvs=interacting_patches_csvs,
        )

        ############################################################################
        # Get PTM data from PTMD (https://ptmd.biocuckoo.cn/) if include_ptm is True
        ############################################################################
        if self.include_ptm:
            get_ptmd_data(
                savedir=Path(ptmd_dir) / "ptmd_data",
                key="ptmd_total",
                overwrite=False,
                reextract=False,
            )
            self.ptmd_total_df = pd.read_csv(Path(ptmd_dir) / "ptmd_data" / "ptmd_total" / "Total.txt", sep="\t")
        self.enriched_variants = defaultdict(dict)
        self.mutated_residues_per_chain = defaultdict(set)

        self.clinvar_mutations = clinvar_mutations_df["Mutation"].tolist()
        self.genes = clinvar_mutations_df["Gene"].tolist()
        self.disease_assoc = clinvar_mutations_df["Disease association"].tolist()
        self.clinvar_pathogenicity = clinvar_mutations_df["ClinVar clinical significance"].tolist()
        # self.afmissense_pathogenicity = clinvar_mutations_df["AlphaMissense pathogenicity"].tolist()
        self.afmissense_scores = clinvar_mutations_df["AlphaMissense score"].tolist()

        self.alpha_missense_dir = alpha_missense_dir
        self.pairwise_alignments_dir = pairwise_alignments_dir

    def enrich_variants(self):
        """

        ## Returns:

        - **enriched_variants (dict)**:<br />
            A dictionary mapping (gene, residue number) tuples to dictionaries containing
            enriched information about each mutation, including:
            - "Gene": Gene name associated with the mutation.
            - "Residue Number": Residue number of the mutation.
            - "Mutation": Set of mutation annotations (e.g., "R123C").
            - "Affected Interfaces": Set of interfaces potentially affected by the mutation.
            - "Disease association": Set of diseases associated with the mutation.
            - "Rigid Body": Set of rigid bodies that the mutated residue is part of.
            - "Secondary Structure": Set of secondary structure annotations for the mutated residue.
            - "RSA Value": Set of Relative Solvent Accessibility (RSA) values for the mutated residue.
            - "Residue Depth": Set of residue depth values for the mutated residue.
            - "CAB Depth": Set of CAB depth values for the mutated residue.
            - "PTM": Set of post-translational modifications (PTMs) that affect the mutated residue.
        """

        for chain_mapping in self.chain_mappings:

            structure_path = chain_mapping["structure_path"]
            entity_chain_map = chain_mapping["entity_chain_map"]
            modeled_range = chain_mapping.get("modeled_range", None)
            idr_chains = chain_mapping.get("idr_chains", [])

            file_extension = os.path.splitext(structure_path)[1].lower()
            if file_extension == ".cif":
                parser = FastMMCIFParser(QUIET=True)
            elif file_extension == ".pdb":
                parser = PDBParser(QUIET=True)
            else:
                raise ValueError("Unsupported file format. Please provide a .cif or .pdb file.")

            structure = parser.get_structure("structure", structure_path)
            rigid_body_ranges = get_per_chain_residues(structure)

            ########################################################################
            # Identify mutated residues in the rigid bodies
            ########################################################################
            mutated_residues = self.get_mutated_rigid_body_residues(
                structure=structure,
                entity_chain_map=entity_chain_map,
            )

            chain_structures = split_structure_by_chain(
                structure=structure,
            )

            # RSA and Residue depth calculations are done for each chain separately
            for ch_id, ch_structure in chain_structures.items():

                entity_ = entity_chain_map[ch_id]
                uniprot_ = self.protein_uniprot_map.get(entity_, None)
                uniprot_base = uniprot_.split("-")[0] if uniprot_ is not None else None

                psa_map, entity_ptm_df = self.prepare_ptm_df(
                    entity_=entity_,
                    uniprot_base=uniprot_base
                )

                ####################################################################
                # Calculate RSA and Residue depth for the chain
                ####################################################################
                temp_dir = tempfile.mkdtemp()
                temp_structure_path = os.path.join(temp_dir, f"{Path(structure_path).stem}_{ch_id}.pdb")
                save_structure_obj(
                    structure=ch_structure,
                    out_file=temp_structure_path,
                    save_type="pdb",
                )
                burial_info_df = get_burial_info(
                    structure=ch_structure,
                    ignore_chains=idr_chains,
                    structure_path=temp_structure_path,
                    residue_selector={ch_id: list(mutated_residues[ch_id])},
                    include_residue_depth=True,
                    msms_executable=self.msms_executable,
                    entity_chain_map=entity_chain_map,
                )

                os.remove(temp_structure_path)

                burial_info_dict = burial_info_df.set_index(
                    ["chain_id", "res_num"]
                ).to_dict(orient="index")

                ####################################################################
                # Iterate through ClinVar data and fill in the dictionary
                # with per-residue information including:
                # - Affected interfaces
                # - Disease association
                # - Parent rigid body
                # - Burial info (secondary structure, RSA, residue depth, CAB depth)
                ####################################################################
                self.enrich_chain_annotations(
                    modeled_range,
                    rigid_body_ranges,
                    ch_id,
                    entity_,
                    psa_map,
                    entity_ptm_df,
                    temp_structure_path,
                    burial_info_dict,
                    self.af_missense_cutoff,
                )

        return self.enriched_variants

    def get_affected_interfaces(
        self,
        gene: str,
        entity: str,
        res_num: int
    ) -> list:
        """ Identify which interfaces could be affected by a mutation.

        ## Arguments:

        - **gene (str)**:<br />
            Gene name associated with the mutation.

        - **entity (str)**:<br />
            Protein entity name in the structure associated with the mutation.

        - **res_num (int)**:<br />
            Residue number of the mutation.

        ## Returns:

        - **list**:<br />
            List of affected interfaces, where each entry is a string in the format
            "Protein1_Protein2 (res1_start-res1_end:res2_start-res2_end)" indicating
            the interface and the residue ranges involved in the interaction.
        """

        affected_interfaces = []
        potential_affected_interfaces = [
            pair for pair in self.interface_dict.keys()
            # if self.imp_protein_map.get(entity, entity) in pair.split("_")
            if entity in pair.split("_")
        ]
        for pair in potential_affected_interfaces:
            prot1, prot2 = pair.split("_")
            # gene1, gene2 = self.imp_gene_map.get(prot1, prot1), self.imp_gene_map.get(prot2, prot2)
            gene1, gene2 = self.protein_gene_map.get(prot1, prot1), self.protein_gene_map.get(prot2, prot2)
            if gene == gene1 and gene != gene2:
                idx_ = 0
            elif gene == gene2 and gene != gene1:
                idx_ = 1
            elif gene == gene1 and gene == gene2:
                idx_ = [0,1]
            else:
                continue

            affected_res_patches = self.interface_dict[pair]
            for patch_ in affected_res_patches:
                r_range_1 = get_res_range_from_key(patch_[0])
                r_range_2 = get_res_range_from_key(patch_[1])
                if isinstance(idx_, list):
                    if res_num in r_range_1 or res_num in r_range_2:
                        affected_interfaces.append(f"{pair} ({":".join(patch_)})")
                        break
                if idx_ == 0 and res_num in r_range_1:
                    affected_interfaces.append(f"{pair} ({":".join(patch_)})")
                    break
                elif idx_ == 1 and res_num in r_range_2:
                    affected_interfaces.append(f"{pair} ({":".join(patch_)})")
                    break

        return affected_interfaces

    def skip_mutation(
        self,
        entity_: str,
        gene: str,
        clinv_p: str,
        afm_s: float,
        modeled_range: dict,
        ch_id: str,
        res_num: int,
        af_missense_cutoff: float = 0.95, # only considered in case of VUS
    ) -> bool:
        """ Determine whether to skip a mutation based on its ClinVar clinical
        significance and AlphaMissense pathogenicity prediction.

        ## Arguments:

        - **entity_ (str)**:<br />
            Protein entity name in the structure associated with the mutation.

        - **gene (str)**:<br />
            Gene name associated with the mutation.

        - **clinv_p (str)**:<br />
            ClinVar clinical significance of the mutation.

        - **afm_s (str)**:<br />
            AlphaMissense pathogenicity prediction for the mutation.

        - **modeled_range (dict)**:<br />
            Dictionary mapping chain IDs to tuples of (start, end) residue numbers
            representing the modeled range of residues for each chain.

        - **ch_id (str)**:<br />
            Chain ID in the structure associated with the mutation.

        - **res_num (int)**:<br />
            Residue number of the mutation.

        - **af_missense_cutoff (float)**:<br />
            Cutoff value for AlphaMissense pathogenicity prediction.

        ## Returns:

        - **bool**:<br />
            True if the mutation should be skipped based on the pathogenicity criteria,
            False otherwise.
        """

        # if not considering VUS, skip mutations that are not "Pathogenic"
        # or "Likely pathogenic" in ClinVar
        if self.include_vus is False and clinv_p not in [
            "Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"
        ]:
            return True

        # if considering VUS, only take Clinvar "Uncertain significance" mutations
        # if they are "Likely pathogenic" in AlphaMissense
        elif (
            self.include_vus is True and
            clinv_p == "Uncertain significance" and
            float(afm_s) < af_missense_cutoff
        ) or (
            self.include_vus is True and
            clinv_p in ["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]
        ):
            return True

        if self.protein_gene_map.get(entity_) != gene:
            return True

        if modeled_range is not None:
            start, end = modeled_range.get(ch_id, (None, None))
            if res_num < start or res_num > end:
                return True

        return False

    def enrich_chain_annotations(
        self,
        modeled_range: dict,
        rigid_body_ranges: dict,
        ch_id: str,
        entity_: str,
        psa_map: dict,
        entity_ptm_df: pd.DataFrame,
        temp_structure_path: str,
        burial_info_dict: dict,
        af_missense_cutoff: float,
    ):
        """ Enrich the annotations for a specific chain in the structure with
        structural information.

        ## Arguments:

        - **modeled_range (dict)**:<br />
            Dictionary mapping chain IDs to tuples of (start, end) residue numbers
            representing the modeled range of residues for each chain.

        - **rigid_body_ranges (_type_)**:<br />
            Dictionary mapping chain IDs to lists of residue numbers representing
            the rigid body ranges for each chain.

        - **ch_id (str)**:<br />
            Chain ID in the structure associated with the mutation.

        - **entity_ (str)**:<br />
            Protein entity name in the structure associated with the mutation.

        - **psa_map (_type_)**:<br />
            Dictionary mapping residue numbers in the reference sequence to
            residue numbers in the modeled sequence, used for mapping PTM
            annotations.
            This is useful if the PTM is mapped on a different isoform than the
            one that is modeled.

        - **entity_ptm_df (_type_)**:<br />
            DataFrame containing PTM annotations for the protein entity, including
            columns for position, UniProt ID, PTM type, residue, mutation site, and
            experimental verification status.

        - **temp_structure_path (_type_)**:<br />
            Path to a temporary PDB file containing the structure of the chain,
            used for calculating burial information such as secondary structure,
            relative solvent accessibility (RSA), residue depth, and CAB depth.

        - **burial_info_dict (_type_)**:<br />
            Dictionary containing burial information for residues in the chain,
            indexed by (chain_id, residue_number) tuples, with values containing
            secondary structure, RSA value, residue depth, and CAB depth.

        - **af_missense_cutoff (_type_)**:<br />
            Cutoff value for AlphaMissense pathogenicity prediction.
        """

        for (
            mutation,
            gene,
            disease,
            clinv_p,
            afm_s
        ) in zip(
            self.clinvar_mutations,
            self.genes,
            self.disease_assoc,
            self.clinvar_pathogenicity,
            self.afmissense_scores
        ):
            res_num = split_missense_mutation(
                p_mutation=mutation,
                return_type="res_num"
            )

            if self.skip_mutation(
                entity_=entity_,
                gene=gene,
                clinv_p=clinv_p,
                afm_s=afm_s,
                modeled_range=modeled_range,
                ch_id=ch_id,
                res_num=res_num,
                af_missense_cutoff=af_missense_cutoff,
            ):
                continue

            self.enrich_variant_annotations(
                rigid_body_ranges=rigid_body_ranges,
                ch_id=ch_id,
                entity_=entity_,
                psa_map=psa_map,
                entity_ptm_df=entity_ptm_df,
                temp_structure_path=temp_structure_path,
                burial_info_dict=burial_info_dict,
                mutation=mutation,
                gene=gene,
                disease=disease,
            )

    def enrich_variant_annotations(
        self,
        rigid_body_ranges: dict,
        ch_id: str,
        entity_: str,
        psa_map: dict,
        entity_ptm_df: pd.DataFrame,
        temp_structure_path: str,
        burial_info_dict: dict,
        mutation: str,
        gene: str,
        disease: str,
    ):
        """ Enrich the annotations for a specific variant with structural information.

        The structural information includes:
        - affected interfaces
        - disease association
        - rigid body parent
        - secondary structure
        - relative solvent accessibility (RSA)
        - residue depth
        - CAB depth (residue depth calculated for representative CA or CB atom)

        ## Arguments:

        - **rigid_body_ranges (dict)**:<br />
            Dictionary mapping chain IDs to lists of residue numbers representing
            the rigid body ranges for each chain.

        - **ch_id (str)**:<br />
            Chain ID in the structure associated with the mutation.

        - **entity_ (str)**:<br />
            Protein entity name in the structure associated with the mutation.

        - **psa_map (dict)**:<br />
            Dictionary mapping residue numbers in the reference sequence to
            residue numbers in the modeled sequence, used for mapping PTM
            annotations.
            This is useful if the PTM is mapped on a different isoform than the
            one that is modeled.

        - **entity_ptm_df (pd.DataFrame)**:<br />
            DataFrame containing PTM annotations for the protein entity, including
            columns for position, UniProt ID, PTM type, residue, mutation site, and
            experimental verification status.

        - **temp_structure_path (str)**:<br />
            Path to a temporary PDB file containing the structure of the chain,
            used for calculating burial information such as secondary structure,
            relative solvent accessibility (RSA), residue depth, and CAB depth.

        - **burial_info_dict (dict)**:<br />
            Dictionary containing burial information for residues in the chain,
            indexed by (chain_id, residue_number) tuples, with values containing
            secondary structure, RSA value, residue depth, and CAB depth.

        - **mutation (str)**:<br />
            Mutation annotation in the format "R123C", where "R" is the wild-type
            residue, "123" is the residue number, and "C" is the mutant residue

        - **gene (str)**:<br />
            Gene name associated with the mutation.

        - **disease (str)**:<br />
            Disease association for the mutation, which may be a string or a list
            of diseases separated by newline characters.
        """

        wt, res_num, mut = split_missense_mutation(
            p_mutation=mutation,
            return_type="all"
        )
        _key_ = (gene, res_num)
        if _key_ not in self.enriched_variants:
            self.enriched_variants[_key_] = {
                "Gene": gene,
                "Residue Number": res_num,
                "Mutation": set(),
                "Affected Interfaces": set(),
                "Disease association": set(),
                "Rigid Body": set(),
                "Secondary Structure": set(),
                "RSA Value": set(),
                "Residue Depth": set(),
                "CAB Depth": set(),
                "PTM": set(),
            }
        ptm_annot = entity_ptm_df[entity_ptm_df["Position"] == res_num]
        ptm_types = ptm_annot["Type"].tolist()
        ptm_sites = ptm_annot["MutationSite"].tolist()
        ptm_verified = ptm_annot["Is_experimental_verification"].tolist()

        for ptm_type, ptm_site, is_verified in zip(ptm_types, ptm_sites, ptm_verified):
            if isinstance(ptm_site, str):
                ptm_site = eval(ptm_site) if ptm_site.startswith("[") else [ptm_site]
            else:
                ptm_site = []
            if is_verified == 1:
                ptm_type += " (E)"
            for _mut_ in ptm_site:
                _wt, _res_num, _mut = split_missense_mutation(
                    p_mutation=_mut_,
                    return_type="all"
                )
                _res_num_mapped, _ = PairwiseSequenceAlignment.get_mapped_residue(
                    psa_map=psa_map,
                    codon_number=_res_num,
                )
                if (
                    _res_num_mapped == res_num and
                    AA_MAP_ONE_TO_THREE.get(_mut) == mut and
                    AA_MAP_ONE_TO_THREE.get(_wt) == wt
                ):
                    self.enriched_variants[_key_]["PTM"].add(ptm_type)

        affected_interfaces = self.get_affected_interfaces(
            gene=gene,
            entity=entity_,
            res_num=res_num,
        )

        diseases = disease.split("\n") if isinstance(disease, str) else [disease]

        self.enriched_variants[_key_]["Affected Interfaces"].update(affected_interfaces)
        self.enriched_variants[_key_]["Disease association"].update(diseases)
        self.enriched_variants[_key_]["Mutation"].add(mutation)

        if res_num in rigid_body_ranges.get(ch_id, []):
            self.enriched_variants[_key_]["Rigid Body"].add(Path(temp_structure_path).stem)

        burial_info = burial_info_dict.get((ch_id, res_num), None)
        if burial_info is not None:
            self.enriched_variants[_key_]["Secondary Structure"].add(burial_info["secondary_structure"])
            self.enriched_variants[_key_]["RSA Value"].add(str(burial_info["rsa_val"]))
            self.enriched_variants[_key_]["Residue Depth"].add(str(burial_info["residue_depth"]))
            self.enriched_variants[_key_]["CAB Depth"].add(str(burial_info["residue_cab_depth"]))

    def prepare_ptm_df(
        self,
        entity_: str,
        uniprot_base: str,
    ):
        """ Prepare PTM annotations for a given protein entity from the fetched
        PTMD database.

        ## Arguments:

        - **entity_ (str)**:<br />
            Protein entity name in the structure associated with the mutation.

        - **uniprot_base (str)**:<br />
            Base UniProt ID for the protein entity, used to fetch PTM annotations
            from the PTMD database.

        ## Returns:

        - **psa_map (dict)**:<br />
            A dictionary mapping positions in the modeled sequence to positions in the reference sequence.

        - **entity_ptm_df (pd.DataFrame)**:<br />
            A DataFrame containing PTM annotations for the specified entity.
        """

        modeled_seq = self.protein_sequences.get(self.protein_uniprot_map[entity_], None)
            # we'll reuse the af_missense sequences
        fasta_dict = fetch_fasta_dict_for_af_missense(
                Path(self.alpha_missense_dir) / "af_missense_sequences.fasta",
                self.protein_uniprot_map,
                overwrite=False,
            )
        ref_seq = fasta_dict.get(uniprot_base, None)
        _pairwise_alignment = PairwiseSequenceAlignment(
                seq1=modeled_seq,
                seq2=ref_seq,
            )
        psa_map = _pairwise_alignment.fetch_pairwise_alingment_map(
                pairwise_alignment_file=os.path.join(
                    self.pairwise_alignments_dir,
                    f"{entity_}_ptm_vs_modeled.fasta"),
                overwrite=True,
            )

        entity_ptm_df = pd.DataFrame(
                data=[],
                columns=[
                    "Position", "UniProt", "Type",
                    "Residue", "MutationSite", "Is_experimental_verification"
                ]
            )
        if self.include_ptm:
            entity_ptm_df = self.ptmd_total_df[self.ptmd_total_df["UniProt"] == uniprot_base]
            entity_ptm_df["Position"] = entity_ptm_df["Position"].map(psa_map)

        return psa_map, entity_ptm_df

    def get_mutated_rigid_body_residues(
        self,
        structure: Bio.PDB.Structure.Structure,
        entity_chain_map: dict,
    ):
        """ Identify mutated residues that are part of rigid bodies in the structure.

        ## Arguments:

        - **structure (Bio.PDB.Structure.Structure)**:<br />
            Rigid body structure object.

        - **entity_chain_map (dict)**:<br />
            Mapping of chain IDs to protein entities in the structure.

        ## Returns:

        - **dict**:<br />
            A dictionary mapping chain IDs to sets of mutated residue numbers that
            are part of rigid bodies.
        """

        mutated_residues = defaultdict(set)
        model = structure[0]

        residues = [
            (chain.id, res_num.id[1])
            for chain in model.get_chains()
            for res_num in model[chain.id].get_residues()
        ]

        mutations = self.clinvar_mutations_df["Mutation"].tolist()
        if "Protein" in self.clinvar_mutations_df.columns:
            proteins = self.clinvar_mutations_df["Protein"].tolist()
        else:
            genes = self.clinvar_mutations_df["Gene"].tolist()
            proteins = [self.gene_protein_map.get(g, g) for g in genes]

        for p, m in zip(proteins, mutations):
            res_num = split_missense_mutation(
                p_mutation=m,
                return_type="res_num"
            )
            ch_ids = [ch for ch, e in entity_chain_map.items() if e == p]
            for ch_id in ch_ids:
                if (ch_id, res_num) in residues:
                    mutated_residues[ch_id].add(res_num)

        return mutated_residues

def modify_affected_interface_val(row):
    if (
        row["RSA Value"] <= ProteinCoreCriteria.rsa_threshold and
        row["Residue Depth"] >= ProteinCoreCriteria.residue_depth_threshold
    ):
        return "core"
    else:
        return row["Affected Interfaces"]
