from __future__ import annotations
import os
import warnings
import pandas as pd
from string import Template
from itertools import product
from IMP_Toolbox.utils_imp_toolbox.file_helpers import write_json
from IMP_Toolbox.analysis.cocomap.cocomap_constants import (
    REP_ATOMS,
    TEMPLATE_CONFIG,
    DOCKER_BASE_COMMAND,
    ATOM_COL_NAMES,
    VALID_INTERACTIONS,
)
from IMP_Toolbox.pre_processing.mutations.utils_mutation import (
    split_missense_mutation,
)

def run_cocomap_docker(
    processed_struct_path: str,
    docker_result_dir: str,
    result_head: str,
    result_metadata: dict,
    docker_base_command: Template = Template(DOCKER_BASE_COMMAND),
    dry_run: bool = False,
    overwrite: bool = False,
):
    """ Run COCOMAP analysis using its docker image.

    Args:

        processed_struct_path (str):
            Path to processed structure file (PDB) to analyze.

        docker_result_dir (str):
            Directory where COCOMAP results will be stored.

        result_head (str):
            Unique identifier for the result (used for naming output directory).

        result_metadata (dict):
            Metadata dictionary containing information about the result,

        docker_base_command (Template, optional):
            Docker command template to run COCOMAP.

        dry_run (bool, optional):
            If True, only print the docker command without executing it.

        overwrite (bool, optional):
            If True, overwrite existing results.
    """

    cocomap_output_dir = os.path.join(docker_result_dir, result_head)
    cocomap_output_dir = os.path.abspath(cocomap_output_dir)

    if os.path.exists(cocomap_output_dir) and overwrite is False:
        warnings.warn(
            f"""COCOMAP results already exist for {result_head} at
            {cocomap_output_dir}. Skipping.
            To overwrite, set overwrite=True."""
        )
        return  # already processed

    processed_struct_path = os.path.abspath(processed_struct_path)
    if not os.path.isfile(processed_struct_path):
        raise Exception(
            f"Processed structure file not found for \
            {result_head}: {processed_struct_path}."
        )

    struct_type = os.path.splitext(processed_struct_path)[1].lower()

    os.makedirs(cocomap_output_dir, exist_ok=True)

    # copy processed structure to the result directory
    target_pdb_path = os.path.join(cocomap_output_dir, f"{result_head}{struct_type}")
    os.system(f"cp {processed_struct_path} {target_pdb_path}")

    docker_command = docker_base_command.substitute(
        path_to_mount=cocomap_output_dir,
        config_path=os.path.join(cocomap_output_dir, "config.json")
    )

    # write config file required by COCOMAP to the result directory
    config_params = TEMPLATE_CONFIG.copy()
    config_params.update(
        {
            "chains_set_1": result_metadata["chains_set_1"],
            "chains_set_2": result_metadata["chains_set_2"],
            "ranges_1": result_metadata["ranges_1"],
            "ranges_2": result_metadata["ranges_2"],
            "pdb_file": target_pdb_path,
        }
    )
    config_path = os.path.join(cocomap_output_dir, "config.json")
    write_json(config_path, config_params)

    print(f"Running Docker command for {result_head}:\n{docker_command}\n")
    if not dry_run:
        os.system(docker_command)
        os.system(f"sudo chown -R $USER:$USER {cocomap_output_dir}")
    print("=" * 100 + "\n")

def postprocess_cocomap_results(
    result_metadata: dict,
    docker_results_dir: str,
    result_head: str,
):
    """ Rename COCOMAP CSVs columns by protein names & explode by interaction types.

    Args:

        result_metadata (dict):
            Metadata dictionary containing information about the result,
            includes 'proteins' key with list of protein names for chain1 & 2.

        docker_results_dir (str):
            Directory where COCOMAP results are stored.

        result_head (str):
            Unique identifier for the result (used for naming output directory).
    """

    cocomap_output_dir = os.path.join(docker_results_dir, result_head)
    cocomap_output_dir = os.path.abspath(cocomap_output_dir)

    proteins = result_metadata.get("proteins", [])
    assert len(proteins) in [0, 2], "proteins must 0 or 2"
    prot_1, prot_2 = proteins[0], proteins[1]
    new_cols = {
        "Chain 1": prot_1,
        "Chain 2": prot_2,
    }

    csv_files = [
        f for f in os.listdir(cocomap_output_dir) if f.endswith(".csv")
    ]
    # Rename chain columns to protein names
    for csv_file in csv_files:
        csv_path = os.path.join(cocomap_output_dir, csv_file)
        df = pd.read_csv(csv_path, index_col=0, encoding='utf-8')
        df = df.rename(columns=new_cols)
        df.to_csv(csv_path)

    final_csvs = [
        f for f in os.listdir(cocomap_output_dir)
        if f.endswith("final_file.csv")
    ]

    # Explode final CSVs by interaction types to keep one interaction per row
    for final_csv in final_csvs:
        final_csv_path = os.path.join(cocomap_output_dir, final_csv)
        df = pd.read_csv(final_csv_path, index_col=0, encoding='utf-8')
        df["Type of Interactions"] = df["Type of Interactions"].str.split("; ")
        df_exploded = df.explode("Type of Interactions")
        df_exploded = df_exploded[
            df_exploded["Type of Interactions"].str.strip() != ""
        ]
        df_exploded = df_exploded.drop_duplicates()
        df_exploded.to_csv(final_csv_path, index=False)

def get_mutation_annotated_df(
    cocomap_df: pd.DataFrame,
    mutation1_df: pd.DataFrame,
    mutation2_df: pd.DataFrame,
) -> pd.DataFrame:
    """ Add mutation annotations to COCOMAP dataframe.

    mutations1_df, mutations2_df are dataframes containing mutation information
    in the output format of clinvar_mutations1.py

    NOTE: This is temporary, might be changed later.

    Args:
        cocomap_df (pd.DataFrame): _description_
        mutation1_df (pd.DataFrame): _description_
        mutation2_df (pd.DataFrame): _description_

    Returns:
        pd.DataFrame: _description_
    """

    df_res_1 = list(cocomap_df["Res. Number 1"].unique())
    df_res_2 = list(cocomap_df["Res. Number 2"].unique())

    mutations1 = mutation1_df["Mutation"].tolist()
    mutations2 = mutation2_df["Mutation"].tolist()

    af_missense1_score = mutation1_df["AlphaMissense score"].tolist()
    af_missense2_score = mutation2_df["AlphaMissense score"].tolist()

    af_missense1_patho = mutation1_df["AlphaMissense pathogenicity"].tolist()
    af_missense2_patho = mutation2_df["AlphaMissense pathogenicity"].tolist()

    clinvar1_patho = mutation1_df["ClinVar clinical significance"].tolist()
    clinvar2_patho = mutation2_df["ClinVar clinical significance"].tolist()

    mutated_res_1 = {}
    mutated_res_2 = {}
    afm_score_res_1 = {}
    afm_score_res_2 = {}

    for res1, res2 in zip(df_res_1, df_res_2):
        mutated_res_1[res1] = [
            mut for idx, mut in enumerate(mutations1)
            if str(res1) == split_missense_mutation(mut)[1]
            and af_missense1_patho[idx] == "Likely pathogenic"
            and clinvar1_patho[idx] in ["Likely pathogenic", "Pathogenic", "Conflicting classifications of pathogenicity", "Pathogenic/Likely pathogenic"]
        ]
        afm_score_res_1[res1] = [
            af_missense1_score[idx] for idx, mut in enumerate(mutations1)
            if str(res1) == split_missense_mutation(mut)[1]
            and af_missense1_patho[idx] == "Likely pathogenic"
        ]

        mutated_res_2[res2] = [
            mut for idx, mut in enumerate(mutations2)
            if str(res2) == split_missense_mutation(mut)[1]
            and af_missense2_patho[idx] == "Likely pathogenic"
            and clinvar2_patho[idx] in ["Likely pathogenic", "Pathogenic", "Conflicting classifications of pathogenicity", "Pathogenic/Likely pathogenic"]
        ]
        afm_score_res_2[res2] = [
            af_missense2_score[idx] for idx, mut in enumerate(mutations2)
            if str(res2) == split_missense_mutation(mut)[1]
            and af_missense2_patho[idx] == "Likely pathogenic"
        ]

    cocomap_df["Mutated 1"] = cocomap_df["Res. Number 1"].map(mutated_res_1)
    cocomap_df["Mutated 2"] = cocomap_df["Res. Number 2"].map(mutated_res_2)
    cocomap_df["AF Missense 1 Score"] = cocomap_df["Res. Number 1"].map(afm_score_res_1)
    cocomap_df["AF Missense 2 Score"] = cocomap_df["Res. Number 2"].map(afm_score_res_2)

    # convert lists to strings for better readability
    cocomap_df["Mutated 1"] = cocomap_df["Mutated 1"].apply(
        lambda x: ", ".join(x) if isinstance(x, list) and len(x) > 0 else ""
    )
    cocomap_df["Mutated 2"] = cocomap_df["Mutated 2"].apply(
        lambda x: ", ".join(x) if isinstance(x, list) and len(x) > 0 else ""
    )
    cocomap_df["AF Missense 1 Score"] = cocomap_df["AF Missense 1 Score"].apply(
        lambda x: ", ".join(map(str, x)) if isinstance(x, list) and len(x) > 0 else ""
    )
    cocomap_df["AF Missense 2 Score"] = cocomap_df["AF Missense 2 Score"].apply(
        lambda x: ", ".join(map(str, x)) if isinstance(x, list) and len(x) > 0 else ""
    )

    return cocomap_df

def add_af_metrics(
    result_metadata: dict,
    cocomap_output_dir: str,
    af_metadata: dict,
    rep_atom_dict: dict = {},
    average_token_pae: bool = False,
    average_token_plddt: bool = False,
    metric_level: str = "per_token",
    mark_clashes: bool = False,
):
    """ Add and save COCOMAP results with AlphaFold metrics (PAE, pLDDT) columns.

    Args:

        result_metadata (dict):
            Metadata dictionary containing information about the result.

        cocomap_output_dir (str):
            Directory where COCOMAP results are stored.

        af_metadata (dict):
            Metadata dictionary containing information about the AlphaFold result,
            includes 'structure_path' and 'data_path' keys.
            Optionally includes 'af_offset' key for residue numbering offsets.

        rep_atom_dict (dict, optional):
            Dictionary mapping residue names to representative atom names.
            Optional, only need if non-standard representative atoms are used.

        average_token_pae (bool, optional):
            If True, use average PAE over all tokens for a residue.

        average_token_plddt (bool, optional):
            If True, use average pLDDT over all tokens for a residue.

        metric_level (str, optional):
            Level at which to compute metrics: 'per_token' or 'representative_token'.

        mark_clashes (bool, optional):
            If True, mark interactions that are clashes as defined by COCOMAPS2.
            COCOMAPS2 annotates interaction types with "*" if the distance between
            the interacting atoms is less than the sum of their Van der Waals radii.
            BUG: Currently, not matching to the final_file.csv clashes

    """

    try:
        from af_pipeline._initialize import _Initialize

    except ImportError:
        raise ImportError(
            "af_pipeline module is required for adding AF metrics. "
        )

    af_struct_path = af_metadata["structure_path"]
    af_data_path = af_metadata["data_path"]
    af_offset = af_metadata.get("af_offset", {})

    _default_chains = ["Chain 1", "Chain 2"]
    chains_set_1_prot = result_metadata.get("proteins", _default_chains)[0]
    chains_set_2_prot = result_metadata.get("proteins", _default_chains)[1]

    grouped_csvs_by_chain_pairs = group_csv_files_by_chain_pairs(
        result_metadata,
        cocomap_output_dir
    )

    processed_struct_path = result_metadata["structure_path"]
    struct_type = os.path.splitext(processed_struct_path)[1].lower()[1:]

    initializer = _Initialize( #TODO: make public in af_pipeline
        data_file_path=af_data_path,
        structure_file_path=af_struct_path,
        af_offset=af_offset,
        rep_atom_dict=rep_atom_dict,
        average_token_pae=average_token_pae,
        average_token_plddt=average_token_plddt,
        metric_level=metric_level,
    )

    for chain_pair, csv_files in grouped_csvs_by_chain_pairs.items():

        chain_1, chain_2 = chain_pair

        final_csvs = [f for f in csv_files if f.endswith("final_file.csv")]
        assert len(final_csvs) <= 1, (
            f"""Multiple final CSVs found for chain pair {chain_1},{chain_2}"""
        )

        master_df_list = []
        for csv_f in csv_files:
            print("-"*100)
            csv_identifier = (
                csv_f.split(".")[2]
                .replace(f"{struct_type}_{chain_1}_{chain_2}_", "")
                .replace(f"{struct_type}_{chain_2}_{chain_1}_", "")
            )

            if csv_identifier not in VALID_INTERACTIONS:
                print(f"Skipping {csv_f} as it is not a recognized contact type.")
                continue

            print(f"Processing {csv_f}.")

            csv_path = os.path.join(cocomap_output_dir, csv_f)
            df = pd.read_csv(csv_path, index_col=0, encoding='utf-8')
            df_chains_1 = list(df[chains_set_1_prot])
            df_chains_2 = list(df[chains_set_2_prot])
            df_res_1 = list(df["Res. Number 1"])
            df_res_2 = list(df["Res. Number 2"])
            df_res_names_1 = list(df["Res. Name 1"])
            df_res_names_2 = list(df["Res. Name 2"])

            interaction_type = VALID_INTERACTIONS[csv_identifier]
            df_interaction_types = [interaction_type]*len(df)

            if mark_clashes:
                df_interaction_types = mark_clashes_in_interaction_types(
                    df=df,
                    interaction_type=interaction_type,
                )

            df_atom_1, df_atom_2 = get_interacting_atoms(df, csv_identifier)

            zipped_cra = zip(
                df_chains_1, df_res_1, df_atom_1,
                df_chains_2, df_res_2, df_atom_2,
            )

            af3_metrics = extract_af3_metrics(
                rep_atom_dict,
                initializer,
                zipped_cra
            )

            atom1_plddts = af3_metrics["atom1_plddts"]
            atom2_plddts = af3_metrics["atom2_plddts"]
            atom1_names = af3_metrics["atom1_names"]
            atom2_names = af3_metrics["atom2_names"]
            metric_atom1_names = af3_metrics["metric_atom1_names"]
            metric_atom2_names = af3_metrics["metric_atom2_names"]
            pae_values = af3_metrics["pae_values"]
            contact_prob_values = af3_metrics["contact_prob_values"]

            assert len(pae_values) == len(atom1_plddts) == len(atom2_plddts), (
                "Length mismatch in computed AF metrics."
            )

            df["AF PAE"] = pae_values
            df["Atom 1 pLDDT"] = atom1_plddts
            df["Atom 2 pLDDT"] = atom2_plddts

            if len(contact_prob_values) == len(pae_values):
                df["AF Contact Probability"] = contact_prob_values

            new_csv_path = os.path.join(
                cocomap_output_dir,
                "af_metrics",
                csv_f.replace(".csv", "_with_AF_metrics.csv")
            )
            os.makedirs(os.path.dirname(new_csv_path), exist_ok=True)
            df.to_csv(new_csv_path, index=False, encoding='utf-8')

            df_dict = {
                "Res. Name 1": df_res_names_1,
                "Res. Number 1": df_res_1,
                chains_set_1_prot: df_chains_1,
                "Res. Name 2": df_res_names_2,
                "Res. Number 2": df_res_2,
                chains_set_2_prot: df_chains_2,
                "Type of Interactions": df_interaction_types,
                "AF PAE": pae_values,
                "Atom 1 pLDDT": atom1_plddts,
                "Atom 2 pLDDT": atom2_plddts,
                "Atom 1": atom1_names,
                "Atom 2": atom2_names,
                "Metric Atom 1": metric_atom1_names,
                "Metric Atom 2": metric_atom2_names,
            }

            if len(contact_prob_values) == len(pae_values):
                df_dict["AF Contact Probability"] = contact_prob_values

            df_to_concat = pd.DataFrame(df_dict)
            master_df_list.append(df_to_concat)

        af_final_csv = pd.concat(master_df_list, ignore_index=True)
        af_final_csv_path = os.path.join(
            cocomap_output_dir,
            "af_metrics",
            final_csvs[0].replace(".csv", "_with_AF_metrics.csv")
        )
        os.makedirs(os.path.dirname(af_final_csv_path), exist_ok=True)
        af_final_csv.to_excel(af_final_csv_path.replace(".csv", ".xlsx"), index=False)

def extract_af3_metrics(
    rep_atom_dict: dict,
    initializer: af_pipeline._initialize._Initialize,
    zipped_cra: zip,
):
    """ Extract AlphaFold 3 metrics for interacting atoms.

    Args:

        rep_atom_dict (dict):
            Dictionary mapping residue names to representative atom names.

        initializer (af_pipeline._initialize._Initialize):
            Initialized AlphaFold parser object.

        zipped_cra (zip):
            Zipped chain, residue, atom information for interacting pairs.

    Returns:
        dict:
            Dictionary containing extracted AF3 metrics.
    """

    atom1_idxs = []
    atom2_idxs = []
    atom1_plddts = []
    atom2_plddts = []
    atom1_names = []
    atom2_names = []
    metric_atom1_names = []
    metric_atom2_names = []

    for (chain_1, res_1, atom_1, chain_2, res_2, atom_2) in zipped_cra:
        atom1_names.append(atom_1)
        atom2_names.append(atom_2)

        atom1_idx = get_atom_idx(
            num_to_idx=initializer.num_to_idx,
            chain=chain_1,
            res=res_1,
            atom=atom_1,
        )

        atom2_idx = get_atom_idx(
            num_to_idx=initializer.num_to_idx,
            chain=chain_2,
            res=res_2,
            atom=atom_2,
        )

        atom1_idxs.append(atom1_idx)
        atom2_idxs.append(atom2_idx)

        atom1_plddt, metric_atom1 = get_atom_plddt(
            rep_atom_dict,
            initializer,
            chain_1,
            res_1,
            atom_1
        )

        atom2_plddt, metric_atom2 = get_atom_plddt(
            rep_atom_dict,
            initializer,
            chain_2,
            res_2,
            atom_2
        )

        metric_atom1_names.append(metric_atom1)
        metric_atom2_names.append(metric_atom2)

        atom1_plddts.append(atom1_plddt)
        atom2_plddts.append(atom2_plddt)

    pae_values = []
    for atoms1, atoms2 in zip(atom1_idxs, atom2_idxs):
        pae_submatrix = initializer.avg_pae[atoms1, :][:, atoms2]
        pae_values.append(pae_submatrix.mean().item())

    contact_prob_values = []
    if initializer.contact_probs is not None:
        for atoms1, atoms2 in zip(atom1_idxs, atom2_idxs):
            _ij = initializer.contact_probs[atoms1, :][:, atoms2]
            _ji = initializer.contact_probs[atoms2, :][:, atoms1]
            contact_prob_value = (_ij.mean().item() + _ji.mean().item()) / 2
            contact_prob_values.append(contact_prob_value)

    af3_metrics = {
        "atom1_plddts": atom1_plddts,
        "atom2_plddts": atom2_plddts,
        "atom1_names": atom1_names,
        "atom2_names": atom2_names,
        "metric_atom1_names": metric_atom1_names,
        "metric_atom2_names": metric_atom2_names,
        "pae_values": pae_values,
        "contact_prob_values": contact_prob_values,
    }

    return af3_metrics

def get_atom_plddt(
    rep_atom_dict: dict,
    initializer: af_pipeline._initialize._Initialize,
    chain: str,
    res: int,
    atom: str,
):
    """ Get pLDDT score for a specific atom.

    Args:

        rep_atom_dict (dict):
            Dictionary mapping residue names to representative atom names.

        initializer (af_pipeline._initialize._Initialize):
            Initialized AlphaFold parser object.

        chain (str):
            Chain identifier.

        res (int):
            Residue number.

        atom (str):
            Atom name.

    Returns:
        tuple:
            pLDDT score and metric atom name.
    """

    try:
        from af_pipeline.parser.structure_parser import StructureParser
        import Bio.PDB
        import Bio.PDB.Atom

    except ImportError:
        raise ImportError(
            "af_pipeline module is required for adding AF metrics. "
        )

    metric_dict = {
        str: lambda x: x,
        Bio.PDB.Atom.Atom: lambda x: x.get_name(),
    }

    orig_res = initializer.renumber.original_chain_res_num(
        chain_res_num=res,
        chain_id=chain,
    )
    res_obj = initializer.structure[0][chain][orig_res]

    if "-" in atom: # Ring From column has residue as (PHE-125)
        res1_name = atom.split("-")[0]
        assert res1_name == res_obj.get_resname(), (
            f"Residue name mismatch: {res1_name} vs {res_obj.get_resname()}"
        )

        rep_atom = rep_atom_dict.get(
            res_obj.get_resname(),
            StructureParser.get_rep_atom(residue=res_obj)
        )

        metric_atom1 = metric_dict[type(rep_atom)](rep_atom)

    elif atom in REP_ATOMS:
        rep_atom = rep_atom_dict.get(
            res_obj.get_resname(),
            REP_ATOMS.get(atom, StructureParser.get_rep_atom(residue=res_obj))
        )

        metric_atom1 = metric_dict[type(rep_atom)](rep_atom)

    else:
        metric_atom1 = atom

    quants = StructureParser.extract_perresidue_quantities(
        residue=res_obj,
        quantities=["plddt"],
        rep_atom=metric_atom1,
    )
    atom1_plddt = quants["plddt"]

    return atom1_plddt, metric_atom1

def get_atom_idx(
    num_to_idx: dict,
    chain: str,
    res: int,
    atom: str,
) -> list:
    """ Get atom index/indices from num_to_idx mapping.

    Args:

        num_to_idx (dict):
            Residue number to residue index mapping.

        chain (str):
            Chain identifier.

        res (int):
            Residue number.

        atom (str):
            Atom name.

    Returns:
        list:
            List of atom indices corresponding to the specified atom.
    """

    assert chain in num_to_idx, f"Chain {chain} not found in num_to_idx."
    assert res in num_to_idx[chain], (
        f"Residue {res} not found in chain {chain} in num_to_idx."
    )

    if atom in num_to_idx[chain][res]:
        atom_idx = [num_to_idx[chain][res][atom]]

    else:
        atom_idx = list(num_to_idx[chain][res].values())

    return atom_idx

def get_interacting_atoms(
    df: pd.DataFrame,
    csv_identifier: str,
):
    """ Get interacting atoms from COCOMAP DataFrame based on CSV identifier.

    Args:

        df (pd.DataFrame):
            DataFrame containing COCOMAP interaction data.

        csv_identifier (str):
            Identifier for the type of interaction (used to determine atom columns).

    Returns:
        tuple:
            Two lists containing interacting atoms for the two residues.
    """

    atom_col_names = ATOM_COL_NAMES.get(csv_identifier, [])

    if len(atom_col_names) == 2:
        atom1_col_name = atom_col_names[0]
        atom2_col_name = atom_col_names[1]
        df_atom_1 = list(df[atom1_col_name])
        df_atom_2 = list(df[atom2_col_name])

    elif len(atom_col_names) == 3: # for pi interactions
        atom_name_col = atom_col_names[0]
        atom_from_col = atom_col_names[1] # residue of the atom (LYS)
        ring_from_col = atom_col_names[2] # ring residue (PHE-125)
        df_atom_1 = []
        df_atom_2 = []

        for _, row in df.iterrows():
            _res1 = f'{row["Res. Name 1"]}-{row["Res. Number 1"]}'
            _res2 = f'{row["Res. Name 2"]}-{row["Res. Number 2"]}'
            atom_from = row[atom_from_col]
            ring_from = row[ring_from_col]

            if atom_from in _res1 and ring_from == _res2:
                df_atom_1.append(row[atom_name_col])
                df_atom_2.append(ring_from)

            elif atom_from in _res2 and ring_from == _res1:
                df_atom_1.append(ring_from)
                df_atom_2.append(row[atom_name_col])

    else: # for pi-pi interactions fall back to residue level
        df_atom_1 = list(df["Res. Name 1"])
        df_atom_2 = list(df["Res. Name 2"])

    return df_atom_1, df_atom_2

def mark_clashes_in_interaction_types(
    df: pd.DataFrame,
    interaction_type: str,
):
    """ Mark clashes in interaction types based on distance column.

    Args:

        df (pd.DataFrame):
            DataFrame containing COCOMAP interaction data.

        df_interaction_types (list):
            List of interaction types corresponding to each row in the DataFrame.

    Returns:
        list:
            List of interaction types with clashes marked by '*'.
    """

    df_interaction_types_w_clash = [] # to mark clashes with '*'
    df_cols = df.columns.tolist()
    df_interaction_types = [interaction_type]*len(df)

    if "Distance (Å)" in df_cols:
        df_distances = list(df["Distance (Å)"])

        for dist_, int_type in zip(df_distances, df_interaction_types):

            if "*" in str(dist_):
                df_interaction_types_w_clash.append(int_type + "*")

            else:
                df_interaction_types_w_clash.append(int_type)

    else:
        df_interaction_types_w_clash = df_interaction_types

    return df_interaction_types_w_clash

def group_csv_files_by_chain_pairs(
    result_metadata: dict,
    cocomap_output_dir: str,
):
    """ Group COCOMAP CSV files by chain pairs.

    Args:

        result_metadata (dict):
            Metadata dictionary containing information about the result.

        cocomap_output_dir (str):
            Directory where COCOMAP results are stored.

    Returns:
        dict:
            Dictionary mapping chain pairs to lists of corresponding CSV files.
    """

    chains_set_1 = result_metadata["chains_set_1"]
    chains_set_2 = result_metadata["chains_set_2"]
    chain_pair_combinations = list(product(chains_set_1, chains_set_2))

    csv_files = [f for f in os.listdir(cocomap_output_dir) if f.endswith(".csv")]

    grouped_csvs_by_chain_pairs = {}

    for chain_pair in chain_pair_combinations:

        ch_1, ch_2 = chain_pair
        grouped_csvs_by_chain_pairs[chain_pair] = []

        for csv_f in csv_files:
            if f"_{ch_1}_{ch_2}_" in csv_f or f"_{ch_2}_{ch_1}_" in csv_f:
                grouped_csvs_by_chain_pairs[chain_pair].append(csv_f)

    grouped_csvs_by_chain_pairs = {
        k: v for k, v in grouped_csvs_by_chain_pairs.items() if len(v) > 0
    }

    return grouped_csvs_by_chain_pairs
