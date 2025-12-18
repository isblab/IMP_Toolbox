from math import prod
import os
import warnings
import pandas as pd
from string import Template
from IMP_Toolbox.utils_imp_toolbox.file_helpers import write_json
from IMP_Toolbox.analysis.cocomap.cocomap_constants import (
    TEMPLATE_CONFIG,
    DOCKER_BASE_COMMAND,
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

    cocomap_output_dir = os.path.join(docker_result_dir, result_head)
    cocomap_output_dir = os.path.abspath(cocomap_output_dir)
    if os.path.exists(cocomap_output_dir) and overwrite is False:
        warnings.warn(
            f"""COCOMAP results already exist for {result_head} at \
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

    os.makedirs(cocomap_output_dir, exist_ok=True)
    # copy processed structure to mount path
    target_path = os.path.join(cocomap_output_dir, f"{result_head}.pdb")
    os.system(f"cp {processed_struct_path} {target_path}")

    docker_command = docker_base_command.substitute(
        path_to_mount=cocomap_output_dir,
        config_path=os.path.join(cocomap_output_dir, "config.json")
    )

    chains_set_1 = result_metadata["chains_set_1"]
    chains_set_2 = result_metadata["chains_set_2"]
    ranges_1 = result_metadata["ranges_1"]
    ranges_2 = result_metadata["ranges_2"]
    pdb_file = target_path

    config_params = TEMPLATE_CONFIG.copy()
    config_params.update(
        {
            "chains_set_1": chains_set_1,
            "chains_set_2": chains_set_2,
            "ranges_1": ranges_1,
            "ranges_2": ranges_2,
            "pdb_file": pdb_file,
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
    for csv_file in csv_files:
        csv_path = os.path.join(cocomap_output_dir, csv_file)
        df = pd.read_csv(csv_path, index_col=0, encoding='utf-8')
        df = df.rename(columns=new_cols)
        df.to_csv(csv_path)

    final_csvs = [
        f for f in os.listdir(cocomap_output_dir)
        if f.endswith("final_file.csv")
    ]

    for final_csv in final_csvs:
        final_csv_path = os.path.join(cocomap_output_dir, final_csv)
        df = pd.read_csv(final_csv_path, index_col=0, encoding='utf-8')
        df["Type of Interactions"] = df["Type of Interactions"].str.split("; ")
        df_exploded = df.explode("Type of Interactions")
        df_exploded = df_exploded[
            df_exploded["Type of Interactions"].str.strip() != ""
        ]
        # df_exploded = df_exploded.drop_duplicates()
        df_exploded.to_csv(final_csv_path, index=False)


def add_af_metrics(
    result_metadata: dict,
    cocomap_output_dir: str,
    af_struct_path: str,
    af_data_path: str,
    af_offset: dict = {},
    rep_atom_dict: dict = {},
    average_token_pae: bool = False,
    average_token_plddt: bool = False,
    metric_level: str = "per_token",
):

    chains_set_1 = result_metadata["chains_set_1"]
    chains_set_2 = result_metadata["chains_set_2"]
    chains_set_1_prot = result_metadata.get("proteins", ["Chain 1", "Chain 2"])[0]
    chains_set_2_prot = result_metadata.get("proteins", ["Chain 1", "Chain 2"])[1]

    from itertools import product

    chain_pair_combinations = list(product(chains_set_1, chains_set_2))

    csv_files = [
        f for f in os.listdir(cocomap_output_dir)
        if f.endswith(".csv")
    ]

    grouped_csvs_by_chain_pairs = {}
    for chain_pair in chain_pair_combinations:
        chain_1, chain_2 = chain_pair
        grouped_csvs_by_chain_pairs[chain_pair] = []
        for csv_f in csv_files:
            if (
                f"_{chain_1}_{chain_2}_" in csv_f
                or f"_{chain_2}_{chain_1}_" in csv_f
            ):
                grouped_csvs_by_chain_pairs[chain_pair].append(csv_f)

    grouped_csvs_by_chain_pairs = {
        k: v for k, v in grouped_csvs_by_chain_pairs.items() if len(v) > 0
    }

    from af_pipeline._initialize import _Initialize
    from af_pipeline.parser.structure_parser import StructureParser
    from IMP_Toolbox.analysis.cocomap.cocomap_constants import (
        ATOM_COL_NAMES,
        VALID_INTERACTIONS,
    )
    import Bio.PDB
    import Bio.PDB.Atom

    initializer = _Initialize(
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
        final_csvs = [
            f for f in csv_files if f.endswith("final_file.csv")
        ]
        assert len(final_csvs) <= 1, f"Multiple final CSVs found for chain pair {chain_1}, {chain_2}: {final_csvs}"
        print(final_csvs)

        final_csv = pd.read_csv(os.path.join(cocomap_output_dir, final_csvs[0]), index_col=0, encoding='utf-8')

        af_final_csv = pd.DataFrame(columns=final_csv.columns.tolist() + ["AF_PAE", "Atom 1 pLDDT", "Atom 2 pLDDT"])

        master_rows = []
        for csv_f in csv_files:
            print("-"*100)
            csv_identifier = csv_f.split(".")[2].replace(f"pdb_{chain_1}_{chain_2}_", "").replace(f"pdb_{chain_2}_{chain_1}_", "")
            if csv_identifier not in VALID_INTERACTIONS:
                print(f"Skipping {csv_f} as it is not a recognized contact type.")
                continue
            print(f"Processing {csv_f} for chain pair {chain_1}, {chain_2} and contact type {csv_identifier}.")
            interaction_type = VALID_INTERACTIONS[csv_identifier]
            print(interaction_type)

            csv_path = os.path.join(cocomap_output_dir, csv_f)
            df = pd.read_csv(csv_path, index_col=0, encoding='utf-8')
            df_cols = df.columns.tolist()
            df_chains_1 = list(df[chains_set_1_prot])
            df_chains_2 = list(df[chains_set_2_prot])
            df_res_1 = list(df["Res. Number 1"])
            df_res_2 = list(df["Res. Number 2"])
            df_res_names_1 = list(df["Res. Name 1"])
            df_res_names_2 = list(df["Res. Name 2"])
            df_type_of_interactions = [interaction_type]*len(df)

            df_type_of_interactions_w_clash = []

            if "Distance (Å)" in df_cols:
                df_distances = list(df["Distance (Å)"])
                for dist_, int_type in zip(df_distances, df_type_of_interactions):
                    if "*" in str(dist_):
                        df_type_of_interactions_w_clash.append(int_type + "*")
                    else:
                        df_type_of_interactions_w_clash.append(int_type)

            else:
                df_type_of_interactions_w_clash = df_type_of_interactions

            atom_col_names = ATOM_COL_NAMES.get(csv_identifier, [])

            if len(atom_col_names) == 2:
                atom1_col_name = atom_col_names[0]
                atom2_col_name = atom_col_names[1]
                df_atom_1 = list(df[atom1_col_name]) if atom1_col_name in df_cols else [None]*len(df)
                df_atom_2 = list(df[atom2_col_name]) if atom2_col_name in df_cols else [None]*len(df)

            elif len(atom_col_names) == 3:
                atom_from_col = atom_col_names[1]
                ring_from_col = atom_col_names[2]
                atom_name_col = atom_col_names[0]
                df_atom_1 = []
                df_atom_2 = []
                for idx, row in df.iterrows():
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
            else:
                # fallback to residue level
                df_atom_1 = list(df["Res. Name 1"])
                df_atom_2 = list(df["Res. Name 2"])

            zipped_cra = zip(
                df_chains_1,
                df_res_1,
                df_atom_1,
                df_chains_2,
                df_res_2,
                df_atom_2,
            )

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

                try:
                    atom1_idx = [initializer.num_to_idx[chain_1][res_1][atom_1]]
                except KeyError:
                    atom1_idx = list(initializer.num_to_idx[chain_1][res_1].values())

                try:
                    atom2_idx = [initializer.num_to_idx[chain_2][res_2][atom_2]]
                except KeyError:
                    atom2_idx = list(initializer.num_to_idx[chain_2][res_2].values())

                atom1_idxs.append(atom1_idx)
                atom2_idxs.append(atom2_idx)

                if len(atom1_idx) > 1:
                    atom1_plddt = initializer.token_plddts[atom1_idx]
                    metric_atom1_names.append(atom_1)
                else:
                    orig_res_1 = initializer.renumber.original_chain_res_num(
                        chain_res_num=res_1,
                        chain_id=chain_1,
                    )
                    res_obj = initializer.structure[0][chain_1][orig_res_1]

                    if "-" in atom_1:
                        res1_name = atom_1.split("-")[0]
                        assert res1_name == res_obj.get_resname(), f"Residue name mismatch: {res1_name} vs {res_obj.get_resname()}"

                        rep_atom = rep_atom_dict.get(
                            res_obj.get_resname(),
                            StructureParser.get_rep_atom(residue=res_obj)
                        )
                        if isinstance(rep_atom, str):
                            metric_atom1_names.append(rep_atom)
                        elif isinstance(rep_atom, Bio.PDB.Atom.Atom):
                            metric_atom1_names.append(rep_atom.get_name())
                        else:
                            raise Exception(f"Unexpected rep_atom type: {type(rep_atom)}")
                    else:
                        rep_atom = atom_1
                        metric_atom1_names.append(rep_atom)

                    quants = StructureParser.extract_perresidue_quantities(
                        residue=res_obj,
                        quantities=["plddt"],
                        rep_atom=rep_atom,
                    )
                    atom1_plddt = quants["plddt"]

                if len(atom2_idx) > 1:
                    atom2_plddt = initializer.token_plddts[atom2_idx]
                    metric_atom2_names.append(atom_2)
                else:
                    orig_res_2 = initializer.renumber.original_chain_res_num(
                        chain_res_num=res_2,
                        chain_id=chain_2,
                    )
                    res_obj = initializer.structure[0][chain_2][orig_res_2]

                    if "-" in atom_2:
                        res2_name = atom_2.split("-")[0]
                        assert res2_name == res_obj.get_resname(), f"Residue name mismatch: {res2_name} vs {res_obj.get_resname()}"

                        rep_atom = rep_atom_dict.get(
                            res_obj.get_resname(),
                            StructureParser.get_rep_atom(residue=res_obj)
                        )
                        if isinstance(rep_atom, str):
                            metric_atom2_names.append(rep_atom)
                        elif isinstance(rep_atom, Bio.PDB.Atom.Atom):
                            metric_atom2_names.append(rep_atom.get_name())
                        else:
                            raise Exception(f"Unexpected rep_atom type: {type(rep_atom)}")
                    else:
                        rep_atom = atom_2
                        metric_atom2_names.append(rep_atom)

                    quants = StructureParser.extract_perresidue_quantities(
                        residue=res_obj,
                        quantities=["plddt"],
                        rep_atom=rep_atom,
                    )
                    atom2_plddt = quants["plddt"]

                atom1_plddts.append(atom1_plddt)
                atom2_plddts.append(atom2_plddt)

            pae_values = []
            for atom1_idx_list, atom2_idx_list in zip(atom1_idxs, atom2_idxs):
                pae_submatrix = initializer.avg_pae[atom1_idx_list, :][:, atom2_idx_list]
                pae_values.append(pae_submatrix.mean().item())

            assert len(pae_values) == len(atom1_plddts) == len(atom2_plddts), "Length mismatch in computed AF metrics."
            df["AF_PAE"] = pae_values
            df["Atom 1 pLDDT"] = atom1_plddts
            df["Atom 2 pLDDT"] = atom2_plddts
            new_csv_path = os.path.join(
                cocomap_output_dir,
                "af_metrics",
                csv_f.replace(".csv", "_with_AF_metrics.csv")
            )
            os.makedirs(os.path.dirname(new_csv_path), exist_ok=True)
            df.to_csv(new_csv_path, index=False, encoding='utf-8')

            df_to_concat = pd.DataFrame(
                {
                    "Res. Name 1": df_res_names_1,
                    "Res. Number 1": df_res_1,
                    chains_set_1_prot: df_chains_1,
                    "Res. Name 2": df_res_names_2,
                    "Res. Number 2": df_res_2,
                    chains_set_2_prot: df_chains_2,
                    "Type of Interactions": df_type_of_interactions_w_clash,
                    "AF_PAE": pae_values,
                    "Atom 1 pLDDT": atom1_plddts,
                    "Atom 2 pLDDT": atom2_plddts,
                    "Atom 1": atom1_names,
                    "Atom 2": atom2_names,
                    "Metric Atom 1": metric_atom1_names,
                    "Metric Atom 2": metric_atom2_names,
                }
            )
            master_rows.append(df_to_concat)
        af_final_csv = pd.concat(master_rows, ignore_index=True)
        af_final_csv_path = os.path.join(
            cocomap_output_dir,
            "af_metrics",
            final_csvs[0].replace(".csv", "_with_AF_metrics.csv")
        )
        os.makedirs(os.path.dirname(af_final_csv_path), exist_ok=True)
        af_final_csv.to_excel(af_final_csv_path.replace(".csv", ".xlsx"), index=False)