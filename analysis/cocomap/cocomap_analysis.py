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

    proteins = result_metadata.get("proteins", [])
    if len(proteins) == 2:
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
            df = pd.read_csv(csv_path, index_col=0)
            df = df.rename(columns=new_cols)
            df.to_csv(csv_path)

    final_csvs = [
        f for f in os.listdir(cocomap_output_dir)
        if f.endswith("final_file.csv")
    ]

    for final_csv in final_csvs:
        final_csv_path = os.path.join(cocomap_output_dir, final_csv)
        df = pd.read_csv(final_csv_path, index_col=0)
        df["Type of Interactions"] = df["Type of Interactions"].str.split("; ")
        df_exploded = df.explode("Type of Interactions")
        df_exploded = df_exploded[
            df_exploded["Type of Interactions"].str.strip() != ""
        ]
        df_exploded = df_exploded.drop_duplicates()
        df_exploded.to_csv(final_csv_path, index=False)