import os
import warnings
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