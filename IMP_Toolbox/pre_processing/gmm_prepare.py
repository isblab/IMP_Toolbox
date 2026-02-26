import os
import pandas as pd
import subprocess
from IMP_Toolbox.constants.imp_toolbox_constants import (
    CHIMERAX_RUN_CMD,
    ChimeraXCommands,
    ChimeraXDefaults,
    ChimeraXVolumeLevelType,
    GMMParams,
    CorrelationMetric,
    MiscStrEnum,
    FileFormat,
)

def replace_data_fn(out_file):
    """ Replace the full path in the `data_fn` in the output txt file with the
    basename.

    ## Arguments:

    - **out_file (str)**:
        Path to the output file where the GMM parameters are saved.
    """

    with open(out_file, 'r') as f:
        lines = f.readlines()

    if "# data_fn" in lines[1]:
        _, attr_val = lines[1].split(":")
        attr_val = attr_val.strip()
        new_attr_val = os.path.basename(attr_val)
        lines[1] = f"# data_fn: {new_attr_val}\n"

    with open(out_file, 'w') as f:
        f.writelines(lines)

def get_best_gmm(
    data_file: str,
    centers_list: list,
    gmm_dir: str,
    log_dir: str,
    data_sd_level: float,
    sd_level: float,
    metric: str = GMMParams.correlation_metric,
    chimerax_run_cmd: str = CHIMERAX_RUN_CMD,
    threshold: float = GMMParams.threshold,
):
    """ For a given reference map and candidate model maps corresponding to
    different GMMs, get the candidate model that best fits the reference map
    (correlation >= 0.95) while having the least number of Gaussians.

    Args:
        data_file (str): Description of the layer pattern to match GMM files.
        centers_list (list): List of centers for the GMMs.
        gmm_dir (str): Directory containing GMM files.
        log_dir (str): Directory to save the ChimeraX log file.
        data_sd_level (float): Standard deviation level for the data volume.
        sd_level (float): Standard deviation level for the GMM volumes.
        metric (str): Metric to use for fitting. Default is "cam".
        chimerax_run_cmd (str): Command to run ChimeraX. Default is for the
            flatpak installation. Change this to the appropriate command to run
            ChimeraX on your system.

    Returns:
        commands_run (list): Updated list of commands run in ChimeraX.
    """

    commands_run = [
        ChimeraXCommands.open_cmd.substitute(file_path=data_file),
        ChimeraXCommands.volume_threshold_cmd.substitute(
            model_count=1,
            level_type=ChimeraXVolumeLevelType.SD_LEVEL,
            level_val=data_sd_level
        )
    ]

    model_count = 2
    f_name = os.path.basename(data_file).split(f".{FileFormat.MRC}")[0]
    valid_gmm_files = [
        GMMParams.gmm_mrc_name.substitute(
            f_name=f_name,
            n_gaussian=centers_list[i]
        ) for i in range(len(centers_list))
    ]

    for gmm_file in sorted(os.listdir(gmm_dir)):
        if gmm_file in valid_gmm_files:
            commands_run.extend([
                ChimeraXCommands.open_cmd.substitute(
                    file_path=os.path.join(gmm_dir, gmm_file)
                ),
                ChimeraXCommands.volume_threshold_cmd.substitute(
                    model_count=model_count,
                    level_type=ChimeraXVolumeLevelType.SD_LEVEL,
                    level_val=sd_level
                ),
                ChimeraXCommands.fitmap_cmd.substitute(
                    model_count=model_count,
                    metric=metric,
                    shift_val=MiscStrEnum.FALSE,
                    rotate_val=MiscStrEnum.FALSE,
                    envelop_val=MiscStrEnum.FALSE,
                    fitmap_max_steps=0,
                    zeros_val=MiscStrEnum.FALSE,
                ),
                # f"measure correlation #{model_count} in_map #1 envelope true",
            ])
            model_count += 1

    log_path = os.path.join(
        log_dir,
        GMMParams.log_file_name.substitute(f_name=f_name) + f".{FileFormat.TXT}"
    )
    # os.system(
    #     f"{chimerax_run_cmd} --exit --nogui --cmd " +
    #     f"'{"; ".join(commands_run)}' > {log_path}"
    # )

    subprocess.run(
        f"{chimerax_run_cmd} --exit --nogui --cmd " +
        f"'{"; ".join(commands_run)}' > {log_path}",
        shell=True,
        check=True,
    )

    best_gmm, _ = compare_gmms(
        chimerax_log=log_path,
        metric=metric,
        threshold=threshold
    )

    return best_gmm

def compare_gmms(
    chimerax_log: str,
    metric: str = GMMParams.correlation_metric,
    threshold: float = GMMParams.threshold,
):
    """ Get the best GMM based on the correlation metrics from the
    ChimeraX log file.

    ## Arguments:

    - **chimerax_log (str)**:<br />
        Path to the ChimeraX log file containing correlation results from the
        command `measure correlation`.

    - **metric (str, optional):**:<br />
        Metric to use for selecting the best GMM. Defaults to "cam".

    ## Returns:

    - **tuple**:<br />
        - Base name of the best GMM file without the extension.
        - Correlation of each GMM as a list of dictionaries.
    """

    assert metric in list(CorrelationMetric), (
        f"Metric must be from {list(CorrelationMetric)}. Got {metric}."
    )

    log_lines = []
    with open(chimerax_log, 'r') as f:
        log_lines = f.readlines()

    correlation_list = []

    for idx, line in enumerate(log_lines):
        if (
            idx + 1 >= len(log_lines) or
            line.startswith("Correlation of") is False
        ):
            continue
        line_parts = line.split(" ")
        nxt_line_parts = log_lines[idx + 1].split(" ")
        model_mrc = line_parts[2]
        ref_mrc = line_parts[-2]
        correlation = nxt_line_parts[2].split(",")[0].strip()
        cam = nxt_line_parts[-1].strip()
        correlation_list.append({
            "model_mrc": model_mrc,
            "ref_mrc": ref_mrc,
            "correlation": correlation,
            "cam": cam,
            "gaussians": int(model_mrc.split("_")[-1].split(f".{FileFormat.MRC}")[0]),
        })

    correlation_list = sorted(correlation_list, key=lambda x: x["gaussians"])
    best_gmm = None
    for item in correlation_list:
        if float(item[metric]) >= threshold:
            best_gmm = item["model_mrc"].split(f".{FileFormat.MRC}")[0]
            print(f"Selected GMM: {best_gmm} with {metric}: {item[metric]}")
            break

    df = pd.DataFrame(correlation_list)
    out_html = chimerax_log.replace(f".{FileFormat.TXT}", f"_{FileFormat.HTML}")
    styled_df = df.style.map(
        lambda val: highlight_greater_than(val, threshold=threshold),
        subset=[CorrelationMetric.CORRELATION, CorrelationMetric.CAM]
    )
    styled_df.to_html(
        out_html,
        index=False,
        table_attributes='style="width:100%"',
        classes='table table-striped table-bordered',
    )

    return best_gmm, correlation_list

def highlight_greater_than(val, threshold=GMMParams.threshold):
    color = 'yellow' if float(val) > threshold else ''
    return f'background-color: {color}'