import os
import pandas as pd
import subprocess
from IMP_Toolbox.chimerax.density_map import parse_chimerax_correlation_log
from IMP_Toolbox.constants.imp_toolbox_constants import (
    CHIMERAX_RUN_CMD,
    ChimeraXCommand,
    CHIMERAX_COMMANDS,
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
) -> str:
    """ For a given reference map and candidate model maps corresponding to
    different GMMs, get the candidate model that best fits the reference map
    (correlation >= 0.95) while having the least number of Gaussians.

    ## Arguments:

    - **data_file (str)**:<br />
        Path to the reference map file.

    - **centers_list (list)**:<br />
        List of centers for the GMMs.

    - **gmm_dir (str)**:<br />
        Directory containing GMM files.

    - **log_dir (str)**:<br />
        Directory to save the ChimeraX log file.

    - **data_sd_level (float)**:<br />
        Standard deviation level for the data volume.

    - **sd_level (float)**:<br />
        Standard deviation level for the GMM volumes.

    - **metric (str, optional):**:<br />
        Metric to use for fitting. Defaults to "cam".

    - **chimerax_run_cmd (str, optional):**:<br />
        Command to run ChimeraX. Defaults to the command for the flatpak

    - **threshold (float, optional):**:<br />
        Correlation threshold for selecting the best GMM. Defaults to 0.95.

    ## Returns:

    - **str**:<br />
        Base name of the best GMM file without the extension.
    """

    commands_run = [
        CHIMERAX_COMMANDS[ChimeraXCommand.OPEN].substitute(file_path=data_file),
        CHIMERAX_COMMANDS[ChimeraXCommand.VOLUME_THRESHOLD].substitute(
            model_idxs="#1",
            level_type=ChimeraXVolumeLevelType.SD_LEVEL,
            level_val=data_sd_level,
        ),
    ]

    model_count = 2
    f_name = os.path.basename(data_file).split(f".{FileFormat.MRC}")[0]
    valid_gmm_files = [
        f"{GMMParams.gmm_mrc_name.substitute(f_name=f_name, n_gaussian=centers_list[i])}.{FileFormat.MRC}"
        for i in range(len(centers_list))
    ]

    for gmm_file in sorted(os.listdir(gmm_dir)):
        if gmm_file in valid_gmm_files:
            commands_run.extend([
                CHIMERAX_COMMANDS[ChimeraXCommand.OPEN].substitute(
                    file_path=os.path.join(gmm_dir, gmm_file)
                ),
                CHIMERAX_COMMANDS[ChimeraXCommand.VOLUME_THRESHOLD].substitute(
                    model_idxs=f"#{model_count}",
                    level_type=ChimeraXVolumeLevelType.SD_LEVEL,
                    level_val=sd_level,
                ),
                CHIMERAX_COMMANDS[ChimeraXCommand.FITMAP].substitute(
                    model_idxs=f"#{model_count}",
                    ref_idxs="#1",
                    metric=metric,
                    shift_val=MiscStrEnum.FALSE,
                    rotate_val=MiscStrEnum.FALSE,
                    envelop_val=MiscStrEnum.FALSE,
                    fitmap_max_steps="0",
                    zeros_val=MiscStrEnum.FALSE,
                ),
                # f"measure correlation #{model_count} in_map #1 envelope true",
            ])
            model_count += 1

    log_path = os.path.join(
        log_dir,
        GMMParams.log_file_name.substitute(f_name=f_name) + f".{FileFormat.TXT}"
    )

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
) -> tuple[str, list]:
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

    correlation_list = parse_chimerax_correlation_log(
        chimerax_log=chimerax_log,
        command=ChimeraXCommand.FITMAP
    )

    correlation_list = [
        {
            MiscStrEnum.MODEL_MRC: item[MiscStrEnum.MODEL_MRC],
            MiscStrEnum.REF_MRC: item[MiscStrEnum.REF_MRC],
            CorrelationMetric.CORRELATION: item[CorrelationMetric.CORRELATION],
            CorrelationMetric.CAM: item[CorrelationMetric.CAM],
            CorrelationMetric.OVERLAP: item[CorrelationMetric.OVERLAP],
            MiscStrEnum.NUM_PTS: item[MiscStrEnum.NUM_PTS],
            MiscStrEnum.NUM_GAUSSIANS: int(item[MiscStrEnum.MODEL_MRC].split("_")[-1]),
        } for item in correlation_list
    ]

    correlation_list = sorted(correlation_list, key=lambda x: x[MiscStrEnum.NUM_GAUSSIANS])
    best_gmm = None
    for item in correlation_list:
        if float(item[metric]) >= threshold:
            best_gmm = item[MiscStrEnum.MODEL_MRC].split(f".{FileFormat.MRC}")[0]
            print(f"Selected GMM: {best_gmm} with {metric}: {item[metric]}")
            break

    df = pd.DataFrame(correlation_list)
    out_html = chimerax_log.replace(f".{FileFormat.TXT}", f".{FileFormat.HTML}")
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

def highlight_greater_than(val: float, threshold: float=GMMParams.threshold) -> str:
    """ Highlight values if greater than threshold.

    ## Arguments:

    - **val (float)**:<br />
        Value to compare against the threshold.

    - **threshold (float, optional):**:<br />
        Threshold to compare the value against. Defaults to 0.95.

    ## Returns:

    - **str**:<br />
        CSS style string for highlighting the value if it is greater than the threshold.
    """

    color = 'yellow' if float(val) > threshold else ''
    return f'background-color: {color}'