import os
import yaml
import argparse
import subprocess
import pandas as pd
from IMP_Toolbox.density.compare import (
    extract_voxel_data,
    get_correlation_metrics,
)
from IMP_Toolbox.utils.file_helpers import read_json
from IMP_Toolbox.constants.imp_toolbox_constants import (
    CLUSTER_MRC_PREFIX,
    CHIMERAX_RUN_CMD,
    CHIMERAX_COMMANDS,
    FileFormat,
    MRCComparisonMode,
    CorrelationMetric,
    ChimeraXCommand,
    MiscStrEnum,
)
import getpass
_user = getpass.getuser()

def get_voxel_correlation(
    key: str,
    model_maps: list[str],
    exp_maps: list[str],
    voxel_size: float | None = None,
) -> pd.DataFrame:
    """ Get correlation metrics between the LPDs and the experimental density maps.

    ## Arguments:

    - **key (str)**:<br />
        Name of the LPD.

    - **model_maps (list[str])**:<br />
        List of paths to the model MRC files (LPDs).

    - **exp_maps (list[str])**:<br />
        List of paths to the experimental MRC files.

    ## Returns:

    - **pd.DataFrame**:<br />
        DataFrame containing the correlation metrics between the model and experimental maps.
    """
    voxel_data_a = extract_voxel_data(mrc_files=model_maps)
    voxel_data_b = extract_voxel_data(mrc_files=exp_maps)

    overlap, corr, corr_over_mean, pts = get_correlation_metrics(
        voxel_data1=voxel_data_a,
        voxel_data2=voxel_data_b,
        voxel_size=voxel_size,
        zeros=False,
    )

    correlation_list = [{
        MiscStrEnum.MODEL_MRC: key,
        MiscStrEnum.REF_MRC: key.replace("_lpd", "_exp"),
        CorrelationMetric.OVERLAP: overlap,
        CorrelationMetric.CORRELATION: corr,
        CorrelationMetric.CAM: corr_over_mean,
        MiscStrEnum.NUM_PTS: pts,
    }]

    df = pd.DataFrame(correlation_list)

    return df

def main(
    em_data_info: dict,
    output_dir: str,
) -> None:
    """ main function to compute fit to EM data.

    ## Arguments:

    - **em_data_info (dict)**:<br />
        Dictionary containing the EM data information, including model and
        experimental maps. For e.g.:

        {"pg_layer_lpd": {
            "model_maps": ["/path/to/LPD_DP-S.mrc", "/path/to/LPD_PG-S.mrc"],
            "exp_maps": ["/path/to/pg_modeled_layer.mrc"]
            },
        "pkp_layer_lpd": {
            "model_maps": ["/path/to/LPD_PKP-C.mrc"],
            "exp_maps": ["/path/to/pkp_modeled_layer.mrc"]
        }}

    - **output_dir (str)**:<br />
        Directory to save the correlation results.
    """

    if em_data_info.endswith(".json"):
        em_data_info: dict = read_json(em_data_info)
    elif em_data_info.endswith(".yaml") or em_data_info.endswith(".yml"):
        em_data_info: dict = yaml.load(open(em_data_info), Loader=yaml.FullLoader)
    else:
        raise ValueError("Input file must be JSON or YAML format.")

    meta_df = pd.DataFrame()

    for key, em_data in em_data_info.items():
        df = get_voxel_correlation(
            key=key,
            model_maps=em_data["model_maps"],
            exp_maps=em_data["exp_maps"],
            voxel_size=em_data.get("voxel_size", None)
        )

        meta_df = pd.concat([meta_df, df], ignore_index=True)

    out_csv = os.path.join(output_dir, 'correlation_LPD_experimentalMap.csv')
    meta_df.to_csv(out_csv, index=False)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the JSON/YAML file containing EM data information.",
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to save the correlation results.",
        default=f"/home/{_user}/Projects/cardiac_desmosome/analysis/prod_runs/output_trans_dsc_5716/set5/fit_to_em",
    )

    args = parser.parse_args()

    main(
        em_data_info=args.input,
        output_dir=args.output_dir,
    )