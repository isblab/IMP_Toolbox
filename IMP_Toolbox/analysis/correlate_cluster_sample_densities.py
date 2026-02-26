import os
import argparse
import subprocess
import pandas as pd
from IMP_Toolbox.structure.density_map import (
    extract_voxel_data,
    get_correlation_metrics,
)
from IMP_Toolbox.chimerax.density_map import parse_chimerax_correlation_log
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

if __name__ == "__main__":

    parser = argparse.ArgumentParser("""
    Correlate LPD maps from Sample_A and Sample_B directories.
    1. If --use_combined_map is set, combines all LPD maps in Sample_A and Sample_B
       respectively and computes correlation between the combined maps.
    2. If --use_combined_map is not set, computes correlation for each pair of LPD maps
       from Sample_A and Sample_B.
    3. Supports two modes for correlation calculation: 'pasani' and 'chimerax'.
       - 'pasani': Uses custom interpolation and correlation calculation.
       - 'chimerax': Uses ChimeraX commands to compute correlations.
    """)

    parser.add_argument(
        '--sampcon_cluster_path',
        type=str,
        required=True,
        help='Path to the directory containing Sample_A and Sample_B'
    )

    parser.add_argument(
        "--mode",
        type=str,
        choices=list(MRCComparisonMode),
        default=MRCComparisonMode.PASANI,
        help="Method to use for correlation calculation. Default is 'pasani'."
    )

    parser.add_argument(
        "--use_combined_map",
        action='store_true',
        help="Whether to use the combined map for correlation calculation. If not set, individual LPD maps will be compared."
    )

    parser.add_argument(
        "--voxel_size",
        type=float,
        default=5.0,
        help="Voxel size to use for correlation calculation. Default is 5.0."
    )

    parser.add_argument(
        "--zeros",
        action='store_true',
        help="Whether to include zero valued grid points in correlation calculation."
    )

    args = parser.parse_args()

    voxel_size = args.voxel_size
    zeros = bool(args.zeros)

    pa = os.path.join(args.sampcon_cluster_path, MiscStrEnum.SAMPLE_A)
    pb = os.path.join(args.sampcon_cluster_path, MiscStrEnum.SAMPLE_B)

    fa = sorted([os.path.join(pa,f"{f}") for f in os.listdir(pa) if CLUSTER_MRC_PREFIX in f])
    fb = sorted([os.path.join(pb,f"{f}") for f in os.listdir(pb) if CLUSTER_MRC_PREFIX in f])

    if args.mode == MRCComparisonMode.PASANI:
        # this is equivalent to the following chimerax command:
        # fitmap #model_mrc inmap #ref_mrc shift false rotate false envelope false maxSteps 0 zeros <true/false>
        # the difference is in the interpolation method used for resampling the maps
        # hence, the results may differ slightly from the ChimeraX results

        if args.use_combined_map:

            voxel_data_a = extract_voxel_data(mrc_files=fa)
            voxel_data_b = extract_voxel_data(mrc_files=fb)

            overlap, corr, corr_over_mean, pts = get_correlation_metrics(
                voxel_data1=voxel_data_a,
                voxel_data2=voxel_data_b,
                voxel_size=voxel_size,
                zeros=zeros,
            )

            correlation_list = [{
                MiscStrEnum.MODEL_MRC: MiscStrEnum.SAMPLE_A,
                MiscStrEnum.REF_MRC: MiscStrEnum.SAMPLE_B,
                CorrelationMetric.OVERLAP: overlap,
                CorrelationMetric.CORRELATION: corr,
                CorrelationMetric.CAM: corr_over_mean,
                MiscStrEnum.NUM_PTS: pts,
            }]

        else:

            correlation_list = []
            for f_a, f_b in zip(fa, fb):

                voxel_data_a = extract_voxel_data(mrc_files=[f_a])
                voxel_data_b = extract_voxel_data(mrc_files=[f_b])

                overlap, corr, corr_over_mean, pts = get_correlation_metrics(
                    voxel_data1=voxel_data_a,
                    voxel_data2=voxel_data_b,
                    voxel_size=voxel_size,
                    zeros=zeros,
                )

                correlation_list.append({
                    MiscStrEnum.MODEL_MRC: os.path.basename(f_b).replace(f".{FileFormat.MRC}", ""),
                    MiscStrEnum.REF_MRC: os.path.basename(f_a).replace(f".{FileFormat.MRC}", ""),
                    CorrelationMetric.OVERLAP: overlap,
                    CorrelationMetric.CORRELATION: corr,
                    CorrelationMetric.CAM: corr_over_mean,
                    MiscStrEnum.NUM_PTS: pts,
                })

    elif args.mode == MRCComparisonMode.CHIMERAX:

        zeros_ = MiscStrEnum.TRUE if zeros else MiscStrEnum.FALSE

        log_path = os.path.join(args.sampcon_cluster_path, 'correlation_logs.txt')

        if args.use_combined_map:

            commands_run = []

            for f1 in fa:
                commands_run += [CHIMERAX_COMMANDS[ChimeraXCommand.OPEN].substitute(file_path=f1)]

            for f2 in fb:
                commands_run += [CHIMERAX_COMMANDS[ChimeraXCommand.OPEN].substitute(file_path=f2)]

            commands_run += [
                CHIMERAX_COMMANDS[ChimeraXCommand.VOLUME_ADD].substitute(model_range=f"#1-{len(fa)}"),
                CHIMERAX_COMMANDS[ChimeraXCommand.VOLUME_ADD].substitute(model_range=f"#{len(fa)+1}-{len(fa)+len(fb)}"),
                CHIMERAX_COMMANDS[ChimeraXCommand.RENAME].substitute(
                    model_idx=f"#{len(fa)+len(fb)+1}",
                    new_name=MiscStrEnum.SAMPLE_A,
                ),
                CHIMERAX_COMMANDS[ChimeraXCommand.RENAME].substitute(
                    model_idx=f"#{len(fa)+len(fb)+2}",
                    new_name=MiscStrEnum.SAMPLE_B,
                ),
                CHIMERAX_COMMANDS[ChimeraXCommand.FITMAP].substitute(
                    model_idxs=f"#{len(fa)+len(fb)+1}",
                    ref_idxs=f"#{len(fa)+len(fb)+2}",
                    metric=CorrelationMetric.CORRELATION, # doesn't matter which metric is used as shift and rotate are false
                    shift_val=MiscStrEnum.FALSE,
                    rotate_val=MiscStrEnum.FALSE,
                    envelop_val=MiscStrEnum.FALSE,
                    fitmap_max_steps=0,
                    zeros_val=zeros_,
                ),
                CHIMERAX_COMMANDS[ChimeraXCommand.FITMAP].substitute(
                    model_idxs=f"#{len(fa)+len(fb)+2}",
                    ref_idxs=f"#{len(fa)+len(fb)+1}",
                    metric=CorrelationMetric.CORRELATION, # doesn't matter which metric is used as shift and rotate are false
                    shift_val=MiscStrEnum.FALSE,
                    rotate_val=MiscStrEnum.FALSE,
                    envelop_val=MiscStrEnum.FALSE,
                    fitmap_max_steps=0,
                    zeros_val=zeros_,
                ),
                CHIMERAX_COMMANDS[ChimeraXCommand.CLOSE].substitute(model_idxs="all"),
            ]

            subprocess.run(
                f"{CHIMERAX_RUN_CMD} --exit --nogui --cmd " +
                f"'{"; ".join(commands_run)}' > {log_path}",
                shell=True,
                check=True,
            )

        else:

            for idx, (f1, f2) in enumerate(zip(fa, fb)):

                commands_run = [
                    CHIMERAX_COMMANDS[ChimeraXCommand.OPEN].substitute(file_path=f1),
                    CHIMERAX_COMMANDS[ChimeraXCommand.OPEN].substitute(file_path=f2),
                    CHIMERAX_COMMANDS[ChimeraXCommand.FITMAP].substitute(
                        model_idxs=f"#1",
                        ref_idxs=f"#2",
                        metric=CorrelationMetric.CORRELATION, # doesn't matter which metric is used as shift and rotate are false
                        shift_val=MiscStrEnum.FALSE,
                        rotate_val=MiscStrEnum.FALSE,
                        envelop_val=MiscStrEnum.FALSE,
                        fitmap_max_steps=0,
                        zeros_val=zeros_,
                    ),
                    CHIMERAX_COMMANDS[ChimeraXCommand.FITMAP].substitute(
                        model_idxs=f"#2",
                        ref_idxs=f"#1",
                        metric=CorrelationMetric.CORRELATION, # doesn't matter which metric is used as shift and rotate are false
                        shift_val=MiscStrEnum.FALSE,
                        rotate_val=MiscStrEnum.FALSE,
                        envelop_val=MiscStrEnum.FALSE,
                        fitmap_max_steps=0,
                        zeros_val=zeros_,
                    ),
                    CHIMERAX_COMMANDS[ChimeraXCommand.CLOSE].substitute(model_idxs="all"),
                ]

                log_mode = ">>" if idx > 0 else ">"

                subprocess.run(
                    f"{CHIMERAX_RUN_CMD} --exit --nogui --cmd " +
                    f"'{"; ".join(commands_run)}' {log_mode} {log_path}",
                    shell=True,
                    check=True,
                )

        correlation_list = parse_chimerax_correlation_log(
            chimerax_log=log_path,
            command=ChimeraXCommand.FITMAP,
        )

    df = pd.DataFrame(correlation_list)

    corr_ = df[CorrelationMetric.CORRELATION].astype(float)
    cam_ = df[CorrelationMetric.CAM].astype(float)
    overlap_ = df[CorrelationMetric.OVERLAP].astype(float)
    pts_ = df[MiscStrEnum.NUM_PTS].astype(int)

    print(f"Overlap: {overlap_.mean()} ± {overlap_.std()}")
    print(f"Average Correlation: {corr_.mean()} ± {corr_.std()}")
    print(f"Average CAM: {cam_.mean()} ± {cam_.std()}")
    print(f"Number of points: {pts_.sum()}")

    out_csv = os.path.join(args.sampcon_cluster_path, 'correlation_results.csv')
    df.to_csv(out_csv, index=False)
    print(f"Correlation results saved to {out_csv}")