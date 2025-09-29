import os
import time
import argparse
import getpass
try:
    from analysis_trajectories import AnalysisTrajectories
except ImportError as e:
    if "No module named `analysis_trajectories`" in str(e):
        print(
            """The `analysis_trajectories` module is part of PMI_analysis.
            Please ensure that PMI_analysis is installed & added to PYTHONPATH.
            You can do this by running:
            export PYTHONPATH=$PYTHONPATH:/path/to/PMI_analysis/pyext/src
            or by adding the above line to your ~/.bash_profile."""
        )

_user = getpass.getuser()
start_t = time.time()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract good scoring models"
    )
    parser.add_argument(
        "--modeling_output_path",
        type=str,
        default=f"/data/{_user}/imp_toolbox_test/modeling",
        help="Path to the modeling output directory. "
    )
    parser.add_argument(
        "--analysis_output_path",
        type=str,
        default=f"/data/{_user}/imp_toolbox_test/analysis",
        help="Path to the analysis output directory. "
    )
    parser.add_argument(
        "--run_start",
        type=int,
        default=1,
        help="Starting run number (default: 1)"
    )
    parser.add_argument(
        "--cluster_id",
        type=int,
        default=1,
        help="Cluster ID to extract models from (default: 1)"
    )
    parser.add_argument(
        "--filter_applied",
        type=str,
        default='False',
        help="Whether variable_filter was run before this script"
    )
    parser.add_argument(
        "--variable_filter_output_dir",
        type=str,
        default=f"/data/{_user}/imp_toolbox_test/analysis/variable_filter_output",
        help="Directory where variable filter outputs are stored"
    )
    parser.add_argument(
        "--run_end",
        type=int,
        default=30,
        help="Ending run number (default: 30)"
    )
    parser.add_argument(
        "--run_interval",
        type=int,
        default=1,
        help="Interval between runs (default: 1)"
    )
    parser.add_argument(
        "--traj_dir_prefix",
        type=str,
        default="run_",
        help="Prefix for trajectory directories (default: 'run_')"
    )
    parser.add_argument(
        "--nproc",
        type=int,
        default=1,
        help="Number of cores to use for analysis (default: 1)"
    )
    parser.add_argument(
        "--burn_in_fraction",
        type=float,
        default=0.1,
        help="Fraction of data to discard as burn-in (default: 0.1)"
    )
    parser.add_argument(
        "--nskip",
        type=int,
        default=1,
        help="Number of consecutive frames to skip (default: 1)"
    )
    args = parser.parse_args()

    traj_dirs = [
        f"{args.modeling_output_path}/{args.traj_dir_prefix}{i}"
        for i in range(args.run_start, args.run_end + 1, args.run_interval)
    ]

    assert all(os.path.exists(traj_dir) for traj_dir in traj_dirs), \
        "One or more trajectory directories do not exist."

    # Load module
    at_obj = AnalysisTrajectories(
        out_dirs=traj_dirs,
        dir_name=args.traj_dir_prefix,
        nproc=args.nproc,
        analysis_dir=args.analysis_output_path,
        burn_in_fraction=args.burn_in_fraction,
        nskip=args.nskip,
    )

    if args.filter_applied == "True":

        print(
            "Variable filter was applied. \
            Extracting models from good_scoring_models files."
        )
        df_A = at_obj.get_models_to_extract(
            os.path.join(
                args.variable_filter_output_dir,
                f"good_scoring_models_A_cluster{str(args.cluster_id)}_detailed.csv"
            )
        )
        df_B = at_obj.get_models_to_extract(
            os.path.join(
                args.variable_filter_output_dir,
                f"good_scoring_models_B_cluster{str(args.cluster_id)}_detailed.csv"
            )
        )

    else:

        print("Variable filter was NOT applied. Extracting models directly.")
        df_A = at_obj.get_models_to_extract(
            os.path.join(
                args.analysis_output_path,
                f"selected_models_A_cluster{str(args.cluster_id)}_detailed.csv"
            )
        )
        df_B = at_obj.get_models_to_extract(
            os.path.join(
                args.analysis_output_path,
                f"selected_models_B_cluster{str(args.cluster_id)}_detailed.csv"
            )
        )

    rmf_file_out_A = f"A_models_clust{str(args.cluster_id)}.rmf3"
    rmf_file_out_B = f"B_models_clust{str(args.cluster_id)}.rmf3"

    # Extract models from the RMF files into a single RMF file for each sample
    at_obj.do_extract_models_single_rmf(
        df_A,
        rmf_file_out_A,
        args.modeling_output_path,
        args.analysis_output_path,
        scores_prefix="A_models_clust" + str(args.cluster_id),
    )

    at_obj.do_extract_models_single_rmf(
        df_B,
        rmf_file_out_B,
        args.modeling_output_path,
        args.analysis_output_path,
        scores_prefix="B_models_clust" + str(args.cluster_id)
    )

    end_t = time.time()
    print(f"Extracted models in {end_t - start_t:.2f} seconds.")