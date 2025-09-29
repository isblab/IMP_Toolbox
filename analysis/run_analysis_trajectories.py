import os
import random
import time
import argparse
import getpass
from utils import generate_cmap, sanity_check_cores
try:
    import analysis_trajectories
    from analysis_trajectories import AnalysisTrajectories
except ImportError as e:
    if "No module named 'analysis_trajectories'" in str(e):
        print(
            """The 'analysis_trajectories' module is part of PMI_analysis.
            Please ensure that PMI_analysis is installed & added to PYTHONPATH.
            You can do this by running:
            `export PYTHONPATH=$PYTHONPATH:/path/to/PMI_analysis/pyext/src`
            or by adding the above line to your ~/.bash_profile."""
        )

_user = getpass.getuser()
start_t = time.time()

random.seed(47)
analysis_trajectories.color_palette = generate_cmap(
    n=40,
    scheme="contrasting"
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="PMI analysis on modeling trajectories."
    )
    parser.add_argument(
        "-m",
        "--modeling_output_path",
        type=str,
        default=f"/data/{_user}/imp_toolbox_test/modeling",
        help="Path to the modeling output directory. "
    )
    parser.add_argument(
        "-a",
        "--analysis_output_path",
        type=str,
        default=f"/data/{_user}/imp_toolbox_test/analysis",
        help="Path to the analysis output directory. "
    )
    parser.add_argument(
        "-s",
        "--run_start",
        type=int,
        default=1,
        help="Starting run number (default: 1)"
    )
    parser.add_argument(
        "-e",
        "--run_end",
        type=int,
        default=30,
        help="Ending run number (default: 30)"
    )
    parser.add_argument(
        "-i",
        "--run_interval",
        type=int,
        default=1,
        help="Interval between runs (default: 1)"
    )
    parser.add_argument(
        "-p",
        "--traj_dir_prefix",
        type=str,
        default="run_",
        help="Prefix for trajectory directories (default: 'run_')"
    )
    parser.add_argument(
        "-c",
        "--nproc",
        type=int,
        default=1,
        help="Number of cores to use for analysis (default: 1)"
    )
    parser.add_argument(
        "-b",
        "--burn_in_fraction",
        type=float,
        default=0.1,
        help="Fraction of data to discard as burn-in (default: 0.1)"
    )
    parser.add_argument(
        "-n",
        "--nskip",
        type=int,
        default=1,
        help="No. of consecutive frames to skip in the `AnalysisTrajectories`"
    )
    parser.add_argument(
        "-r",
        "--restraint_handles",
        type=str,
        nargs="+",
        default=[
            "GaussianEMRestraint:EM",
            "SingleAxisMinGaussianRestraint:SAMGR",
        ],
        help="List of restraint handles to analyze (default: \
            ['GaussianEMRestraint:EM', 'SingleAxisMinGaussianRestraint:SAMGR'])"
    )
    parser.add_argument(
        "-hr",
        "--hdbscan_restraint_handles",
        type=str,
        nargs="+",
        default=["EV_sum", "CR_sum", "EM_sum"],
        help="List of restraint handles to use for HDBSCAN clustering \
            (default: ['EV_sum', 'CR_sum', 'EM_sum'])"
    )
    parser.add_argument(
        "--min_cluster_size",
        type=int,
        default=150,
        help="Minimum cluster size for HDBSCAN (default: 150)"
    )
    parser.add_argument(
        "--min_samples",
        type=int,
        default=5,
        help="Minimum samples for HDBSCAN (default: 5)"
    )
    args = parser.parse_args()

sanity_check_cores(args.nproc)

os.makedirs(args.analysis_output_path, exist_ok=True)

traj_dirs = [
    f"{args.modeling_output_path}/{args.traj_dir_prefix}{i}"
    for i in range(args.run_start, args.run_end + 1, args.run_interval)
]

assert all(os.path.exists(traj_dir) for traj_dir in traj_dirs), \
    "One or more trajectory directories do not exist."

res_handles1 = [h.split(":")[1] for h in args.restraint_handles]
res_handles2 = [
    h.replace("_sum", "") for h in args.hdbscan_restraint_handles
    if h not in ["EV_sum", "CR_sum"]
]

assert set(res_handles2).issubset(set(res_handles1)), \
    f"""
    All HDBSCAN restraint handles must be included in the restraint handles.
    Restraint handles: {res_handles1}
    HDBSCAN restraint handles: {res_handles2}
    """

at_obj = AnalysisTrajectories(
    out_dirs=traj_dirs,
    dir_name=args.traj_dir_prefix,
    nproc=args.nproc,
    analysis_dir=args.analysis_output_path,
    burn_in_fraction=args.burn_in_fraction,
    nskip=args.nskip,
    detect_equilibration=True,
)

# Usual restraints
at_obj.set_analyze_Connectivity_restraint()
at_obj.set_analyze_Excluded_volume_restraint()
# at_obj.set_analyze_EM_restraint()

for handle in args.restraint_handles:
    restraint, short_name = handle.split(':')
    at_obj.set_analyze_score_only_restraint(
        handle=restraint,
        short_name=short_name,
        do_sum=True,
    )

at_obj.read_stat_files()
at_obj.write_models_info()

at_obj.hdbscan_clustering(
    args.hdbscan_restraint_handles,
    min_cluster_size=args.min_cluster_size,
    min_samples=args.min_samples,
    skip=args.nskip,
)

end_time = time.time()
print(f"Analysis completed in {end_time - start_t:.2f} seconds.")