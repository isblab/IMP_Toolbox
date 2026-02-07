import os
import sys
import time
import random
import argparse
import logging
import getpass
_user = getpass.getuser()

logging.basicConfig(
    filename=os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'end_to_end_analysis.log'
    ),
    filemode="a+",
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO
)

#NOTE: Change these as per your setup and requirements
# do not to add CR and EV to RESTRAINT_HANDLES as they are analyzed by default
RESTRAINT_HANDLES = [
    "GaussianEMRestraint:EM",
    "SingleAxisMinGaussianRestraint:SAMGR",
]
# restaint sums to be used for HDBSCAN clustering
# you need to add CR_sum and EV_sum if you wish to use them for clustering
HDBSCAN_RESTRAINT_HANDLES = [
    "EV_sum",
    "EM_sum",
]
TRAJ_DIR_PREFIX = "run_"
PRISM_PATH = "/home/$USER/IMP_OMG/prism"
SAMPCON_PATH = "/home/$USER/IMP_OMG/imp-sampcon"
IMP_TOOLOBX_PATH = "/home/$USER/Omkar/IMP_TOOLBOX/IMP_Toolbox"
SYSNAME = "cardiac_desmosome"
SAMPCON_DENSITY_TXT = f"/data/{_user}/imp_toolbox_test/input/density_sampcon.txt"
MODEL_CAP = 30000 # for variable filter

random.seed(47)

def return_major_cluster(hdbscan_log_path: str):
    """
    Return the cluster with the maximum number of models from HDBSCAN
    clustering summary file.

    Args:
        hdbscan_log_path (str): Path to the HDBSCAN clustering summary file.

    Returns:
        tuple: (major_cluster_idx, number_of_models)
    """

    with open (hdbscan_log_path, "r") as f:
        cluster_summary = f.readlines()

    models_count = {}
    #Get the index for the N_models header
    index = cluster_summary[0].split(",").index("N_models")

    # Get cluster no. and model counts for all clusters.
    for i in range(1,len(cluster_summary)):
        line = cluster_summary[i].split(",")
        if line[0] != "-1":
            models_count[int(line[0])] = int(line[index])
        else:
            continue

    if len(models_count) == 0:
        raise ValueError(
            f"""No clusters found in HDBSCAN summary file {hdbscan_log_path}.
            Change the parameters for HDBSCAN clustering and try again."""
        )

    major_clust = [
        (i,models_count[i]) for i in range(len(models_count))
        if models_count[i] == max(models_count.values())
    ][0]

    return major_clust

def run_analysis_trajectories(
    script_path: str,
    modeling_dir: str,
    analysis_dir: str,
    traj_dir_prefix: str,
    run_start: int,
    run_end: int,
    run_interval: int,
    nproc: int,
    burn_in_fraction: float,
    nskip: int,
    restraint_handles: list,
    hdbscan_restraint_handles: list,
    min_cluster_size: int,
    min_samples: int,
    logger: logging.Logger | None = None,
    save_log: bool = True,
    log_dir: str = None,
):
    """ Run the script `run_analysis_trajectories.py`

    - This uses `AnalysisTrajectories` from PMI_analysis to analyze the
      trajectories generated from the modeling script.
    - It performs HDBSCAN clustering on the selected restraint handles
      provided in `hdbscan_restraint_handles`.
    - The odd/even trajectories are split into sample A and B respectively

    The output includes:
    - Summary of HDBSCAN clustering
    - Plots for each trajectory denoting the trajectory of restraint scores
    - Plot indicating the contributions of each trajectory in each cluster
        towards sample A and B
    - Pairwise plots for the restraint scores colored by clusters
    - Distribution of total score and score convergence for each cluster

    Args:
        script_path (str): Path to the `run_analysis_trajectories.py` script
        modeling_dir (str): Path to the modeling output directory
        analysis_dir (str): Path to the pmi_analysis output directory
        traj_dir_prefix (str): Prefix for trajectory directories (e.g., 'run_')
        run_start (int): Starting run number
        run_end (int): Ending run number
        run_interval (int): Interval between runs (e.g., 1 for every run)
        nproc (int): Number of cores to use for analysis
        burn_in_fraction (float): Fraction of data to discard as burn-in
        nskip (int): Number of consecutive frames to skip
        hdbscan_restraint_handles (list): List of restraint handles to use for
            HDBSCAN clustering
        logger (logging.Logger | None, optional): Logger for logging messages.
            If None, logging is skipped. Defaults to None.
        save_log (bool, optional): Whether to save the log file.
    """

    command = [
        "python", script_path,
        "--modeling_dir", modeling_dir,
        "--analysis_dir", analysis_dir,
        "--traj_dir_prefix", traj_dir_prefix,
        "--run_start", run_start,
        "--run_end", run_end,
        "--run_interval", run_interval,
        "--nproc", nproc,
        "--burn_in_fraction", burn_in_fraction,
        "--nskip", nskip,
        "--restraint_handles", *restraint_handles,
        "--hdbscan_restraint_handles", *hdbscan_restraint_handles,
        "--min_cluster_size", min_cluster_size,
        "--min_samples", min_samples
    ]

    if save_log:
        if log_dir is None:
            log_dir = os.path.join(os.getcwd(), 'logs')
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.join(log_dir, 'run_analysis_trajectories.log')
        command.append(f"> {log_path} 2>&1")

    if logger is not None:
        logger.info("Running analysis of trajectories with command:")
        logger.info(" ".join(map(str, command)))

    os.system(" ".join(map(str, command)))

def variable_filter(
    script_path: str,
    pmi_clust_idx: int,
    lowest_cutoff: float,
    highest_cutoff: float,
    step_size: float,
    model_cap: int,
    gsmsel_dir: str,
    output_dir: str,
    restraint_handles: list,
    logger: logging.Logger | None = None,
    save_log: bool = True,
    log_dir: str = None,
):
    """ Run the script `variable_filter.py`

    - This uses the `variable_filter.py` script from PMI_analysis to
      filter models from the major cluster obtained from
      `run_analysis_trajectories.py` in case greater than `model_cap`.

    Args:
        script_path (str): Path to the `variable_filter.py` script
        pmi_clust_idx (int): Index of the major cluster from HDBSCAN
        lowest_cutoff (float): stdev multiplier for the most stringent cutoff
        highest_cutoff (float): stdev multiplier for the most lenient cutoff
        step_size (float): step size for the cutoff
        model_cap (int): Maximum number of models to retain
        gsmsel_dir (str): directory where pmi_analysis csv files are stored
        output_dir (str): path to store the output of the script
        logger (logging.Logger): Logger for logging messages.
        save_log (bool, optional): Whether to save the log file. Defaults to True.
    """

    command = [
        "python", script_path,
        "--cluster_num", str(pmi_clust_idx),
        "--lowest_cutoff", str(lowest_cutoff),
        "--highest_cutoff", str(highest_cutoff),
        "--step_size", str(step_size),
        "--num_models", str(model_cap),
        "--gsmsel", gsmsel_dir,
        "--output_dir", output_dir,
        "--restraint_handles", *restraint_handles
    ]

    if save_log:
        if log_dir is None:
            log_dir = os.path.join(os.getcwd(), 'logs')
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.join(log_dir, 'variable_filter.log')
        command.append(f">> {log_path} 2>&1")

    if logger is not None:
        logger.info("Running variable filter with command:")
        logger.info(" ".join(map(str, command)))

    os.system(" ".join(map(str, command)))

def run_extract_models(
    script_path: str,
    modeling_dir: str,
    analysis_dir: str,
    traj_dir_prefix: str,
    run_start: int,
    run_end: int,
    run_interval: int,
    nproc: int,
    burn_in_fraction: float,
    nskip: int,
    cluster_id: int,
    filter_applied: bool,
    variable_filter_output_dir: str,
    logger: logging.Logger | None = None,
    save_log: bool = True,
    log_dir: str = None,
):
    """ Run the script `run_extract_models.py`

    - This extracts good scoring models mentioned in the csv files.
    - The csv files are the output of `variable_filter.py` if run before
      otherwise the output of `run_analysis_trajectories.py`
    - The extracted rmf3 and txt files indicating scores per frame are stored
      in `analysis_dir` for sample A and B respectively.

    Args:
        script_path (str): Path to the `run_extract_models.py` script
        modeling_dir (str): Path to the modeling output directory
        analysis_dir (str): Path to the analysis output directory
        traj_dir_prefix (str): Prefix for trajectory directories (e.g., 'run_')
        run_start (int): Starting run number
        run_end (int): Ending run number
        run_interval (int): Interval between runs (e.g., 1 for every run)
        nproc (int): Number of cores to use for analysis
        burn_in_fraction (float): Fraction of data to discard as burn-in
        nskip (int): Number of consecutive frames to skip
        cluster_id (int): Cluster ID to extract models from for the major cluster
            from the HDBSCAN clustering
        filter_applied (bool): Whether `variable_filter.py` was run or not
        logger (logging.Logger): Logger for logging messages.
        save_log (bool, optional): Whether to save the log file.
    """

    command = [
        "python", script_path,
        "--modeling_dir", modeling_dir,
        "--analysis_dir", analysis_dir,
        "--traj_dir_prefix", traj_dir_prefix,
        "--run_start", run_start,
        "--run_end", run_end,
        "--run_interval", run_interval,
        "--nproc", nproc,
        "--burn_in_fraction", burn_in_fraction,
        "--nskip", nskip,
        "--cluster_id", cluster_id
    ]

    if filter_applied and variable_filter_output_dir is not None:
        command.append("--filter_applied")
        command.append(
            f"--variable_filter_output_dir {variable_filter_output_dir}"
        )

    if save_log:
        if log_dir is None:
            log_dir = os.path.join(os.getcwd(), 'logs')
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.join(log_dir, 'run_extract_models.log')
        command.append(f"> {log_path} 2>&1")

    logger.info("Running extract models with command:")
    logger.info(" ".join(map(str, command)))
    os.system(" ".join(map(str, command)))

def exhaust(
    script_path: str,
    pmi_analysis_dir: str,
    sysname: str,
    scoreA: str,
    scoreB: str,
    rmfA: str,
    rmfB: str,
    density: str,
    gnuplot: bool,
    prism: bool,
    align: bool,
    mode: str,
    matrix_cores: int,
    cluster_cores: int,
    gridsize: int,
    logger: logging.Logger | None = None,
    save_log: bool = True,
    log_dir: str = None,
):
    """ Run the script `exhaust.py` from imp-sampcon

    - This checks for sampling exhaustiveness
    - Output files are stored in the directory from which the script is run
    - The output includes:
        - plots for cluster populations
        - txt files indicating model and cluster precision
        - score distribution of sample 1 and 2
        - score convergence plot
        - "cluster.<idx>" directories containing localization probability
          densities for each cluster
        - prism input `.npz` files for each cluster if `prism` is True

    Args:
        script_path (str): Path to the `exhaust.py` script
        pmi_analysis_dir (str): Path to the pmi_analysis output directory
        sysname (str): System name for the analysis (e.g., 'cardiac_desmosome')
        scoreA (str): txt file containing scores for sample A
        scoreB (str): txt file containing scores for sample B
        rmfA (str): rmf3 file containing models for sample A
        rmfB (str): rmf3 file containing models for sample B
        density (str): Path to the density txt file for imp-sampcon
        gnuplot (bool): whether to generate gnuplot files
        prism (bool): whether to generate prism input files
        align (bool): whether to align the models before analysis
        mode (str): mode for imp-sampcon (e.g., 'cpu_omp')
        matrix_cores (int): number of cores for rmsd calculation
        cluster_cores (int): number of cores for clustering
        gridsize (int): grid size for clustering
        logger (logging.Logger | None, optional): Logger for logging messages.
        save_log (bool, optional): Whether to save the log file.
    """

    command = [
        "python", script_path,
        "--sysname", sysname,
        "--mode", mode,
        "--matrix-cores", matrix_cores,
        "--cluster-cores", cluster_cores,
        "--density", density,
        "--gridsize", gridsize,
        "--scoreA", scoreA,
        "--scoreB", scoreB,
        "--rmfA", rmfA,
        "--rmfB", rmfB,
        "--path", pmi_analysis_dir
    ]

    if gnuplot:
        command.append("--gnuplot")
    if prism:
        command.append("--prism")
    if align:
        command.append("--align")

    if save_log:
        if log_dir is None:
            log_dir = os.path.join(os.getcwd(), 'logs')
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.join(log_dir, 'exhaust.log')
        command.append(f">> {log_path} 2>&1")

    if logger is not None:
        logger.info("Running exhaust with command:")
        logger.info(" ".join(map(str, command)))

    os.system(" ".join(map(str, command)))

def correlate_cluster_sample_LPDs(
    script_path: str,
    sampcon_cluster_path: str,
    mode: str,
    use_combined_map: bool,
    logger: logging.Logger | None = None,
    save_log: bool = True,
    log_dir: str = None,
):
    command = [
        "python", script_path,
        "--sampcon_cluster_path", sampcon_cluster_path,
        "--mode", mode,
    ]

    if use_combined_map:
        command.append("--use_combined_map")

    if save_log:
        if log_dir is None:
            log_dir = os.path.join(os.getcwd(), 'logs')
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.join(log_dir, 'correlate_cluster_sample_LPDs.log')
        command.append(f">> {log_path} 2>&1")

    if logger is not None:
        logger.info("Running correlation of cluster sample LPDs with command:")
        logger.info(" ".join(map(str, command)))

    os.system(" ".join(map(str, command)))

def fit_pdb_to_ccm(
    script_path: str,
    ccm_file: str,
    input_config: str,
    output_dir: str,
    logger: logging.Logger | None = None,
):
    command = [
        "python", script_path,
        "--ccm_file", ccm_file,
        "--input", input_config,
        "--output_dir", output_dir,
    ]

    if logger is not None:
        logging.info("Running fit PDB to CCM with command:")
        logging.info(" ".join(map(str, command)))

    os.system(" ".join(map(str, command)))

def prism_annotate(
    script_path: str,
    input: str,
    input_type: str,
    voxel_size: int,
    return_spread: bool,
    output: str,
    classes: int,
    cores: int,
    models: float,
    n_breaks: int,
    resolution: int,
    subunit: str | None,
    selection: None,
    logger: logging.Logger | None = None,
    save_log: bool = True,
    log_dir: str = None,
):
    """ Run the script `main.py` from prism

    - This generates annotations for indicating precision on the bead models

    Args:
        script_path (str): Path to the `main.py` script from prism
        input (str): Path to the input file (e.g., .npz file from imp-sampcon)
        input_type (str): Type of the input file (e.g., 'npz')
        voxel_size (int): Voxel size in Angstroms
        return_spread (bool): Whether to return bead spread information
        output (str): Path to the output file
        classes (int): Number of classes for annotation (2, or 3)
        cores (int): Number of cores to use
        models (float): Fraction of models to use (0-1)
        n_breaks (int): Number of breaks for jenkspy
        resolution (int): Resolution as number of residues per bead
        subunit (str | None): name of the subunit to analyze (None for all)
        selection (None): selection within the subunit (None for all)
        logger (logging.Logger): Logger for logging messages.
        save_log (bool, optional): Whether to save the log file.
    """

    command = [
        "python", script_path,
        "--input", input,
        "--input_type", input_type,
        "--voxel_size", str(voxel_size),
        "--output", output,
        "--classes", classes,
        "--cores", cores,
        "--models", models,
        "--n_breaks", n_breaks,
        "--resolution", resolution,
    ]

    if subunit is not None:
        command.extend(["--subunit", subunit])

    if selection is not None:
        command.extend(["--selection", selection])

    if return_spread:
        command.append("--return_spread")

    if save_log:
        if log_dir is None:
            log_dir = os.path.join(os.getcwd(), 'logs')
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.join(log_dir, 'prism_annotate.log')
        command.append(f"> {log_path}")

    if logger is not None:
        logger.info("Running prism annotate with command:")
        logger.info(" ".join(map(str, command)))

    os.system(" ".join(map(str, command)))

def prism_color(
    script_path: str,
    input: str,
    frame_index: int,
    subunit: str | None,
    resolution: int,
    selection: None,
    annotations_file: str,
    output: str,
    logger: logging.Logger | None = None,
):
    """ Run the script `color_precision.py` from prism
    - This colors the bead models based on the annotations generated
    - The output is a rmf3 file with the colored bead models

    Args:
        script_path (str): Path to the `color_precision.py` script from prism
        input (str): Path to the input rmf3 file
        frame_index (int): Frame index of the model to color from the rmf3 file
        subunit (str | None): name of the subunit to analyze (None for all)
        resolution (int): Resolution as number of residues per bead
        selection (None): selection within the subunit (None for all)
        annotations_file (str): Path to the annotations file from prism
        output (str): Path to the output rmf3 file
        logger (logging.Logger): Logger for logging messages.
    """

    command = [
        "python", script_path,
        "--input", input,
        "--frame_index", frame_index,
        "--resolution", resolution,
        "--annotations_file", annotations_file,
        "--output", output
    ]

    if subunit is not None:
        command.extend(["--subunit", subunit])

    if selection is not None:
        command.extend(["--selection", selection])

    logger.info("Running prism color with command:")
    logger.info(" ".join(map(str, command)))
    os.system(" ".join(map(str, command)))

def extract_sampcon(
    script_path: str,
    rmf1: str,
    list1: str,
    rmf2: str,
    list2: str,
    rmf_out: str,
    logger: logging.Logger | None = None,
):
    """ Run the script `extract_sampcon.py` from imp-sampcon

    - This extracts the models from the major cluster obtained from
        `exhaust.py` using the rmf3 and referring to the txt files from
        `run_extract_models.py`
    - The output is a single rmf3 file containing the extracted models

    Args:
        script_path (str): Path to the `extract_sampcon.py` script
        rmf1 (str): path to the rmf3 file for sample A
        list1 (str): path to the txt file for sample A
        rmf2 (str): path to the rmf3 file for sample B
        list2 (str): path to the txt file for sample B
        rmf_out (str): path to the output rmf3 file
        logger (logging.Logger | None, optional): Logger for logging messages.
        save_log (bool, optional): Whether to save the log file. Defaults to True.
    """

    command = [
        "python", script_path,
        "--rmf_out", rmf_out,
        "--rmf1", rmf1,
        "--list1", list1,
        "--rmf2", rmf2,
        "--list2", list2
    ]

    if logger is None:
        logger.info("Running extract_sampcon with command:")
        logger.info(" ".join(map(str, command)))

    os.system(" ".join(map(str, command)))

def rmf_to_xyzr(
    script_path: str,
    rmf_path: str,
    output_path: str,
    frame_subset: str | None = None,
    float_dtype: int = 64,
    nproc: int = 24,
    logger: logging.Logger | None = None,
):
    """ Run the script `rmf_to_xyzr.py`

    - This extracts the bead coordinates and radii from the input rmf3 file
    - The output is a hdf5 file containing the XYZR data for each molecule
      organized in a dictionary
    - key: molecule name with copy index (e.g. "mol1_0")
    - value: dictionary of fragments with their XYZR data
      (e.g. {"1-10": [[x, y, z, r], ...], "11": [[x, y, z, r], ...]})

    Args:
        script_path (str): Path to the `rmf_to_xyzr.py` script
        rmf_path (str): Path to the input rmf3 file
        output_path (str): Path to the output hdf5 file
        frame_subset (str | None, optional): Subset of frames to process
            (e.g., "0-9,17" for first 10 and 16th frame).
            Defaults to None (all frames).
        nproc (int, optional): Number of cores to use. Defaults to 24.
        logger (logging.Logger | None, optional): Logger for logging messages.
    """

    command = [
        "python", script_path,
        "--rmf_path", rmf_path,
        "--output_path", output_path,
        "--float_dtype", float_dtype,
        "--nproc", nproc,
    ]

    if frame_subset is not None:
        command.extend(["--frame_subset", frame_subset])

    if logger is not None:
        logger.info("Running rmf_to_xyzr with command:")
        logger.info(" ".join(map(str, command)))

    os.system(" ".join(map(str, command)))

def contact_map(
    script_path: str,
    xyzr_file: str,
    contact_map_dir: str,
    nproc: int = 24,
    dist_cutoff: float = 10.0,
    frac_cutoff: float = 0.25,
    plotting: str = "matplotlib",
    merge_copies: bool = False,
    float_dtype: int = 64,
    int_dtype: int = 32,
    logger: logging.Logger | None = None,
):
    """ Run the script `contact_map.py` and get pairwise distance/contact maps

    Args:
        script_path (str): Path to the `contact_map.py` script
        xyzr_file (str): Path to the input hdf5 file containing XYZR data
        nproc (int): Number of cores to use
        dist_cutoff (float): dist_cutoff distance for contact map (in Angstroms)
        frac_cutoff (float): Minimum fraction of frames for a contact
        contact_map_dir (str): Directory to save contact map outputs
        plotting (str): Type of plotting to perform ('matplotlib', 'plotly')
        merge_copies (bool): Whether to merge maps across copies for protein
            pairs
        logger (logging.Logger | None, optional): Logger for logging messages.
    """

    command = [
        "python", script_path,
        "--xyzr_file", xyzr_file,
        "--nproc", nproc,
        "--dist_cutoff", dist_cutoff,
        "--frac_cutoff", frac_cutoff,
        "--contact_map_dir", contact_map_dir,
        "--plotting", plotting,
        "--float_dtype", float_dtype,
        "--int_dtype", int_dtype,
    ]

    if merge_copies:
        command.append("--merge_copies")

    if logger is not None:
        logger.info("Running contact_map with command:")
        logger.info(" ".join(map(str, command)))

    os.system(" ".join(map(str, command)))

def fit_to_binding_data(
    script_path: str,
    xyzr_file: str,
    input_config: str,
    output_dir: str,
    nproc: int = 24,
    merge_copies: bool = True,
    float_dtype: int = 64,
    logger: logging.Logger | None = None,
):
    """ Run the script `fit_to_binding_data.py` and fit the models to binding data

    Args:
        script_path (str): Path to the `fit_to_binding_data.py` script
        xyzr_file (str): Path to the input hdf5 file containing XYZR data
        input_config (str): Path to the input config file containing binding data
        output_dir (str): Directory to save outputs
        nproc (int): Number of cores to use
        merge_copies (bool): Whether to merge maps across copies for protein
            pairs
        logger (logging.Logger | None, optional): Logger for logging messages.
    """

    command = [
        "python", script_path,
        "--xyzr_file", xyzr_file,
        "--input", input_config,
        "--output_dir", output_dir,
        "--nproc", nproc,
        "--float_dtype", float_dtype,
    ]

    if merge_copies:
        command.append("--merge_copies")

    if logger is not None:
        logger.info("Running fit_to_binding_data with command:")
        logger.info(" ".join(map(str, command)))

    os.system(" ".join(map(str, command)))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Scripts to analyze modeling output in an end-to-end manner"
    )
    parser.add_argument(
        "--analysis_dir",
        type=str,
        default=f"/data/{_user}/imp_toolbox_test/analysis",
        help="Path to the analysis output directory"
    )
    parser.add_argument(
        "--modeling_dir",
        type=str,
        default=f"/data/{_user}/imp_toolbox_test/modeling",
        help="Path to the modeling output directory"
    )
    parser.add_argument(
        "--keep_logs",
        action='store_true',
        help="Whether to keep intermediate log files (default: False)"
    )
    parser.add_argument(
        "--compound_log_mode",
        type=str,
        default="a+",
        choices=["a+", "w"],
        help="Mode for compound log file (default: a+)"
    )
    parser.add_argument(
        "--scripts_to_run",
        nargs='+',
        default=[
            "run_analysis_trajectories",
            "run_extract_models",
            "exhaust",
            "extract_sampcon",
            "prism_annotate",
            "prism_color",
            "rmf_to_xyzr",
        ],
        help="List of scripts to run in sequence (default: all) \
            (default: [run_analysis_trajectories, run_extract_models \
            exhaust, extract_sampcon, prism_annotate, prism_color])"
    )
    args = parser.parse_args()

    assert all([s in [
        "run_analysis_trajectories",
        "variable_filter",
        "run_extract_models",
        "exhaust",
        "extract_sampcon",
        "prism_annotate",
        "prism_color",
        "rmf_to_xyzr",
        "contact_map",
    ] for s in args.scripts_to_run]), (
        f"""
        Invalid script name in scripts_to_run.
        Valid options are:
        run_analysis_trajectories
        variable_filter
        run_extract_models
        exhaust
        extract_sampcon
        prism_annotate
        prism_color
        rmf_to_xyzr
        contact_map
        """
    )

    ###########################################################################

    start_t = time.perf_counter()

    ANALYSIS_OUTPUT_PATH = args.analysis_dir
    modeling_dir = args.modeling_dir
    LOG_DIR = os.path.join(ANALYSIS_OUTPUT_PATH, "logs")

    assert os.path.exists(modeling_dir), (
        f"""
        Modeling output path {modeling_dir} does not exist.
        Please check run modeling script first.
        """
    )

    os.makedirs(LOG_DIR, exist_ok=True)
    logger = logging.getLogger(__name__)
    file_handler = logging.FileHandler(
        os.path.join(LOG_DIR, "end_to_end_analysis.log"),
        mode=args.compound_log_mode
    )
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    logger.info("\nStarting end-to-end analysis...")

    os.makedirs(ANALYSIS_OUTPUT_PATH, exist_ok=True)

    ###########################################################################
    # run analysis trajectories
    ###########################################################################
    pmi_analysis_output_path = os.path.join(
        ANALYSIS_OUTPUT_PATH, "pmi_analysis"
    )
    os.makedirs(pmi_analysis_output_path, exist_ok=True)

    if "run_analysis_trajectories" in args.scripts_to_run:

        run_analysis_trajectories(
            script_path=f"{IMP_TOOLOBX_PATH}/analysis/run_analysis_trajectories.py",
            modeling_dir=modeling_dir,
            analysis_output_path=pmi_analysis_output_path,
            traj_dir_prefix=TRAJ_DIR_PREFIX,
            run_start=1,
            run_end=2,
            run_interval=1,
            nproc=4,
            burn_in_fraction=0.1,
            nskip=1,
            restraint_handles=RESTRAINT_HANDLES,
            hdbscan_restraint_handles=HDBSCAN_RESTRAINT_HANDLES,
            min_cluster_size=150,
            min_samples=5,
            logger=logger,
            save_log=args.keep_logs
        )
        lap = time.perf_counter()
        logger.info(
            f"Completed run_analysis_trajectories in {lap - start_t:0.4f} seconds"
        )

    else:
        logger.info("Skipping run_analysis_trajectories as per user request.")

    ###########################################################################
    # get major cluster
    ###########################################################################
    hdbscan_summary_path = os.path.join(
        pmi_analysis_output_path, 'summary_hdbscan_clustering.dat',
    )
    assert os.path.exists(hdbscan_summary_path), (
        f"""
        HDBSCAN summary file {hdbscan_summary_path} does not exist.
        Please check if run_analysis_trajectories has been run successfully.
        And there is at least one cluster found.
        """
    )
    major_cluster_idx, major_cluster_size = return_major_cluster(
        hdbscan_log_path=hdbscan_summary_path
    )
    pmi_cluster_idx = major_cluster_idx
    logger.info(
        f"""
        Major cluster from PMI_analysis is {major_cluster_idx}
        with {major_cluster_size} models.
        """
    )

    ###########################################################################
    # variable filter
    ###########################################################################

    if "variable_filter" in args.scripts_to_run:

        if major_cluster_size <= MODEL_CAP:
            logger.info(
                f"""
                Major cluster size {major_cluster_size} is less than
                model cap {MODEL_CAP}. Skipping variable filter.
                """
            )

        else:
            logger.info(
                f"""
                Major cluster size {major_cluster_size} is greater than
                model cap {MODEL_CAP}. Proceeding with variable filter.
                """
            )
            var_filter_output_path = os.path.join(
                ANALYSIS_OUTPUT_PATH, 'variable_filter_output'
            )
            os.makedirs(var_filter_output_path, exist_ok=True)
            variable_filter(
                script_path=f"{IMP_TOOLOBX_PATH}/analysis/variable_filter.py",
                major_cluster_idx=pmi_cluster_idx,
                lowest_cutoff=-2.0,
                highest_cutoff=3.0,
                step_size=0.01,
                model_cap=MODEL_CAP,
                gsmsel_dir=pmi_analysis_output_path,
                output_dir=var_filter_output_path,
                restraint_handles=HDBSCAN_RESTRAINT_HANDLES,
                logger=logger,
                save_log=args.keep_logs
            )
            lap = time.perf_counter()
            logger.info(
                f"Completed variable_filter in {lap - start_t:0.4f} seconds"
            )

    else:
        logger.info("Skipping variable_filter as per user request.")
        var_filter_output_path = None

    ###########################################################################
    # extract models
    ###########################################################################
    if "run_extract_models" in args.scripts_to_run:

        run_extract_models(
            script_path=f"{IMP_TOOLOBX_PATH}/analysis/run_extract_models.py",
            modeling_dir=modeling_dir,
            analysis_output_path=pmi_analysis_output_path,
            traj_dir_prefix=TRAJ_DIR_PREFIX,
            run_start=1,
            run_end=4,
            run_interval=1,
            nproc=4,
            burn_in_fraction=0.1,
            nskip=1,
            cluster_id=pmi_cluster_idx,
            filter_applied=major_cluster_size > MODEL_CAP,
            variable_filter_output_dir=var_filter_output_path,
            logger=logger,
            save_log=args.keep_logs
        )
        lap = time.perf_counter()
        logger.info(
            f"Completed run_extract_models in {lap - start_t:0.4f} seconds"
        )
    else:
        logger.info("Skipping run_extract_models as per user request.")

    ###########################################################################
    # exhaust
    ###########################################################################
    sampcon_output_dir = os.path.join(
        ANALYSIS_OUTPUT_PATH, 'sampcon_output'
    )

    if "exhaust" in args.scripts_to_run:

        os.makedirs(sampcon_output_dir, exist_ok=True)
        _current_dir = os.getcwd()
        os.chdir(sampcon_output_dir)
        exhaust(
            script_path=f"{SAMPCON_PATH}/pyext/src/exhaust.py",
            pmi_analysis_output_path=pmi_analysis_output_path,
            sysname=SYSNAME,
            scoreA=f"A_models_clust{str(pmi_cluster_idx)}.txt",
            scoreB=f"B_models_clust{str(pmi_cluster_idx)}.txt",
            rmfA=f"A_models_clust{str(pmi_cluster_idx)}.rmf3",
            rmfB=f"B_models_clust{str(pmi_cluster_idx)}.rmf3",
            density=SAMPCON_DENSITY_TXT,
            gnuplot=True,
            prism=True,
            align=True,
            mode="cpu_omp",
            matrix_cores=12,
            cluster_cores=4,
            gridsize="5",
            logger=logger,
            save_log=args.keep_logs
        )
        os.chdir(_current_dir)
        lap = time.perf_counter()
        logger.info(f"Completed exhaust in {lap - start_t:0.4f} seconds")

    else:
        logger.info("Skipping exhaust as per user request.")

    sampcon_cluster_idx = 0 # only analyze the cluster 0 from exhaust.py

    ###########################################################################
    # compare cluster sample LPDs for sample A and B
    ###########################################################################
    if "correlate_cluster_sample_LPDs" in args.scripts_to_run:
        sampcon_cluster_path=os.path.join(
            sampcon_output_dir, f"cluster.{sampcon_cluster_idx}"
        )
        assert os.path.exists(sampcon_cluster_path), (
            f"""Sampcon cluster path {sampcon_cluster_path} does not exist.
            Please check if exhaust has been run successfully and 
            cluster.{sampcon_cluster_idx} directory exists. """
        )

        correlate_cluster_sample_LPDs(
            script_path=f"{IMP_TOOLOBX_PATH}/analysis/correlate_cluster_sample_LPDs.py",
            sampcon_cluster_path=sampcon_cluster_path,
            mode="pasani",
            use_combined_map=True,
            logger=logger,
            save_log=args.keep_logs
        )
        lap = time.perf_counter()
        logger.info(
            f"Completed correlation of cluster sample LPDs in {lap - start_t:0.4f} seconds"
        )

    ###########################################################################
    # extract sampcon frames
    ###########################################################################
    extracted_rmf_path = os.path.join(
        sampcon_output_dir, "sampcon_extracted_frames.rmf3"
    )
    frame_ids_A_txt = os.path.join(
        sampcon_output_dir, f"cluster.{sampcon_cluster_idx}.sample_A.txt"
    )
    frame_ids_B_txt = os.path.join(
        sampcon_output_dir, f"cluster.{sampcon_cluster_idx}.sample_B.txt"
    )
    frames_A_rmf = os.path.join(
        pmi_analysis_output_path, f"A_models_clust{str(pmi_cluster_idx)}.rmf3"
    )
    frames_B_rmf = os.path.join(
        pmi_analysis_output_path, f"B_models_clust{str(pmi_cluster_idx)}.rmf3"
    )

    if "extract_sampcon" in args.scripts_to_run:
        assert (
            os.path.exists(frame_ids_A_txt) and os.path.exists(frame_ids_B_txt)
        ), (
            f"""Sample A/B txt files does not exist in {sampcon_output_dir}.
            Please check if exhaust has been run successfully. """
        )
        extract_sampcon(
            script_path=f"{IMP_TOOLOBX_PATH}/analysis/xtract_sampcon.py",
            rmf1=frames_A_rmf,
            list1=frame_ids_A_txt,
            rmf2=frames_B_rmf,
            list2=frame_ids_B_txt,
            rmf_out=extracted_rmf_path,
            logger=logger
        )
        assert os.path.exists(extracted_rmf_path), (
            f"""Extracted RMF file {extracted_rmf_path} does not exist.
            Please check if extract_sampcon has been run successfully.
            """
        )
        lap = time.perf_counter()
        logger.info(f"Completed extract_sampcon in {lap - start_t:0.4f} seconds")

    else:
        logger.info("Skipping extract_sampcon as per user request.")

    ###########################################################################
    # prism
    ###########################################################################
    prism_output_dir = os.path.join(
        ANALYSIS_OUTPUT_PATH, 'prism_output'
    )
    prism_input_path = os.path.join(
        sampcon_output_dir, f"cluster.{sampcon_cluster_idx}.prism.npz"
    )
    prism_annotation_path = os.path.join(
        prism_output_dir, f"annotations_cl2.txt"
    )

    if "prism_annotate" in args.scripts_to_run:
        os.makedirs(prism_output_dir, exist_ok=True)
        prism_annotate(
            script_path=f"{PRISM_PATH}/src/main.py",
            input=prism_input_path,
            input_type="npz",
            voxel_size=4,
            return_spread=True,
            output=prism_output_dir,
            classes=2,
            cores=16,
            models=1.0,
            n_breaks=50,
            resolution=1,
            subunit=None,
            selection=None,
            logger=logger,
            save_log=args.keep_logs
        )
        lap = time.perf_counter()
        logger.info(f"Completed prism in {lap - start_t:0.4f} seconds")

    else:
        logger.info("Skipping prism as per user request.")

    if "prism_color" in args.scripts_to_run:
        assert os.path.exists(prism_annotation_path), (
            f"""Prism annotation file {prism_annotation_path} does not exist.
            Please check if prism has been run successfully.
            """
        )
        cluster_center_rmf_path = os.path.join(
            sampcon_output_dir, f"cluster.{sampcon_cluster_idx}", "cluster_center_model.rmf3"
        )
        prism_cluster_center_rmf_path = cluster_center_rmf_path.replace(
            ".rmf3", "_prism_colored.rmf3"
        )
        prism_color(
            script_path=f"{PRISM_PATH}/src/color_precision.py",
            input=cluster_center_rmf_path,
            frame_index=0,
            subunit=None,
            resolution=1,
            selection=None,
            annotations_file=prism_annotation_path,
            output=prism_cluster_center_rmf_path,
            logger=logger
        )

        lap = time.perf_counter()
        logger.info(f"Completed prism_color in {lap - start_t:0.4f} seconds")

    else:
        logger.info("Skipping prism_color as per user request.")

    ###########################################################################
    # align pdb to cluster center rmf
    ###########################################################################

    if "align_pdb_to_ccm" in args.scripts_to_run:

        assert os.path.exists(cluster_center_rmf_path), (
            f"""Cluster center RMF file {cluster_center_rmf_path}
            does not exist."""
        )

        fit_pdb_to_ccm(
            script_path=f"{IMP_TOOLOBX_PATH}/analysis/align_pdb_to_ccm.py",
            ccm_file=cluster_center_rmf_path,
            input_config=f"{IMP_TOOLOBX_PATH}/analysis/fits_to_perform.json",
            output_dir=os.path.join(ANALYSIS_OUTPUT_PATH, "aligned_pdb_to_ccm"),
            logger=logger
        )

        lap = time.perf_counter()
        logger.info(f"Completed align_pdb_to_ccm in {lap - start_t:0.4f} seconds")

    else:
        logger.info("Skipping align_pdb_to_ccm as per user request.")

    ###########################################################################
    # rmf to xyzr
    ###########################################################################

    xyzr_output_path = os.path.join(
        ANALYSIS_OUTPUT_PATH, "sampcon_extracted_frames_xyzr.h5"
    )
    if "rmf_to_xyzr" in args.scripts_to_run:

        assert os.path.exists(extracted_rmf_path), (
            f"""Extracted RMF file {extracted_rmf_path} does not exist.
            Please check if extract_sampcon has been run successfully.
            """
        )
        rmf_to_xyzr(
            script_path=f"{IMP_TOOLOBX_PATH}/analysis/rmf_to_xyzr1.py",
            rmf_path=extracted_rmf_path,
            output_path=xyzr_output_path,
            frame_subset=None,
            float_dtype=64,
            nproc=24,
            logger=logger
        )
        assert os.path.exists(xyzr_output_path), (
            f"""XYZR output file {xyzr_output_path} does not exist.
            Please check if rmf_to_xyzr has been run successfully.
            """
        )
        lap = time.perf_counter()
        logger.info(f"Completed rmf_to_xyzr in {lap - start_t:0.4f} seconds")

    else:
        logger.info("Skipping rmf_to_xyzr as per user request.")

    ###########################################################################
    # contact map
    ###########################################################################

    if "contact_map" in args.scripts_to_run:

        contact_map_dir = os.path.join(ANALYSIS_OUTPUT_PATH, "contact_map")
        os.makedirs(contact_map_dir, exist_ok=True)
        assert os.path.exists(xyzr_output_path), (
            f"""XYZR output file {xyzr_output_path} does not exist.
            Please check if rmf_to_xyzr has been run successfully.
            """
        )
        contact_map(
            script_path=f"{IMP_TOOLOBX_PATH}/analysis/contact_map.py",
            xyzr_file=xyzr_output_path,
            contact_map_dir=contact_map_dir,
            nproc=24,
            dist_cutoff=10.0,
            frac_cutoff=0.25,
            plotting="matplotlib",
            merge_copies=False,
            float_dtype=64,
            int_dtype=32,
            logger=logger,
        )
        lap = time.perf_counter()
        logger.info(f"Completed contact_map in {lap - start_t:0.4f} seconds")

    else:
        logger.info("Skipping contact_map as per user request.")

    ###########################################################################
    # Fit to binding data (optional)
    ###########################################################################

    if "fit_to_binding_data" in args.scripts_to_run:

        fit_to_binding_data_dir = os.path.join(
            ANALYSIS_OUTPUT_PATH, "fit_to_binding_data"
        )
        os.makedirs(fit_to_binding_data_dir, exist_ok=True)
        input_config = f"{IMP_TOOLOBX_PATH}/analysis/fit_to_binding.json"
        assert os.path.exists(xyzr_output_path), (
            f"""XYZR output file {xyzr_output_path} does not exist.
            Please check if rmf_to_xyzr has been run successfully.
            """
        )
        assert os.path.exists(input_config), (
            f"""Input config file {input_config} does not exist.
            Please check if the file exists at the specified path. 
            """
        )
        fit_to_binding_data(
            script_path=f"{IMP_TOOLOBX_PATH}/analysis/fit_to_binding_data.py",
            xyzr_file=xyzr_output_path,
            input_config=input_config,
            output_dir=fit_to_binding_data_dir,
            nproc=24,
            merge_copies=True,
            float_dtype=64,
            logger=logger
        )
        lap = time.perf_counter()
        logger.info(f"Completed fit_to_binding_data in {lap - start_t:0.4f} seconds")

    ###########################################################################
    logger.info("End-to-end analysis completed.")
    logger.info("Ran following scripts:")
    for script in args.scripts_to_run:
        logger.info(f" - {script}")
    end_time = time.perf_counter()
    logger.info(f"Total time: {end_time - start_t:0.4f} seconds")