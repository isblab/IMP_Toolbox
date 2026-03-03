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
IMP_TOOLOBX_PATH = "/home/$USER/Projects/IMP_TOOLBOX/IMP_Toolbox"
SYSNAME = "cardiac_desmosome"
SAMPCON_DENSITY_TXT = f"/data/{_user}/imp_toolbox_test/input/density_sampcon.txt"
MODEL_CAP = 30000 # for variable filter

random.seed(47)

def return_major_cluster(hdbscan_log_path: str):
    """ Return the cluster with the maximum number of models from HDBSCAN
    clustering summary file.

    ## Arguments:

    - **hdbscan_log_path (str)**:<br />
        Path to the HDBSCAN clustering summary file.

    ## Returns:

    - **tuple**:<br />
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

    ## Arguments:

    - **script_path (str)**:<br />
        Path to the `run_analysis_trajectories.py` script.

    - **modeling_dir (str)**:<br />
        Path to the modeling output directory containing the trajectory folders.

    - **analysis_dir (str)**:<br />
        Path to the analysis output directory where the results will be stored.

    - **traj_dir_prefix (str)**:<br />
        Prefix for trajectory directories (e.g., 'run_') in the modeling output directory.

    - **run_start (int)**:<br />
        Starting run number to analyze (e.g., 1 for 'run_1').

    - **run_end (int)**:<br />
        Ending run number to analyze (e.g., 50 for 'run_50').

    - **run_interval (int)**:<br />
        Interval between runs to analyze (e.g., 1 for every run, 2 for every other run).

    - **nproc (int)**:<br />
        Number of cores to use for analysis.

    - **burn_in_fraction (float)**:<br />
        Fraction of data to discard as burn-in. These are fraction of frames from the start of the trajectory to discard (e.g., 0.1 for 10% burn-in).

    - **nskip (int)**:<br />
        Number of consecutive frames to skip for analysis (e.g., 0 for no skipping, 1 to analyze every other frame).

    - **restraint_handles (list)**:<br />
        List of restraint handles to analyze (e.g., ['GaussianEMRestraint:EM', 'SingleAxisMinGaussianRestraint:SAMGR']).

    - **hdbscan_restraint_handles (list)**:<br />
        List of restraint handles to use for HDBSCAN clustering (e.g., ['EV_sum', 'EM_sum']).

    - **min_cluster_size (int)**:<br />
        Minimum number of samples in a cluster for HDBSCAN.

    - **min_samples (int)**:<br />
        Minimum number of samples in a cluster for HDBSCAN.

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages. If None, logging is skipped. Defaults to None.

    - **save_log (bool, optional):**:<br />
        Whether to save the log file.

    - **log_dir (str, optional):**:<br />
        Directory where the log file is saved.
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

    ## Arguments:

    - **script_path (str)**:<br />
        Path to the `variable_filter.py` script.

    - **pmi_clust_idx (int)**:<br />
        Index of the major cluster from HDBSCAN clustering to filter models from.

    - **lowest_cutoff (float)**:<br />
        Standard deviation multiplier for the most stringent cutoff. Models with restraint scores above this cutoff will be retained.

    - **highest_cutoff (float)**:<br />
        Standard deviation multiplier for the most lenient cutoff. Models with restraint scores above this cutoff will be retained.

    - **step_size (float)**:<br />
        Step size for the cutoff. The script will apply cutoffs starting from `lowest_cutoff` to `highest_cutoff` with this step size until the number of models retained is less than or equal to `model_cap`.

    - **model_cap (int)**:<br />
        Maximum number of models to retain after filtering. The script will apply the cutoffs iteratively until the number of models retained is less than or equal to this cap.

    - **gsmsel_dir (str)**:<br />
        Directory where the pmi_analysis csv files are stored. These csv files contain the restraint scores for each model and are used for filtering.

    - **output_dir (str)**:<br />
        Path to the directory where the output of the script will be stored. The output includes the filtered rmf3 file and a txt file indicating the scores of the retained models.

    - **restraint_handles (list)**:<br />
        List of restraint handles to consider for filtering. The script will use the scores from these restraint handles for applying the cutoffs.

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages. If None, logging is skipped. Defaults to None.

    - **save_log (bool, optional):**:<br />
        Whether to save the log file.

    - **log_dir (str, optional):**:<br />
        Directory where the log file will be saved.
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

    ## Arguments:

    - **script_path (str)**:<br />
        Path to the `run_extract_models.py` script.

    - **modeling_dir (str)**:<br />
        Path to the modeling output directory containing the trajectory folders.

    - **analysis_dir (str)**:<br />
        Path to the analysis output directory where the extracted models will be stored.

    - **traj_dir_prefix (str)**:<br />
        Prefix for trajectory directories (e.g., 'run_')

    - **run_start (int)**:<br />
        Starting run number

    - **run_end (int)**:<br />
        Ending run number

    - **run_interval (int)**:<br />
        Interval between runs (e.g., 1 for every run)

    - **nproc (int)**:<br />
        Number of cores to use for analysis

    - **burn_in_fraction (float)**:<br />
        Fraction of data to discard as burn-in

    - **nskip (int)**:<br />
        Number of consecutive frames to skip (e.g., 0 for no skipping, 1 to analyze every other frame)

    - **cluster_id (int)**:<br />
        Cluster ID to extract models from for the major cluster from the HDBSCAN clustering. This is used to identify the csv files containing the scores of the models to be extracted.

    - **filter_applied (bool)**:<br />
        Whether `variable_filter.py` was run or not. This is used to determine whether to look for the csv files from `variable_filter.py` or from `run_analysis_trajectories.py` for extracting the models.

    - **variable_filter_output_dir (str)**:<br />
        Path to the output directory of `variable_filter.py` if `filter_applied` is True. This is used to identify the csv files containing the scores of the models to be extracted. If `filter_applied` is False, this argument is ignored and the script looks for the csv files from `run_analysis_trajectories.py` in the `analysis_dir`.

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages. If None, logging is skipped. Defaults to None.

    - **save_log (bool, optional):**:<br />
        Whether to save the log file.

    - **log_dir (str, optional):**:<br />
        Path to the directory where log files are saved. If None, a default path is used.
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
    voxel_size: float,
    selection: str | None = None,
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

    ## Arguments:

    - **script_path (str)**:<br />
        Path to the `exhaust.py` script from imp-sampcon.

    - **pmi_analysis_dir (str)**:<br />
        Path to the pmi_analysis output directory containing the scores and rmf3 files for sample A and B. This is used as input for the exhaustiveness analysis.

    - **sysname (str)**:<br />
        System name for the analysis (e.g., 'cardiac_desmosome'). This is used as a prefix for the output files generated from the analysis.

    - **scoreA (str)**:<br />
        txt file containing scores for sample A. This is used as input for the exhaustiveness analysis.

    - **scoreB (str)**:<br />
        txt file containing scores for sample B. This is used as input for the exhaustiveness analysis.

    - **rmfA (str)**:<br />
        rmf3 file containing models for sample A. This is used as input for the exhaustiveness analysis.

    - **rmfB (str)**:<br />
        rmf3 file containing models for sample B. This is used as input for the exhaustiveness analysis.

    - **density (str)**:<br />
        Path to the density txt file for imp-sampcon. This file contains the density of models in the score space and is used for calculating the localization probability densities.

    - **gnuplot (bool)**:<br />
        Whether to generate gnuplot files for visualizing the score distributions and convergence.

    - **prism (bool)**:<br />
        Whether to generate prism input `.npz` files for each cluster. These files can be used for further analysis in prism.

    - **align (bool)**:<br />
        Whether to align the models before analysis. If True, the models will be aligned to a reference structure before calculating the localization probability densities.

    - **mode (str)**:<br />
        Mode for imp-sampcon (e.g., 'cpu_omp'). This is used to specify the mode for running the analysis.

    - **matrix_cores (int)**:<br />
        Number of cores to use for RMSD matrix calculation.

    - **cluster_cores (int)**:<br />
        Number of cores to use for clustering. Do not set too high, might lead
        to memory issues.

    - **gridsize (int)**:<br />
        Grid size for clustering. This is used to specify the grid size for calculating the localization probability densities.

    - **voxel_size (float)**:<br />
        Voxel size for calculating localization probability densities. This is used to specify the resolution of the density maps generated from the analysis.

    - **selection (str | None, optional):**:<br />
        Selection string for analyzing a subset of the model (e.g., "chain A and resnum 1-100"). If None, the entire model is analyzed.

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages. If None, logging is skipped. Defaults to None.

    - **save_log (bool, optional):**:<br />
        Whether to save the log output to a file.

    - **log_dir (str, optional):**:<br />
        Directory where the log file is saved. If None, the log is saved in the current working directory.
    """

    command = [
        "python", script_path,
        "--sysname", sysname,
        "--mode", mode,
        "--matrix-cores", matrix_cores,
        "--cluster-cores", cluster_cores,
        "--density", density,
        "--gridsize", gridsize,
        "--voxel", voxel_size,
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
    if selection is not None:
        command.append("--selection")
        command.append(selection)

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

def correlate_cluster_sample_densities(
    script_path: str,
    sampcon_cluster_path: str,
    mode: str,
    use_combined_map: bool,
    logger: logging.Logger | None = None,
    save_log: bool = True,
    log_dir: str = None,
):
    """ Run the script `correlate_cluster_sample_densities.py` from IMP_Toolbox

    ## Arguments:

    - **script_path (str)**:<br />
        Path to the `correlate_cluster_sample_densities.py` script.

    - **sampcon_cluster_path (str)**:<br />
        Path to the directory containing the cluster directories generated from
        `exhaust.py` in imp-sampcon. Each cluster directory should contain the
        localization probability density maps for sample A and B.

    - **mode (str)**:<br />
        Mode for imp-sampcon (e.g., 'cpu_omp'). This is used to specify the mode
        for running the analysis.

    - **use_combined_map (bool)**:<br />
        Whether to use the combined map for correlation analysis.

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages. If None, logging is skipped. Defaults to None.

    - **save_log (bool, optional):**:<br />
        Whether to save the log output to a file.

    - **log_dir (str, optional):**:<br />
        Directory where the log file is saved. If None, the log is saved in the
        current working directory.
    """

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
    """ Run the script `align_pdb_to_ccm.py` from IMP_Toolbox.

    ## Arguments:

    - **script_path (str)**:<br />
        Path to the `align_pdb_to_ccm.py` script.

    - **ccm_file (str)**:<br />
        Path to the input CCM file. This file contains the contact information
        that will be used for fitting the PDB structure.

    - **input_config (str)**:<br />
        Path to the input configuration file. This file contains the parameters
        for fitting the PDB structure to the CCM.

    - **output_dir (str)**:<br />
        Path to the output directory where the fitted PDB structure and related
        files will be saved.

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages. If None, logging is skipped.
    """

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

    ## Arguments:

    - **script_path (str)**:<br />
        Path to the `main.py` script from prism.

    - **input (str)**:<br />
        Path to the input file for annotation. These are ".npz" files generated
        from `exhaust.py` in imp-sampcon.

    - **input_type (str)**:<br />
        Type of the input file (e.g., 'npz', 'rmf3'). This is used to specify
        the type of the input file for annotation.

    - **voxel_size (int)**:<br />
        Voxel size in Angstroms for calculating the precision. This is used to
        specify the resolution of the density maps generated for annotation.

    - **return_spread (bool)**:<br />
        Whether to return bead spread information. If True, the script will
        return the spread of the beads in addition to the precision annotation.

    - **output (str)**:<br />
        Path to the output file for the annotations. The output is a rmf3 file
        containing the bead models with precision annotations.

    - **classes (int)**:<br />
        Number of classes for annotation.

    - **cores (int)**:<br />
        Number of cores.

    - **models (float)**:<br />
        Fraction of models to use for annotation (between 0 and 1). This is used
        to specify the fraction of models to be used for generating the annotations.

    - **n_breaks (int)**:<br />
        Number of breaks for jenkspy.

    - **resolution (int)**:<br />
        Resolution as number of residues per bead. This is used to specify the
        resolution at which the precision annotation is generated.

    - **subunit (str | None)**:<br />
        Name of the subunit to analyze. If None, all subunits are analyzed.

    - **selection (None)**:<br />
        Selection within the subunit. If None, the entire subunit is analyzed.

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages. If None, logging is skipped. Defaults to None.

    - **save_log (bool, optional):**:<br />
        Whether to save the log output to a file.

    - **log_dir (str, optional):**:<br />
        Directory where the log file is saved.
        If None, the log is saved in the current working directory.
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

    ## Arguments:

    - **script_path (str)**:<br />
        Path to the `color_precision.py` script from prism.

    - **input (str)**:<br />
        Path to the input rmf3 file containing the bead models to be colored.

    - **frame_index (int)**:<br />
        Frame index of the model to color from the rmf3 file. This is used to
        specify which frame of the rmf3 file to use for coloring the models.

    - **subunit (str | None)**:<br />
        Name of the subunit to analyze. If None, all subunits are analyzed. This is
        used to specify which subunit's models to color based on the annotations.

    - **resolution (int)**:<br />
        Resolution as number of residues per bead. This is used to specify the
        resolution at which the models are colored based on the annotations.

    - **selection (None)**:<br />
        Selection within the subunit. If None, the entire subunit is analyzed. This
        is used to specify which part of the subunit's models to color based on the annotations.

    - **annotations_file (str)**:<br />
        Path to the annotations file generated from the `prism_annotate.py` script.

    - **output (str)**:<br />
        Path to the output rmf3 file with the colored bead models.

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages. If None, logging is skipped. Defaults to None.
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

    ## Arguments:

    - **script_path (str)**:<br />
        Path to the `extract_sampcon.py` script from imp-sampcon.

    - **rmf1 (str)**:<br />
        Path to the rmf3 file for sample A. This file contains the models for
        sample A and is used as input for extracting the models from the major cluster.

    - **list1 (str)**:<br />
        Path to the txt file for sample A. This file contains the list of models
        in the major cluster for sample A.

    - **rmf2 (str)**:<br />
        Path to the rmf3 file for sample B. This file contains the models for
        sample B and is used as input for extracting the models from the major cluster.

    - **list2 (str)**:<br />
        Path to the txt file for sample B. This file contains the list of models
        in the major cluster for sample B.

    - **rmf_out (str)**:<br />
        Path to the output rmf3 file. The output is a single rmf3 file containing
        the extracted models from the major cluster for both sample A and B.

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages. If None, logging is skipped. Defaults to None.
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
    overwrite: bool = False,
    logger: logging.Logger | None = None,
):
    """ Run the script `rmf_to_xyzr.py`

    - This extracts the bead coordinates and radii from the input rmf3 file
    - The output is a hdf5 file containing the XYZR data for each molecule
      organized in a dictionary
    - key: molecule name with copy index (e.g. "mol1_0")
    - value: dictionary of fragments with their XYZR data
      (e.g. {"1-10": [[x, y, z, r], ...], "11": [[x, y, z, r], ...]})

    ## Arguments:

    - **script_path (str)**:<br />
        Path to the `rmf_to_xyzr.py` script.

    - **rmf_path (str)**:<br />
        Path to the input rmf3 file containing the bead models.

    - **output_path (str)**:<br />
        Path to the output hdf5 file. The output is a hdf5 file containing the
        XYZR data for each molecule organized in a dictionary.

    - **frame_subset (str | None, optional):**:<br />
        Subset of frames to process (e.g., "0-9,17" for first 10 and 16th frame).
        If None, all frames are processed.

    - **float_dtype (int, optional):**:<br />
        Floating point precision to use for the XYZR data in the output hdf5 file (e.g., 32 or 64). Defaults to 64.

    - **nproc (int, optional):**:<br />
        Number of cores to use for processing.

    - **overwrite (bool, optional):**:<br />
        Whether to overwrite existing output files.

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages.
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

    if overwrite:
        command.append("--overwrite")

    if logger is not None:
        logger.info("Running rmf_to_xyzr with command:")
        logger.info(" ".join(map(str, command)))

    os.system(" ".join(map(str, command)))

def interaction_map(
    script_path: str,
    xyzr_file: str,
    interaction_map_dir: str,
    nproc: int = 24,
    dist_cutoff: float = 10.0,
    frac_cutoff: float = 0.25,
    plotting: str = "matplotlib",
    merge_copies: bool = False,
    binarize_cmap: bool = False,
    binarize_dmap: bool = False,
    float_dtype: int = 64,
    int_dtype: int = 32,
    overwrite: bool = False,
    logger: logging.Logger | None = None,
):
    """ Run the script `interaction_map.py` and get pairwise distance/contact maps

    ## Arguments:

    - **script_path (str)**:<br />
        Path to the `interaction_map.py` script.

    - **xyzr_file (str)**:<br />
        Path to the input hdf5 file containing XYZR data. This file is generated
        from the `rmf_to_xyzr.py` script and contains the coordinates and radii of
        the beads for each molecule in the system.

    - **interaction_map_dir (str)**:<br />
        Directory to save the outputs from the interaction map analysis.
        The outputs include the pairwise distance maps and contact maps for each
        frame.

    - **nproc (int, optional):**:<br />
        Number of cores to use for parallel processing. Defaults to 24.

    - **dist_cutoff (float, optional):**:<br />
        Distance cutoff (in Angstroms) for defining contacts in the contact map.
        Beads that are within this distance cutoff are considered to be in contact.
        Defaults to 10.0 Angstroms.

    - **frac_cutoff (float, optional):**:<br />
        Minimum fraction of frames for a contact to be included in the final contact map.
        For example, if frac_cutoff is 0.25, then only contacts that are
        present in at least 25% of the frames will be included in the final contact map.
        Defaults to 0.25.

    - **plotting (str, optional):**:<br />
        Type of plotting to perform for visualizing the interaction maps.
        Options include 'matplotlib' for static plots and 'plotly' for interactive plots.
        Defaults to 'matplotlib'.

    - **merge_copies (bool, optional):**:<br />
        Whether to merge maps across copies for protein pairs. If True, the script will
        merge the distance and contact maps across different copies of the same protein
        pair to generate a single map for each unique protein pair. Defaults to False.

    - **binarize_cmap (bool, optional):**:<br />
        Whether to binarize the contact maps based on the distance cutoff. If True, the
        contact maps will be binarized such that contacts within the distance cutoff are
        assigned a value of 1 and those outside the cutoff are assigned a value of 0.
        Defaults to False.

    - **binarize_dmap (bool, optional):**:<br />
        Whether to binarize the distance maps based on the distance cutoff. If True, the
        distance maps will be binarized such that distances within the cutoff are assigned a
        value of 1 and those outside the cutoff are assigned a value of 0. Defaults to False.

    - **float_dtype (int, optional):**:<br />
        Floating point precision to use for the distance maps in the output files (e.g., 32 or 64).
        Defaults to 64.

    - **int_dtype (int, optional):**:<br />
        Integer precision to use for the contact maps in the output files (e.g., 32 or 64).
        Defaults to 32.

    - **overwrite (bool, optional):**:<br />
        Whether to overwrite existing output files. If False, the script will
        skip processing if output files already exist.
        Defaults to False.

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages.
    """

    command = [
        "python", script_path,
        "--xyzr_file", xyzr_file,
        "--nproc", nproc,
        "--dist_cutoff", dist_cutoff,
        "--frac_cutoff", frac_cutoff,
        "--interaction_map_dir", interaction_map_dir,
        "--plotting", plotting,
        "--float_dtype", float_dtype,
        "--int_dtype", int_dtype,
    ]

    if merge_copies:
        command.append("--merge_copies")

    if binarize_cmap:
        command.append("--binarize_cmap")

    if binarize_dmap:
        command.append("--binarize_dmap")

    if overwrite:
        command.append("--overwrite")

    if logger is not None:
        logger.info("Running interaction_map with command:")
        logger.info(" ".join(map(str, command)))

    os.system(" ".join(map(str, command)))

def interaction_metapatches(
    script_path: str,
    interaction_map_dir: str,
    threshold: int = 40,
    plotting: str = "plotly",
    logger: logging.Logger | None = None,
):
    """ Run the script `interaction_metapatches.py` to identify coarse-level
    interacting regions from the interaction map data.

    ## Arguments:

    - **script_path (str)**:<br />
        Path to the `interaction_metapatches.py` script.

    - **interaction_map_dir (str)**:<br />
        Directory containing the interaction map data generated from the `interaction_map.py` script.
        This directory should contain the distance and contact maps for each frame.

    - **threshold (int, optional):**:<br />
        Gap threshold (in number of residues) for merging neighboring patches into metapatches.
        If the gap between two patches is less than or equal to this threshold, they will be
        merged into a single metapatch. Defaults to 40 residues.

    - **plotting (str, optional):**:<br />
        Plotting library to use for visualizations. Options include "matplotlib" for static plots,
        "plotly" for interactive plots, or "no_plot" to skip plotting. Defaults
        to "plotly".

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages.
    """

    command = [
        "python", script_path,
        "--interaction_map_dir", interaction_map_dir,
        "--threshold", threshold,
        "--plotting", plotting,
    ]

    if logger is not None:
        logger.info("Running interacting_metapatches with command:")
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

    ## Arguments:

    - **script_path (str)**:<br />
        Path to the `fit_to_binding_data.py` script.

    - **xyzr_file (str)**:<br />
        Path to the input hdf5 file containing XYZR data. This file is generated
        from the `rmf_to_xyzr.py` script and contains the coordinates and radii
        of the beads for each molecule in the system.

    - **input_config (str)**:<br />
        Path to the input configuration file containing the binding data and
        parameters for fitting the models to the binding data. This file should
        specify the binding data to be used for fitting, as well as any parameters
        needed for the fitting process.

    - **output_dir (str)**:<br />
        Directory to save the outputs from the fitting process. The outputs may
        include the fitted models, fit scores, and any visualizations generated from
        the fitting process.

    - **nproc (int, optional):**:<br />
        Number of processors to use for the fitting process.

    - **merge_copies (bool, optional):**:<br />
        Whether to merge maps across copies for protein pairs.

    - **float_dtype (int, optional):**:<br />
        Float dtype for calculations (e.g., 32 or 64)

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages.
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

def fit_to_immunoem_data(
    script_path: str,
    xyzr_file: str,
    input_config: str,
    output_dir: str,
    nproc: int = 24,
    merge_copies: bool = True,
    float_dtype: int = 64,
    logger: logging.Logger | None = None,
):
    """ Run the script `fit_to_immunoem_data.py` and fit the models to immuno-EM data.

    ## Arguments:

    - **script_path (str)**:<br />
        Path to the `fit_to_immunoem_data.py` script.

    - **xyzr_file (str)**:<br />
        Path to the input hdf5 file containing XYZR data. This file is generated
        from the `rmf_to_xyzr.py` script and contains the coordinates and radii
        of the beads for each molecule in the system.

    - **input_config (str)**:<br />
        Path to the input configuration file containing the immuno-EM data and
        parameters for fitting the models to the immuno-EM data. This file should
        specify the immuno-EM data to be used for fitting, as well as any parameters
        needed for the fitting process.

    - **output_dir (str)**:<br />
        Directory to save the outputs from the fitting process. The outputs may
        include the fitted models, fit scores, and any visualizations generated from
        the fitting process.

    - **nproc (int, optional):**:<br />
        Number of processors to use for the fitting process.

    - **merge_copies (bool, optional):**:<br />
        Whether to merge maps across copies for protein pairs.

    - **float_dtype (int, optional):**:<br />
        Float dtype for calculations (e.g., 32 or 64)

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages.
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
        logger.info("Running fit_to_immunoem_data with command:")
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
        "interaction_map",
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
        interaction_map
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
            voxel_size=5.0,
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

        correlate_cluster_sample_densities(
            script_path=f"{IMP_TOOLOBX_PATH}/analysis/correlate_cluster_sample_densities.py",
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
            overwrite=False,
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
    # Interaction maps
    ###########################################################################

    if "interaction_map" in args.scripts_to_run:

        interaction_map_dir = os.path.join(ANALYSIS_OUTPUT_PATH, "interaction_map")
        os.makedirs(interaction_map_dir, exist_ok=True)
        assert os.path.exists(xyzr_output_path), (
            f"""XYZR output file {xyzr_output_path} does not exist.
            Please check if rmf_to_xyzr has been run successfully.
            """
        )
        interaction_map(
            script_path=f"{IMP_TOOLOBX_PATH}/analysis/interaction_map.py",
            xyzr_file=xyzr_output_path,
            interaction_map_dir=interaction_map_dir,
            nproc=24,
            dist_cutoff=10.0,
            frac_cutoff=0.25,
            plotting="matplotlib",
            merge_copies=False,
            binarize_cmap=False,
            binarize_dmap=False,
            float_dtype=64,
            int_dtype=32,
            overwrite=False,
            logger=logger,
        )
        lap = time.perf_counter()
        logger.info(f"Completed interaction_map in {lap - start_t:0.4f} seconds")

    else:
        logger.info("Skipping interaction_map as per user request.")

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
    # Fit to immunoem data (optional)
    ###########################################################################

    if "fit_to_immunoem_data" in args.scripts_to_run:

        fit_to_immunoem_data_dir = os.path.join(
            ANALYSIS_OUTPUT_PATH, "fit_to_immunoem_data"
        )
        os.makedirs(fit_to_immunoem_data_dir, exist_ok=True)
        input_config = f"{IMP_TOOLOBX_PATH}/analysis/fit_to_immunoem.json"
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
        fit_to_immunoem_data(
            script_path=f"{IMP_TOOLOBX_PATH}/analysis/fit_to_immunoem_data.py",
            xyzr_file=xyzr_output_path,
            input_config=input_config,
            output_dir=fit_to_immunoem_data_dir,
            nproc=24,
            merge_copies=True,
            float_dtype=64,
            logger=logger
        )
        lap = time.perf_counter()
        logger.info(f"Completed fit_to_immunoem_data in {lap - start_t:0.4f} seconds")

    ###########################################################################
    logger.info("End-to-end analysis completed.")
    logger.info("Ran following scripts:")
    for script in args.scripts_to_run:
        logger.info(f" - {script}")
    end_time = time.perf_counter()
    logger.info(f"Total time: {end_time - start_t:0.4f} seconds")