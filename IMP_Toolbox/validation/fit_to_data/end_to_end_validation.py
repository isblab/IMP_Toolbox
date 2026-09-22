import os
import sys
import time
import random
import argparse
import logging
import getpass
from pathlib import Path
_user = getpass.getuser()

here = os.path.dirname(os.path.abspath(__file__))

logging.basicConfig(
    filename=os.path.join(here, 'end_to_end_analysis.log'),
    filemode="a+",
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO
)

random.seed(47)

def fit_to_binding_data(
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

    script_path = Path(here) / "fit_to_data" / "binding_data.py"

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

def fit_to_em_data(
    input_config: str,
    output_dir: str,
    logger: logging.Logger | None = None,
):
    """ Run the script `fit_to_data/em_data.py` and obtain cross-correlation
    of the localization probability densities witht the experimental maps.

    ## Arguments:

    - **input_config (str)**:<br />
        Path to the input configuration file containing the parameters for fitting
        the models to the EM data. This file should specify the paths to the
        localization probability density maps, and the experimental EM maps.

    - **output_dir (str)**:<br />
        Directory to save the outputs from the fitting process. The outputs may
        include the fit scores, correlation maps, and any visualizations generated from
        the fitting process.

    - **logger (logging.Logger | None, optional):**:<br />
        Logger for logging messages.
    """

    script_path = Path(here) / "fit_to_data" / "em_data.py"

    command = [
        "python", script_path,
        "--input", input_config,
        "--output_dir", output_dir,
    ]

    if logger is not None:
        logger.info("Running fit_to_em_data with command:")
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
            "fit_to_binding_data",
            "fit_to_em_data",
        ],
        help="List of scripts to run in sequence (default: all) \
            (default: [fit_to_binding_data, fit_to_em_data])"
    )
    args = parser.parse_args()

    assert all([s in [
        "fit_to_binding_data",
        "fit_to_em_data",
    ] for s in args.scripts_to_run]), (
        f"""
        Invalid script name in scripts_to_run.
        Valid options are:
        fit_to_binding_data
        fit_to_em_data
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
    # Fit to binding data (optional)
    ###########################################################################

    xyzr_output_path = os.path.join(
        ANALYSIS_OUTPUT_PATH, "sampcon_extracted_frames_xyzr.h5"
    )
    if "fit_to_binding_data" in args.scripts_to_run:

        fit_to_binding_data_dir = os.path.join(
            ANALYSIS_OUTPUT_PATH, "fit_to_binding_data"
        )
        os.makedirs(fit_to_binding_data_dir, exist_ok=True)
        input_config = f"{here}/analysis/fit_to_binding.json"
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