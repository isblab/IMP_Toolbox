from collections import defaultdict
import warnings
import os
from utils import read_json

class RankAF3Predictions:

    def __init__(
        self,
        predictions: dict,
    ):
        self.predictions = predictions
        self.af_key = None
        self.af_pred_key = None
        self.af_target_key = None
        self.job_id = None
        self.try_af_offset_from_path = False

    @staticmethod
    def get_data_path_from_structure_path(structure_path: str) -> str:
        """ Get the data path from the structure path

        Specific to AF3 predictions, where the structure file is in .cif format \n
        and the data file is in .json format with a specific naming convention.

        example structure path
            "/path/to/AF3_pred/p1_1_1to100_p2_2_101to200/model_0001.cif"

        The output will be a path to the data file in the format
            "/path/to/AF3_pred/p1_1_1to100_p2_2_101to200/full_data_0001.json"

        Args:
            structure_path (str): Path to the structure file

        Returns:
            str: Path to the data file
        """

        return structure_path.replace(".cif", ".json").replace("model_", "full_data_")


    def add_key(self, af_key: str):
        """ Add key to the AF3 predictions

        Args:
            af_key (str): Key for the AF3 predictions
        """

        self.af_key = af_key
        self.af_pred_key = f"{self.af_key}_af3_predictions"
        self.af_target_key = f"{self.af_key}_targets"


    def add_job_id(self, job_id: str):
        """ Add job id to the AF3 predictions

        Args:
            job_id (str): Job id for the AF3 predictions
        """

        self.job_id = job_id


    @staticmethod
    def get_af3_model_metrics(af3_summary_confidence_file: str) -> dict:
        """Get AF3 model metrics from summary confidence file
        Reads the summary confidence file and extracts the required metrics.
        The summary confidence file is expected to be in JSON format with the
        following keys:

            - fraction_disordered
            - has_clash
            - iptm
            - num_recycles
            - ptm
            - ranking_score

        Args:
            af3_summary_confidence_file (str): Path to AF3 summary confidence file

        Returns:
            dict: Dictionary containing AF3 model metrics
        """

        metric_data = read_json(af3_summary_confidence_file)

        required_data = {
            "fraction_disordered": metric_data.get("fraction_disordered"),
            "has_clash": metric_data.get("has_clash"),
            "iptm": metric_data.get("iptm"),
            "num_recycles": metric_data.get("num_recycles"),
            "ptm": metric_data.get("ptm"),
            "ranking_score": metric_data.get("ranking_score"),
        }

        return required_data


    @staticmethod
    def get_seed_from_job_request(af3_job_request_file: str) -> int:
        """Get seed from job request file

        Args:
            af3_job_request_file (str): Path to job request file

        Returns:
            int: Seed of the model
        """

        job_request_data = read_json(af3_job_request_file)

        return job_request_data[0].get("modelSeeds")[0]


    def get_af3_model_metrics_per_seed(self) -> dict:
        """ Get AF3 model metrics per seed

        Reads the AF3 predictions directory and extracts the model metrics for
        each seed.
        The directory structure is expected to be as follows:

            job_cycle_pred_dir/
              └── job1_pred_dir/
                    └── seed1/
                        ├── summary_confidences_0.json
                        ├── summary_confidences_1.json
                        ├── summary_confidences_2.json
                        ├── summary_confidences_3.json
                        ├── summary_confidences_4.json
                        ├── job_request.json
                        ├── model_0001.cif
                        ├── model_0002.cif
                        ├── model_0003.cif
                        ├── model_0004.cif
                        └── model_0005.cif

        Each seed directory will contain 5 summary confidence files, one for each model.

        The function will return a dictionary with the seed as the key and a list of
        model metrics as the value. Each model metric will be a dictionary containing
        the following
        keys:

            - ranking_score
            - iptm
            - ptm
            - fraction_disordered
            - model_path
            - model_idx

        Args:
            af_pred_key (str): Key for the predictions in the input dictionary

        Returns:
            dict: Dictionary containing the best AF3 model metrics per seed
        """

        job_cycle = self.predictions.get(self.af_pred_key, None)

        if job_cycle is None:
            warnings.warn(
                f"No predictions found for key {self.af_pred_key}. "
                "Skipping ranking of AF3 predictions."
            )
            return None

        pred_dir = (
            self.predictions[self.af_pred_key].get(self.job_id, None)
        )

        if not os.path.exists(pred_dir):
            warnings.warn(
                f"Prediction directory {pred_dir} does not exist. "
                "Skipping ranking of AF3 predictions."
            )
            return None

        seed_prediction_dirs = [
            os.path.join(pred_dir, pred)
            for pred in os.listdir(pred_dir)
            if os.path.isdir(os.path.join(pred_dir, pred))
        ]

        if len(seed_prediction_dirs) == 0:
            warnings.warn(
                f"No seed prediction directories found in {pred_dir}. "
                "Skipping ranking of AF3 predictions."
            )
            return None

        af3_metrics_per_seed = defaultdict(list)

        for af3_seed_prediction_dir in seed_prediction_dirs:

            model_metrics = {}

            summary_confidences_files = [
                os.path.join(af3_seed_prediction_dir, af3_file)
                for af3_file in os.listdir(af3_seed_prediction_dir)
                if "summary_confidences" in af3_file
            ]

            job_request_file = [
                os.path.join(af3_seed_prediction_dir, af3_file)
                for af3_file in os.listdir(af3_seed_prediction_dir)
                if "job_request" in af3_file
            ]

            assert len(summary_confidences_files) == 5, \
            f"There should be 5 summary_confidences files. \
                Found: {len(summary_confidences_files)}"

            assert len(job_request_file) == 1
            "There should be 1 job_request file"

            for summary_confidence_path in summary_confidences_files:

                model_idx = summary_confidence_path.split("_")[-1].split(".")[0]
                af3_model_metrics = self.get_af3_model_metrics(
                    af3_summary_confidence_file=summary_confidence_path
                )
                model_metrics[int(model_idx)] = af3_model_metrics

            # choosing the best model based on ranking score
            model_seed = self.get_seed_from_job_request(job_request_file[0])

            for i, model_metric in model_metrics.items():

                af3_metrics_per_seed[model_seed].append(
                    {
                        "ranking_score": model_metric["ranking_score"],
                        "iptm": model_metric["iptm"],
                        "ptm": model_metric["ptm"],
                        "fraction_disordered": model_metric["fraction_disordered"],
                        "model_path": [
                            os.path.join(af3_seed_prediction_dir, af3_file)
                            for af3_file in os.listdir(af3_seed_prediction_dir)
                            if f"model_{i}" in af3_file and af3_file.endswith("cif")
                        ][0],
                        "model_idx": i,
                    }
                )

        return af3_metrics_per_seed


    def rank_seeds(self) -> tuple:
        """ Rank the AF3 predictions based on the ranking score

        Args:
            pred_key (str): Key for the predictions in the input dictionary

        Returns:
            tuple: Model seed, ranking score, model path, model index
            corresponding to the best model
        """

        af3_metrics_per_seed = self.get_af3_model_metrics_per_seed()

        if af3_metrics_per_seed is None:
            return None, None, None, None

        ranking = []
        for seed, metrics_list in af3_metrics_per_seed.items():
            for metrics in metrics_list:
                ranking.append(
                    (
                        seed,
                        metrics["ranking_score"],
                        metrics["iptm"],
                        metrics["ptm"],
                        metrics["fraction_disordered"],
                        metrics["model_path"],
                        metrics["model_idx"],
                    )
                )

        # Sort by ranking_score
        ranking.sort(
            key=lambda item: (
                item[1],  # ranking_score
                item[6]*(-1), # model index
                item[2], # iptm
                item[3],  # ptm
                item[4],  # fraction_disordered
                item[0],  # seed
            ),
            reverse=True,
        )

        best_model = ranking[0]
        best_model_seed = best_model[0]
        best_ranking_score = best_model[1]
        best_model_path = best_model[5]
        best_model_idx = best_model[6]

        for idx, model in enumerate(ranking):
            print(
                f"Seed: {model[0]}, "
                f"Ranking Score: {model[1]:.2f}, "
                f"iptm: {model[2]}, "
                f"ptm: {model[3]}, "
                f"Fraction Disordered: {model[4]}, "
                f"Model Index: {model[6]}"
            )

        return best_model_seed, best_ranking_score, best_model_path, best_model_idx


    @staticmethod
    def extract_af_offset_from_path(structure_path: str) -> dict:
        """ Extract the offset for AF3 prediction from the structure path

        example structure path
            "/path/to/AF3_pred/p1_1_1to100_p2_2_101to200_1234/model_0001.cif" or
            "/path/to/AF3_pred/p1_1_1to100_p2_2_101to200/model_0001.cif"

        The directory name is expected to be in the format
            p1_copy1_1to100_p2_copy2_101to200_seed

        where p1, p2 are the chain ids, copy1, copy2 are the number of copies,
        and 1to100, 101to200 are the residue ranges.

        The output will be a dictionary with the chain ids as keys and the residue
        ranges as values.

        output:

            {
                "A": [1, 100],
                "B": [101, 200],
                "C": [101, 200]
            }

        Args:
            structure_path (str): Path to the structure file

        Returns:
            af_offset (dict): Dictionary containing the offset for each chain in the
            format
            {
                "A": [start, end],
                "B": [start, end]
            }
        """

        dirname = os.path.basename(os.path.dirname(structure_path))
        af_offset = {}

        assert len(dirname.split("_")) >= 3
        "Invalid directory name for AF3 prediction"

        assert len(dirname.split("_")) % 3 < 27
        "Invalid directory name for AF3 prediction, too many chains"

        if len(dirname.split("_")[0:-1]) % 3 == 0:
            seed_is_present = True

        elif len(dirname.split("_")[0:-1]) % 3 == 2:
            seed_is_present = False

        else:
            raise ValueError(
                "Invalid directory name for AF3 prediction, "
                "should be in the format p1_copy1_1to100_p2_copy2_101to200_seed "
                "or p1_copy1_1to100_p2_copy2_101to200"
            )

        if seed_is_present:
            p_c_r = dirname.split("_")[0:-1]  # Exclude the seed

        else:
            p_c_r = dirname.split("_") # Seed is not present in the path

        atoz = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        chain_count = 0

        for i in range(len(p_c_r) // 3):

            copy_number = p_c_r[i * 3 + 1]
            res_range = p_c_r[i * 3 + 2]

            for _ in range(int(copy_number)):

                chain_id = atoz[chain_count]
                start, end = res_range.split("to")
                af_offset[chain_id] = [int(start), int(end)]
                chain_count += 1

        return af_offset


    def extract_af_offset_from_af_jobs(self, af_jobs: dict| None = None) -> dict:
        """ Extract the offsets for AF3 predictions from the AF jobs
        dictionary in input yaml file

        Args:
            af_jobs (dict): Dictionary containing AF input data

        Returns:
            af_offset (dict): Dictionary containing the offsets for each chain in
            the format
            {
                "A": [start, end],
                "B": [start, end]
            }
        """

        af_offset = {}

        if af_jobs is None:
            warnings.warn("No AF jobs found. Skipping extraction of offsets.")
            return af_offset

        atoz = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        chain_count = 0

        af_jobs = af_jobs.get(self.af_target_key, None)
        af_entities = af_jobs[self.job_id-1].get("entities", None)

        if af_entities is None:
            raise ValueError(
                f"No entities found for job {self.job_id} in AF jobs. "
                "Please check the input yaml file."
            )

        for entity in af_entities:

            for _entity_count in range(entity.get("count", 1)):

                res_range = entity.get("range", None)

                if res_range:
                    af_offset[atoz[chain_count]] = [
                        entity["range"][0],
                        entity["range"][1],
                    ]

                else:
                    warnings.warn(f"Entity {entity['name']} does not have a range.")
                    af_offset[atoz[chain_count]] = []
                    self.try_af_offset_from_path = True

                chain_count += 1

        if len(af_offset) == 0:
            self.try_af_offset_from_path = True

        return af_offset


    def extract_af3_best_pred_data(self, af_jobs: dict| None = None) -> list:
        """ Extract AF3 model metrics and paths to the best model and data for a
        given prediction directory

        Args:
            af_jobs (dict): Dictionary containing AF input data
            af_predictions (dict): Dictionary containing paths to AF3 predictions
            pred_key (str): Type of prediction to find best model for. This is the key
                in the af_predictions dictionary
            pred_yaml (list): List of dictionaries containing paths to the best model
                and data for a given prediction directory
            prediction_dir (str): Path to the prediction directory

        Returns:
            pred_yaml (list): List of dictionaries containing paths to the best model and
                data for a given prediction directory along with offsets
        """

        pred_yaml = []

        _best_seed, _ranking_score, structure_path, _best_model_idx = self.rank_seeds()

        if structure_path is None:
            warnings.warn(
                f"No valid AF3 predictions found for {self.af_pred_key} job {self.job_id}. "
                "Skipping ranking of AF3 predictions."
            )
            return []

        data_path = self.get_data_path_from_structure_path(structure_path)

        af_offset = self.extract_af_offset_from_af_jobs(af_jobs)

        if self.try_af_offset_from_path:

            try:
                af_offset = self.extract_af_offset_from_path(structure_path)

            except ValueError as e:
                print(f"Error extracting AF offset from path: {e}")

        pred_yaml.append(
            {
                "structure_path": structure_path,
                "data_path": data_path,
                "af_offset": af_offset
            }
        )

        return pred_yaml