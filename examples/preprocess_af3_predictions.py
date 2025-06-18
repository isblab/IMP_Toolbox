from af_pipeline.RankPredictions import RankAF3Predictions
from argparse import ArgumentParser
import yaml
from utils import update_config
import os
from pprint import pprint

if __name__ == "__main__":

    args = ArgumentParser()

    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="./input/config.yaml",
        help="Path to input yaml file containing the target proteins and their uniprot ids",
    )

    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="./output/ranked_af3_predictions",
        help="Output directory for ranked predictions",
    )

    args.add_argument(
        "-k",
        "--pred_keys",
        nargs="+",
        default=["Lb2Cas12a_RNA_DNA_complex_8I54:1", "Actin_profilin:1,2"],
        help="Keys of AF3 predictions to rank. If not provided, all predictions will be ranked.",
    )

    args = args.parse_args()

    config_yaml = yaml.load(open(args.input), Loader=yaml.FullLoader)

    input_yml = config_yaml.get("af_jobs", None)

    af3_predictions = config_yaml.get("af_predictions", None)
    if af3_predictions is None:
        raise ValueError("No AF3 predictions found in the input yaml file.")

    pred_keys = args.pred_keys

    ranker = RankAF3Predictions(af3_predictions)

    for pred_key in pred_keys:

        af_key, job_ids = pred_key.split(":")
        job_ids = job_ids.split(",")
        for job_id in job_ids:

            ranker.add_key(af_key)
            ranker.add_job_id(int(job_id))
            print(af_key, job_id)
            print(ranker.af_key, ranker.job_id)
            print(ranker.af_pred_key)
            print(ranker.af_target_key)

            pred_dir = af3_predictions[ranker.af_pred_key][ranker.job_id]
            if not os.path.exists(pred_dir):
                print(f"Prediction directory {pred_dir} does not exist. Skipping.")
                continue
            pred_yaml = ranker.extract_af3_best_pred_data(af_jobs=input_yml)
            pprint(pred_yaml)
            print("*"*150)

            if len(pred_yaml) > 0:
                print("Updating config yaml with ranked predictions")
                update_config(
                    input_file=args.input,
                    updates={f"{ranker.af_pred_key}_job_{str(ranker.job_id)}": pred_yaml},
                    mode="replace",
                )
            else:
                print(f"No predictions found for {ranker.af_pred_key} in job {ranker.job_id}. Skipping.")
                print("You should update the config file manually.")