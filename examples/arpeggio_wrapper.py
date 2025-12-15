from string import Template
import argparse
import yaml
import os
import pandas as pd
import warnings
from IMP_Toolbox.analysis.arpeggio.arpeggio_constants import DOCKER_BASE_COMMAND
from IMP_Toolbox.analysis.arpeggio.arpeggio_analysis import (
    run_arpeggio_docker,
    get_interactions,
    ArpeggioParser,
    ArpeggioChimeraX,
)
from IMP_Toolbox.analysis.arpeggio.arpeggio_constants import (
    CHOSEN_ARI_TYPES,
    CHOSEN_CONTACT_TYPES,
    CHOSEN_CLASHES,
    PI_PI_STACKING_TYPES,
    INTERACTION_DF_COLS
)

if __name__ == "__main__":

    args = argparse.ArgumentParser(
        description="Run Arpeggio Docker container on a given PDB file."
    )
    args.add_argument(
        "-i",
        "--input_config",
        type=str,
        required=False,
        default="input/config.yaml",
        help="Path to the YAML configuration file."
    )
    args.add_argument(
        "-d",
        "--docker_results_dir",
        type=str,
        required=False,
        default=f"output/arpeggio_docker_results",
        help="Path to the output directory to store results."
    )
    args.add_argument(
        "-a",
        "--analysis_dir",
        type=str,
        required=False,
        default=f"output/arpeggio_analysis",
        help="Path to the directory to store analysis results."
    )
    args.add_argument(
        "-p",
        "--processed_structures_dir",
        type=str,
        required=False,
        default=f"output/processed_structures",
        help="Path to the directory containing processed PDB structures."
    )
    args.add_argument(
        "-l",
        "--level",
        type=str,
        choices=["residue", "atom"],
        default="residue",
        help="Level at which to compute interactions: 'residue' or 'atom'."
    )
    args.add_argument(
        "--chimerax_commands",
        action="store_true",
        help="If set, generate ChimeraX command files for visualization."
    )
    args.add_argument(
        "--dry_run",
        action="store_true",
        help="If set, only print the Docker commands without executing them."
    )
    args = args.parse_args()

    docker_base_command = Template(DOCKER_BASE_COMMAND)
    container_path = os.path.abspath(args.processed_structures_dir)
    dry_run = args.dry_run
    chimerax_commands = args.chimerax_commands
    assert isinstance(chimerax_commands, bool)
    assert isinstance(dry_run, bool)

    ###########################################################################
    # Load input configuration
    ###########################################################################
    with open(args.input_config, 'r') as f:
        config = yaml.safe_load(f)

    arpeggio_results = config["arpeggio_results"]
    docker_result_dir = os.path.abspath(args.docker_results_dir)

    for result_head, result_metadata in arpeggio_results.items():

        run_arpeggio_docker(
            docker_base_command,
            container_path,
            docker_result_dir,
            result_head,
            result_metadata,
            dry_run=dry_run,
        )

        # Residue selections that were used for Arpeggio
        selections = []
        for sel in result_metadata.get("selections", []):
            _, chain, res, _ = sel.split("/")
            selections.append((chain, res))
        res_selections = [res for _chain, res in selections]

        contacts_path = os.path.join(
            docker_result_dir, result_head, f"{result_head}.contacts"
        )

        arpeggio_parser = ArpeggioParser(
            result_dir=os.path.join(docker_result_dir, result_head),
        )

        #######################################################################
        # Parse Arpeggio result files
        #######################################################################
        contacts_df = arpeggio_parser.parse_contacts(split_atom_col=True)
        # print(contacts_df)

        ri_df = arpeggio_parser.parse_ri(
            add_marker_id=True,
            split_residue_col=True
        )
        # print(ri_df)

        rings_df = arpeggio_parser.parse_rings(
            add_marker_id=True,
            split_residue_col=True,
        )
        # print(rings_df)

        ari_df = arpeggio_parser.parse_ari(
            split_atom_col=True,
            split_residue_col=True,
        )
        # print(ari_df)

        interactions_df = get_interactions(
            ri_df=ri_df,
            contacts_df=contacts_df,
            ari_df=ari_df,
            contact_level=args.level,
        )

        #######################################################################
        # first residue is always the selected residue
        #######################################################################
        interactions_df_ = pd.DataFrame(columns=interactions_df.columns)

        for idx, row in interactions_df.iterrows():

            if row["res_1"] in res_selections:
                interactions_df_ = pd.concat(
                    [interactions_df_, pd.DataFrame([row])],
                    ignore_index=True
                )

            elif row["res_2"] in res_selections:
                new_row = row.copy()
                new_row["chain_1"] = row["chain_2"]
                new_row["res_1"] = row["res_2"]
                new_row["chain_2"] = row["chain_1"]
                new_row["res_2"] = row["res_1"]
                if "atom_1" in row.index and "atom_2" in row.index:
                    new_row["atom_1"] = row.get("atom_2", None)
                    new_row["atom_2"] = row.get("atom_1", None)
                interactions_df_ = pd.concat(
                    [interactions_df_, pd.DataFrame([new_row])],
                    ignore_index=True
                )

        #######################################################################
        # merge by residue pairs if information is required at residue level
        #######################################################################
        if args.level == "residue":
            row_dict = {}
            for idx, row in interactions_df_.iterrows():

                key = (row["chain_1"], row["res_1"], row["chain_2"], row["res_2"])

                if key not in row_dict:
                    row_dict[key] = set()

                for interaction in row["interaction_type"].split(","):
                    row_dict[key].add(interaction)

            interactions_df = pd.DataFrame(columns=interactions_df_.columns)

            for key, interactions in row_dict.items():
                new_row = {
                    "chain_1": key[0],
                    "res_1": key[1],
                    "chain_2": key[2],
                    "res_2": key[3],
                    "interaction_type": ",".join(sorted(interactions))
                }
                interactions_df = pd.concat(
                    [interactions_df, pd.DataFrame([new_row])],
                    ignore_index=True
                )

        #######################################################################
        # Formatting and saving interactions dataframe
        #######################################################################
        interactions_df = interactions_df.sort_values(
            by=INTERACTION_DF_COLS[args.level][:-1]
        ).reset_index(drop=True)

        total_interactions_type = (
            CHOSEN_CONTACT_TYPES + CHOSEN_CLASHES
            + CHOSEN_ARI_TYPES + ["PI-PI_STACKING"]
        )

        for interaction in total_interactions_type:
            interactions_df[interaction] = (
                interactions_df["interaction_type"].apply(
                    lambda x: 1 if interaction in x.split(",") else 0
                )
            )

        del interactions_df["interaction_type"]

        interactions_df = interactions_df[
            INTERACTION_DF_COLS[args.level][:-1] + sorted(total_interactions_type)
        ]

        outpath = os.path.join(
            args.analysis_dir,
            result_head,
            f"{result_head}_interactions_{args.level}_level.csv"
        )
        os.makedirs(os.path.dirname(outpath), exist_ok=True)
        interactions_df.to_csv(outpath, index=False)

        #######################################################################
        # Generate ChimeraX command files
        #######################################################################
        if chimerax_commands is False:
            continue

        ArpeggioChimeraX.add_rings(
            rings_df=rings_df,
            save_path=outpath.replace(
                f"_interactions_{args.level}_level.csv",
                f"_rings.cxc"
            ),
            add_selections=True,
        )

        ArpeggioChimeraX.add_ri(
            ri_df=ri_df,
            save_path=outpath.replace(
                f"_interactions_{args.level}_level.csv",
                f"_ri.cxc"
            ),
        )

        ArpeggioChimeraX.add_contacts(
            contacts_df=contacts_df,
            save_path=outpath.replace(
                f"_interactions_{args.level}_level.csv",
                f"_contacts.cxc"
            ),
        )

        ArpeggioChimeraX.add_ari(
            ari_df=ari_df,
            save_path=outpath.replace(
                f"_interactions_{args.level}_level.csv",
                f"_ari.cxc"
            ),
        )

        print(f"ChimeraX command files saved in {os.path.abspath(outpath)}.")