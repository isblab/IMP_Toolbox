"""
Module to parse and analyze Arpeggio result files into DataFrames.

Currently supports parsing and analysis of:
- .contacts (Pairwise contacts between two individual atoms)
- .ari (Atom-aromatic ring interactions)
- .ri (Aromatic ring-aromatic ring interactions)
- .rings (Aromatic rings found in the structure)
"""

import os
import warnings
from string import Template
import pandas as pd
from IMP_Toolbox.utils_imp_toolbox.obj_helpers import get_res_range_from_key
from IMP_Toolbox.analysis.arpeggio.arpeggio_constants import (
    CONTACTS_FIELDS,
    RI_FIELDS,
    RING_FIELDS,
    ARI_FIELDS,
    CHOSEN_ARI_TYPES,
    CHOSEN_CONTACT_TYPES,
    CHOSEN_CLASHES,
    PI_PI_STACKING_TYPES,
    INTERACTION_DF_COLS,
    MARKER_COMMAND,
    PBOND_ATTRIBUTES,
    PBOND_COMMAND,
    MARKER_ATTRIBUTES,
    TRANSPARENCY_COMMAND,
    MODEL_ID,
)

def get_interactions(
    ri_df: pd.DataFrame | None = None,
    contacts_df: pd.DataFrame | None = None,
    ari_df: pd.DataFrame | None = None,
    contact_level: str = "residue",
) -> pd.DataFrame:
    """ Combine interactions from contacts, ri, and ari dataframes.

    Args:

        ri_df (pd.DataFrame | None, optional):
            Ring-ring interactions dataframe.

        contacts_df (pd.DataFrame | None, optional):
            Atomic interactions dataframe.

        ari_df (pd.DataFrame | None, optional):
            Atom-ring interactions dataframe.

        contact_level (str, optional):
            Contact level: "atom" or "residue".

    Returns:

        pd.DataFrame:
            DataFrame containing combined interactions.
    """

    assert contact_level in INTERACTION_DF_COLS, (
        f"contact_level must be one of {list(INTERACTION_DF_COLS.keys())}."
    )

    combined_df = pd.DataFrame(columns=INTERACTION_DF_COLS[contact_level])

    if contacts_df is not None:

        contacts_df = contacts_df.reset_index(drop=True)
        combined_df = pd.concat(
            [
                combined_df,
                contacts_df.filter(INTERACTION_DF_COLS[contact_level])
            ],
            ignore_index=True
        )

        for idx, row in contacts_df.iterrows():
            detected_interactions = [
                interaction
                for interaction in CHOSEN_CONTACT_TYPES + CHOSEN_CLASHES
                if row[interaction] == '1'
            ]
            if len(detected_interactions) == 0:
                continue

            detected_interactions = ",".join(detected_interactions)
            combined_df.iloc[idx, -1] = detected_interactions

        combined_df = combined_df.dropna(
            subset=["interaction_type"]
        ).reset_index(drop=True)

    if ri_df is not None:

        ri_df = ri_df.reset_index(drop=True)
        if contact_level == "residue":
            ri_df["interaction_type"] = "PI-PI_STACKING"

        combined_df = pd.concat(
            [
                combined_df,
                ri_df.filter(INTERACTION_DF_COLS[contact_level])
            ],
            ignore_index=True
        )

    if ari_df is not None:

        ari_df = ari_df.reset_index(drop=True)
        combined_df = pd.concat(
            [
                combined_df,
                ari_df.filter(INTERACTION_DF_COLS[contact_level])
            ],
            ignore_index=True
        )

    combined_df = combined_df.drop_duplicates().reset_index(drop=True)

    return combined_df

def run_arpeggio_docker(
    docker_base_command: Template,
    container_path: str,
    docker_result_dir: str,
    result_head: str,
    result_metadata: dict,
    dry_run: bool=False,
    overwrite: bool=False,
):
    path_to_mount = os.path.join(docker_result_dir, result_head)
    if os.path.exists(path_to_mount) and overwrite is False:
        warnings.warn(
            f"""Arpeggio results already exist for {result_head} at \
            {path_to_mount}. Skipping.
            To overwrite, set overwrite=True."""
        )
        return  # already processed

    processed_struct_path = os.path.join(container_path, f"{result_head}.pdb")
    if not os.path.isfile(processed_struct_path):
        raise Exception(
            f"Processed structure file not found for \
            {result_head}: {processed_struct_path}. Skipping."
        )

    os.makedirs(path_to_mount, exist_ok=True)
    # copy processed structure to mount path
    target_path = os.path.join(path_to_mount, f"{result_head}.pdb")
    os.system(f"cp {processed_struct_path} {target_path}")

    docker_command = docker_base_command.substitute(
        path_to_mount=path_to_mount,
        container_path=container_path,
        input_pdb_path=processed_struct_path,
    )

    # arpeggio_sel = " ".join(result_metadata.get("selections", [])) or ""
    arpeggio_sels = result_metadata.get("selections", [""])

    selection_strs = []
    for arpeggio_sel in arpeggio_sels:

        if len(arpeggio_sel) == 0:
            continue

        _, chain, res, _ = arpeggio_sel.split("/")

        if len(res) == 0:
            selection_strs.append(f"/{chain}/")
        else:
            res_list = get_res_range_from_key(res)
            selection_strs.extend([
                f"/{chain}/{res}/" for res in res_list
            ])

        # docker_command = docker_command.strip() + f" -s {arpeggio_sel}"
    docker_command = docker_command.strip() + f" -s {' '.join(selection_strs)}"

    print(f"Running Docker command for {result_head}:\n{docker_command}\n")
    if dry_run is False:
        os.system(docker_command)
    print("=" * 100 + "\n")


class ArpeggioParser:

    def __init__(self, result_dir: str):
        self.result_dir = result_dir
        self.result_head = os.path.basename(result_dir)
        self.contacts_path = os.path.join(result_dir, f"{self.result_head}.contacts")
        self.rings_path = os.path.join(result_dir, f"{self.result_head}.rings")
        self.ri_path = os.path.join(result_dir, f"{self.result_head}.ri")
        self.ari_path = os.path.join(result_dir, f"{self.result_head}.ari")

    def parse_contacts(
        self,
        split_atom_col: bool=True
    ) -> pd.DataFrame | None:
        """ Parses the .contacts file generated by Arpeggio and returns a dataframe.

        NOTE:
        The atoms in the .contacts file are represented as:
            `chain/res_num/atom`
        This function splits these into three separate columns depending on the
        value of `split_atom_col`.

        Args:
            file_path (str): Path to the .contacts file.
            split_atom_col (bool, optional): Whether to split the atom columns into
                chain, res, and atom.

        Returns:
            pd.DataFrame: DataFrame containing the parsed contact information.
        """

        with open(self.contacts_path, 'r') as f:
            lines = f.readlines()

        if len(lines) == 0:
            return None

        df = pd.DataFrame(columns=CONTACTS_FIELDS)

        for line in lines:
            parts = line.split("\t")
            row_to_add = {k: v for k, v in zip(CONTACTS_FIELDS, parts)}
            df = pd.concat([df, pd.DataFrame([row_to_add])], ignore_index=True)

        if not split_atom_col:
            return df

        df[["chain_1", "res_1", "atom_1"]] = df["atom_1"].str.split("/", expand=True)
        df[["chain_2", "res_2", "atom_2"]] = df["atom_2"].str.split("/", expand=True)

        front_cols = [
            "chain_1", "res_1", "atom_1",
            "chain_2", "res_2", "atom_2"
        ]
        df = df[front_cols + [col for col in df.columns if col not in front_cols]]

        return df

    def parse_ri(
        self,
        split_residue_col: bool=True,
        add_marker_id: bool=True
    ) -> pd.DataFrame | None:

        with open(self.ri_path, 'r') as f:
            lines = f.readlines()

        if len(lines) == 0:
            return None

        if add_marker_id:
            rings_df = self.parse_rings(
                add_marker_id=True,
                split_residue_col=False
            )
            r_to_m = {
                row["ring_id"]: row["marker_id"]
                for _, row in rings_df.iterrows()
            }

        df = pd.DataFrame(columns=RI_FIELDS)

        for line in lines:
            parts = line.split("\t")
            row_to_add = {k: v for k, v in zip(RI_FIELDS, parts)}
            row_to_add["ring1_centroid"] = ",".join(
                [str(x) for x in eval(row_to_add["ring1_centroid"])]
            )
            row_to_add["ring2_centroid"] = ",".join(
                [str(x) for x in eval(row_to_add["ring2_centroid"])]
            )
            if add_marker_id:
                row_to_add["marker_id_1"] = r_to_m.get(row_to_add["ring1_id"], "")
                row_to_add["marker_id_2"] = r_to_m.get(row_to_add["ring2_id"], "")

            df = pd.concat([df, pd.DataFrame([row_to_add])], ignore_index=True)

        if not split_residue_col:
            return df

        df[["chain_1", "res_1", "_"]] = df["ring1_residue"].str.split("/", expand=True)
        df[["chain_2", "res_2", "_"]] = df["ring2_residue"].str.split("/", expand=True)

        del df["ring1_residue"]
        del df["ring2_residue"]
        del df["_"]

        front_cols = [
            "chain_1", "res_1",
            "chain_2", "res_2"
        ]
        df = df[front_cols + [col for col in df.columns if col not in front_cols]]

        return df

    def parse_rings(
        self,
        add_marker_id: bool=True,
        split_residue_col: bool=True
    ) -> pd.DataFrame | None:
        """ Parses the .rings file generated by Arpeggio and returns a dataframe.

        Args:
            file_path (str): Path to the .rings file.

        Returns:
            pd.DataFrame: DataFrame containing the parsed ring information.
        """

        with open(self.rings_path, 'r') as f:
            lines = f.readlines()

        if len(lines) == 0:
            return None

        df = pd.DataFrame(columns=RING_FIELDS)

        for idx, line in enumerate(lines):
            parts = line.split("\t")
            row_to_add = {k: v for k, v in zip(RING_FIELDS, parts)}
            row_to_add["ring_centroid"] = ",".join(
                [str(x) for x in eval(row_to_add["ring_centroid"])]
            )

            if add_marker_id:
                row_to_add["marker_id"] = f"{idx+1}"

            df = pd.concat([df, pd.DataFrame([row_to_add])], ignore_index=True)

        if not split_residue_col:
            return df

        df[["chain", "res", "_"]] = df["ring_residue"].str.split("/", expand=True)

        del df["ring_residue"]
        del df["_"]

        front_cols = ["chain", "res", "marker_id"]
        df = df[front_cols + [col for col in df.columns if col not in front_cols]]

        return df

    def parse_ari(
        self,
        split_atom_col: bool=True,
        split_residue_col: bool=True,
        add_marker_id: bool=True,
    ) -> pd.DataFrame | None:
        """ Parses the .ari file generated by Arpeggio and returns a dataframe.

        Args:
            file_path (str): Path to the .ari file.
            split_atom_col (bool, optional): Whether to split the atom columns into
                chain, res, and atom.

        Returns:
            pd.DataFrame | None: DataFrame containing the parsed contact information.
        """

        with open(self.ari_path, 'r') as f:
            lines = f.readlines()

        if len(lines) == 0:
            return None

        if add_marker_id:
            rings_df = self.parse_rings(
                add_marker_id=True,
                split_residue_col=False
            )
            r_to_m = {
                row["ring_id"]: row["marker_id"]
                for _, row in rings_df.iterrows()
            }

        df = pd.DataFrame(columns=ARI_FIELDS)

        for _, line in enumerate(lines):
            parts = line.split("\t")
            row_to_add = {k: v for k, v in zip(ARI_FIELDS, parts)}
            row_to_add["ring_centroid"] = ",".join(
                [str(x) for x in eval(row_to_add["ring_centroid"])]
            )

            row_to_add["interaction_type"] = ",".join(
                [itype.strip() for itype in eval(row_to_add["interaction_type"])]
            )

            if add_marker_id:
                row_to_add["marker_id"] = r_to_m.get(row_to_add["ring_id"], "")

            df = pd.concat([df, pd.DataFrame([row_to_add])], ignore_index=True)

        if not split_atom_col and not split_residue_col:
            return df

        df[["chain_2", "res_2", "_"]] = df["residue"].str.split("/", expand=True)
        del df["residue"]
        del df["_"]
        front_cols = ["chain_2", "res_2", "ring_id"]
        df = df[front_cols + [col for col in df.columns if col not in front_cols]]

        df[["chain_1", "res_1", "atom_1"]] = df["atom"].str.split("/", expand=True)
        del df["atom"]
        front_cols = ["chain_1", "res_1", "atom_1"]
        df = df[front_cols + [col for col in df.columns if col not in front_cols]]

        return df


class ArpeggioChimeraX:

    def __init__(self):
        pass

    @staticmethod
    def save_commands_to_file(
        save_path: str,
        commands: list[str]
    ) -> None:
        """ Save a list of ChimeraX commands to a file.

        Args:

            save_path (str):
                Path to save the commands.

            commands (list[str]):
                List of ChimeraX commands to save.
        """

        if not os.path.exists(os.path.dirname(save_path)):
            warnings.warn(
                f"Directory does not exist: {os.path.dirname(save_path)}. \
                Creating a new directory."
            )
            os.makedirs(os.path.dirname(save_path), exist_ok=True)

        with open(save_path, 'w') as f:
            for command in commands:
                f.write(command + "\n")

    @staticmethod
    def add_contacts(
        contacts_df: pd.DataFrame | None = None,
        save_path: str | None = None,
        add_selections: bool = True,
    ):
        """ Generate ChimeraX commands to add pseudobonds for atomic contacts.

        Args:

            contacts_df (pd.DataFrame):
                DataFrame containing inter-atomic contact information.

            save_path (str | None, optional):
                Path to save the commands. If None, returns the commands as a list.

            add_selections (bool, optional):
                Whether to add selection commands for the interacting residues.

        Returns:

            list[str] | None:
                List of ChimeraX commands if save_path is None, else None.
        """

        if contacts_df is None or contacts_df.empty:
            warnings.warn("No contacts found.")
            return []

        command_template = Template(PBOND_COMMAND)
        transparency_template = Template(TRANSPARENCY_COMMAND)
        # model_id = 1
        commands = []
        sub_model_ids = {}
        sub_model_id = 1

        for _, row in contacts_df.iterrows():

            atom1_spec = f"{MODEL_ID}/{row['chain_1']}:{row['res_1']}@{row['atom_1']}"
            atom2_spec = f"{MODEL_ID}/{row['chain_2']}:{row['res_2']}@{row['atom_2']}"

            if add_selections:
                sel_cmd = f"sel add #{atom1_spec.split("@")[0]}#{atom2_spec.split('@')[0]}"
                if sel_cmd not in commands:
                    commands.append(sel_cmd)

            detected_interactions = [
                interaction for interaction in CHOSEN_CONTACT_TYPES + CHOSEN_CLASHES
                if row[interaction] == '1'
            ]

            for interaction in detected_interactions:

                if interaction not in sub_model_ids:
                    sub_model_ids[interaction] = sub_model_id
                    sub_model_id += 1

                command = command_template.substitute(
                    spec1=atom1_spec,
                    spec2=atom2_spec,
                    bond_color=PBOND_ATTRIBUTES.get(
                        interaction, PBOND_ATTRIBUTES["default"]
                        ).get("color", "yellow"),
                    bond_radius=PBOND_ATTRIBUTES.get(
                        interaction, PBOND_ATTRIBUTES["default"]
                        ).get("radius", 0.1),
                    bond_name=interaction,
                    bond_dashes=PBOND_ATTRIBUTES.get(
                        interaction, PBOND_ATTRIBUTES["default"]
                        ).get("dashes", 6),
                )
                commands.append(command)

                command = transparency_template.substitute(
                    model_id=f"{MODEL_ID}.{sub_model_ids.get(interaction, sub_model_id)}",
                    target_spec="p",
                    transparency=PBOND_ATTRIBUTES.get(
                        interaction, PBOND_ATTRIBUTES["default"]
                        ).get("transparency", 0),
                )
                commands.append(command)

        if save_path:
            ArpeggioChimeraX.save_commands_to_file(save_path, commands)

        return commands

    @staticmethod
    def add_rings(
        rings_df: pd.DataFrame | None = None,
        save_path: str | None = None,
        add_selections: bool = True,
    ) -> list[str] | None:
        """ Generate ChimeraX commands to add markers at ring centroids.

        Args:

            rings_df (pd.DataFrame):
                DataFrame containing ring information.

            save_path (str | None, optional):
                Path to save the commands. If None, returns the commands as a list.

            add_selections (bool, optional):
                Whether to add selection commands aromatic residues.

        Returns:

            list[str] | None:
                List of ChimeraX commands if save_path is None, else None.
        """

        if rings_df is None or rings_df.empty:
            warnings.warn("No rings found.")
            return []

        command_template = Template(MARKER_COMMAND)
        commands = []

        for _, row in rings_df.iterrows():
            x, y, z = row["ring_centroid"].split(",")
            command = command_template.substitute(
                marker_model_id=MARKER_ATTRIBUTES["model_id"],
                x=x,
                y=y,
                z=z,
                marker_color=MARKER_ATTRIBUTES["color"],
                marker_radius=MARKER_ATTRIBUTES["radius"],
            )
            commands.append(command)

            if add_selections:
                chain_id = row["chain"]
                res_id = row["res"]
                sel_cmd = f"sel add #{MODEL_ID}/{chain_id}:{res_id}"
                if sel_cmd not in commands:
                    commands.append(sel_cmd)

        commands.append(
            f"rename #{MARKER_ATTRIBUTES['model_id']} ring_centroids"
        )

        if save_path:
            ArpeggioChimeraX.save_commands_to_file(save_path, commands)

        return commands

    @staticmethod
    def add_ri(
        ri_df: pd.DataFrame | None = None,
        save_path: str | None = None,
        add_selections: bool = True,
    ) -> list[str] | None:
        """ Generate ChimeraX commands to add pseudobonds for ring interactions.

        Args:

            ri_df (pd.DataFrame):
                DataFrame containing ring interaction information.

            save_path (str | None, optional):
                Path to save the commands. If None, returns the commands as a list.

            add_selections (bool, optional):
                Whether to add selection commands for the interacting residues.

        Returns:

            list[str] | None:
                List of ChimeraX commands if save_path is None, else None.
        """

        if ri_df is None or ri_df.empty:
            warnings.warn("No ring interactions found.")
            return []

        command_template = Template(PBOND_COMMAND)
        transparency_template = Template(TRANSPARENCY_COMMAND)
        commands = []
        sub_model_ids = {} # each interaction type gets its own idx
        sub_model_id = 1

        for _, row in ri_df.iterrows():

            marker1_id = row["marker_id_1"]
            marker2_id = row["marker_id_2"]
            marker1_spec = f"{MARKER_ATTRIBUTES["model_id"]}/M:{marker1_id}@M"
            marker2_spec = f"{MARKER_ATTRIBUTES["model_id"]}/M:{marker2_id}@M"
            bond_name = row["interaction_type"]

            if bond_name not in sub_model_ids:
                sub_model_ids[bond_name] = sub_model_id
                sub_model_id += 1

            command = command_template.substitute(
                spec1=marker1_spec,
                spec2=marker2_spec,
                bond_color=PBOND_ATTRIBUTES.get(
                    bond_name, PBOND_ATTRIBUTES["ri"]
                )["color"],
                bond_radius=PBOND_ATTRIBUTES.get(
                    bond_name, PBOND_ATTRIBUTES["ri"]
                )["radius"],
                bond_name=bond_name,
                bond_dashes=PBOND_ATTRIBUTES.get(
                    bond_name, PBOND_ATTRIBUTES["ri"]
                )["dashes"],
            )
            commands.append(command)

            if add_selections:
                chain1 = row["chain_1"]
                res1 = row["res_1"]
                chain2 = row["chain_2"]
                res2 = row["res_2"]
                sel_cmd = f"sel add #{MODEL_ID}/{chain1}:{res1}#{MODEL_ID}/{chain2}:{res2}"
                if sel_cmd not in commands:
                    commands.append(sel_cmd)

            command = transparency_template.substitute(
                model_id=f"{MARKER_ATTRIBUTES["model_id"]}.{sub_model_ids[bond_name]}",
                target_spec="p",
                transparency=PBOND_ATTRIBUTES.get(
                    bond_name, PBOND_ATTRIBUTES["ri"]
                    )["transparency"],
            )
            commands.append(command)

        if save_path:
            ArpeggioChimeraX.save_commands_to_file(save_path, commands)

        return commands

    @staticmethod
    def add_ari(
        ari_df: pd.DataFrame | None = None,
        save_path: str | None = None,
        add_selections: bool = True,
    ):
        """ Generate ChimeraX commands to add pseudobonds for atom-ring interactions.

        Args:

            ari_df (pd.DataFrame):
                DataFrame containing atom-ring interaction information.

            save_path (str | None, optional):
                Path to save the commands. If None, returns the commands as a list.

            add_selections (bool, optional):
                Whether to add selection commands for the interacting residues.

        Returns:

            list[str] | None:
                List of ChimeraX commands if save_path is None, else None.
        """

        if ari_df is None or ari_df.empty:
            warnings.warn("No atom-ring interactions found.")
            return []

        command_template = Template(PBOND_COMMAND)
        transparency_template = Template(TRANSPARENCY_COMMAND)
        model_id = 1
        commands = []
        sub_model_ids = {}
        sub_model_id = 1

        for _, row in ari_df.iterrows():
            atom_spec = f"{model_id}/{row['chain_1']}:{row['res_1']}@{row['atom_1']}"
            ring_marker_id = row["marker_id"]
            ring_spec = f"{MARKER_ATTRIBUTES["model_id"]}/M:{ring_marker_id}@M"
            interaction = row["interaction_type"]
            assert interaction in CHOSEN_ARI_TYPES, (
                f"Unexpected interaction type '{interaction}' found in ARI data."
            )
            if interaction not in sub_model_ids:
                sub_model_ids[interaction] = sub_model_id
                sub_model_id += 1

            command = command_template.substitute(
                spec1=atom_spec,
                spec2=ring_spec,
                bond_name=interaction,
                bond_color=PBOND_ATTRIBUTES.get(
                    interaction, PBOND_ATTRIBUTES["ari"]
                    ).get("color", "yellow"),
                bond_radius=PBOND_ATTRIBUTES.get(
                    interaction, PBOND_ATTRIBUTES["ari"]
                    ).get("radius", 0.1),
                bond_dashes=PBOND_ATTRIBUTES.get(
                    interaction, PBOND_ATTRIBUTES["ari"]
                    ).get("dashes", 6),
            )
            commands.append(command)

            if add_selections:
                chain1 = row["chain_1"]
                res1 = row["res_1"]
                chain2 = row["chain_2"]
                res2 = row["res_2"]
                sel_cmd = f"sel add #{MODEL_ID}/{chain1}:{res1}#{MODEL_ID}/{chain2}:{res2}"
                if sel_cmd not in commands:
                    commands.append(sel_cmd)

            command = transparency_template.substitute(
                model_id=f"{model_id}.{sub_model_ids[interaction]}",
                target_spec="p",
                transparency=PBOND_ATTRIBUTES.get(
                    interaction, PBOND_ATTRIBUTES["ari"]
                    ).get("transparency", 0),
            )
            commands.append(command)

        if save_path:
            ArpeggioChimeraX.save_commands_to_file(save_path, commands)

        return commands