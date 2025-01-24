import argparse
import numpy as np
import os
import sys
import matplotlib.pyplot as plt


class WeAreConfident:
    """ Find confident interactions between two proteins and plot the PAE matrix and the confident contacts.
    """

    def __init__(self):
        self.model = 0
        self.confident_contacts = {}
        self.pae_cutoff = 10
        self.plddt_cutoff = 70
        self.contact_threshold = 8
        self.af_predictions_pairwsie = "../inputs/AF_predictions"
        self.output_dir = "../output/af_contact_maps"
        self.af_pipeline_path = None

    def get_protein_names(self, job_name):
        """Get the names of the two proteins from the job name.

        Args:
            job_name (str): The name of the job.

        Returns:
            tuple: The names of the two proteins.
        """
        if len(job_name.split("_")) == 6:
            p1_name = job_name.split("_")[0]
            p2_name = job_name.split("_")[3]

        else:
            print("Job name not in the expected format")
            print("Expected format: P1_copy_StarttoEnd_P2_copy_StarttoEnd")
            p1_name, p2_name = None, None

        return p1_name, p2_name

    def setup_interaction(self, structure_path, pae_path):
        """Setup the interaction object from afpipeline.

        Args:
            structure_path (str): The path to the structure file.
            pae_path (str): The path to the PAE file.

        Returns:
            Interaction: The interaction object.
        """
        sys.path.append(self.af_pipeline_path)
        from main import Interaction
        interaction = Interaction(structure_path, pae_path)
        interaction.pae_cutoff = self.pae_cutoff
        interaction.plddt_cutoff = self.plddt_cutoff
        interaction.contact_threshold = self.contact_threshold

        return interaction


    def find_interacting_residues(self, interaction_mat):
        """Find the interacting residues from the interaction matrix.

        Args:
            interaction_mat (np.array): The interaction matrix.

        Returns:
            tuple: The interacting residues of the two proteins.
        """

        interacting_res_pairs = np.argwhere(interaction_mat == 1)
        interacting_res_pairs = interacting_res_pairs + 1
        interacting_res_pairs = [list(res_pair) for res_pair in interacting_res_pairs]

        p1_contact_res = list(set([res_pair[0] for res_pair in interacting_res_pairs]))
        p2_contact_res = list(set([res_pair[1] for res_pair in interacting_res_pairs]))

        return p1_contact_res, p2_contact_res

    def get_avg_pae_mat(self, interaction):
        """Get the average PAE matrix.

        Args:
            interaction (Interaction): The interaction object.

        Returns:
            np.array: The average PAE matrix.
        """

        prediction_data = interaction.get_data_dict()
        pae_mat = interaction.get_pae(prediction_data)

        return pae_mat

    def calculate_chain_lengths(self, interaction):
        """Calculate the lengths of the two chains.

        Args:
            interaction (Interaction): The interaction object.
        Returns:
            tuple: The lengths of the two chains.
        """

        lengths_dict = interaction.get_chain_lengths(
            interaction.get_residue_positions()
        )

        p1_length = lengths_dict["A"]
        p2_length = lengths_dict["B"]

        return p1_length, p2_length

    def matrix_cuts(self, p1_contact_res, p2_contact_res):
        """Get the four points to cut the matrix.

        Args:
            p1_contact_res (list): The interacting residues of protein 1.
            p2_contact_res (list): The interacting residues of protein 2.

        Returns:
            list: The four points to cut the matrix.
        """

        h1 = np.min(p1_contact_res)
        h2 = np.max(p1_contact_res)
        v1 = np.min(p2_contact_res)
        v2 = np.max(p2_contact_res)
        four_pts = [h1, h2, v1, v2]

        return four_pts

    def process_pae_mat(self, pae_mat, p1_length, interaction_mat, matrix_cuts_):
        """Process the PAE matrix.

        Args:
            pae_mat (np.array): The PAE matrix.
            p1_length (int): The length of the first protein
            interaction_mat (np.array): The interaction matrix.
            matrix_cuts_ (list): The four points to cut the matrix.

        Returns:
            tuple: The sub interaction matrix and the sub PAE matrix.
        """

        h1, h2, v1, v2 = matrix_cuts_
        sub_interaction_mat = interaction_mat[h1 - 1 : h2 + 1, v1 - 1 : v2 + 1]
        sub_pae_mat = pae_mat[h1 - 1 : h2 + 1, v1 - 1 + p1_length : v2 + 1 + p1_length]

        return sub_interaction_mat, sub_pae_mat

    def get_confident_contacts(self):
        """Get the confident interactions between two proteins and plot the PAE matrix and the confident contacts.
        """

        for job_name in os.listdir(self.af_predictions_pairwsie)[0:]:

            if not os.path.isdir(os.path.join(self.af_predictions_pairwsie, job_name)):
                continue
            print("Processing", job_name)

            p_names = self.get_protein_names(job_name)

            if p_names[0] is None or p_names[1] is None:
                continue

            structure_path = os.path.join(
                self.af_predictions_pairwsie,
                job_name,
                f"fold_{job_name}_model_{self.model}.cif",
            )
            pae_path = os.path.join(
                self.af_predictions_pairwsie,
                job_name,
                f"fold_{job_name}_full_data_{self.model}.json",
            )

            interaction = self.setup_interaction(structure_path, pae_path)
            pae_mat = self.get_avg_pae_mat(interaction)

            p1_length, p2_length = self.calculate_chain_lengths(interaction)

            interacting_regions = {
                "A": [1, p1_length],
                "B": [1, p2_length]
            }

            interaction_mat = interaction.get_confident_interactions(
                interacting_regions
            )

            p1_contact_res, p2_contact_res = self.find_interacting_residues(
                interaction_mat
            )

            if len(p1_contact_res) == 0 or len(p2_contact_res) == 0:
                continue

            four_pts = self.matrix_cuts(p1_contact_res, p2_contact_res)
            sub_interaction_mat, sub_pae_mat = self.process_pae_mat(
                pae_mat, p1_length, interaction_mat, four_pts
            )

            self.plot_interactions(
                p_names, pae_mat, p1_length, four_pts, sub_interaction_mat, sub_pae_mat
            )

    def plot_interactions(
        self, p_names, pae_mat, p1_length, four_pts, sub_interaction_mat, sub_pae_mat
    ):
        """Plot the PAE matrix, interaction matrix and the confident contacts.

        Args:
            p_names (tuple): The names of the two proteins.
            pae_mat (np.array): The PAE matrix.
            p1_length (int): The length of the first protein.
            four_pts (list): The four points to cut the matrix.
            sub_interaction_mat (np.array): The interaction matrix.
            sub_pae_mat (np.array): The PAE matrix of the confident contacts.
        """
        p1_name, p2_name = p_names

        h1, h2, v1, v2 = four_pts

        _, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 5))

        ax1.imshow(pae_mat, cmap="hot")

        # make a square on the interacting region
        ax1.plot(
            [v1 - 1 + p1_length, v2 + p1_length],
            [h1 - 1, h1 - 1],
            color="blue",
            lw=1
        )
        ax1.plot(
            [v1 - 1 + p1_length, v2 + p1_length],
            [h2, h2],
            color="blue",
            lw=1
        )
        ax1.plot(
            [v1 - 1 + p1_length, v1 - 1 + p1_length],
            [h1 - 1, h2],
            color="blue",
            lw=1
        )
        ax1.plot(
            [v2 + p1_length, v2 + p1_length],
            [h1 - 1, h2],
            color="blue",
            lw=1
        )

        # add arrow to show the interacting region
        ax1.annotate(
            "",
            xy=(v2 + p1_length, h1 - 1 - 100),
            xytext=(v1 - 1 + p1_length, h1 - 1),
            arrowprops=dict(arrowstyle="<-", lw=1),
        )

        ax1.set_xticks(
            np.arange(0, len(pae_mat), len(pae_mat) - 1),
            [str(x) for x in np.arange(1, len(pae_mat) + 1, len(pae_mat) - 1)],
        )
        ax1.set_yticks(
            np.arange(0, len(pae_mat), len(pae_mat) - 1),
            [str(x) for x in np.arange(1, len(pae_mat) + 1, len(pae_mat) - 1)],
        )

        ax1.set_title(f"PAE Matrix ({p1_name}, {p2_name})")
        ax1.set_xlabel("Residues")
        ax1.set_ylabel("Residues")
        cbar = plt.colorbar(ax1.imshow(pae_mat, cmap="hot", interpolation="nearest"))
        cbar.set_label("PAE")

        ax2.imshow(sub_pae_mat, cmap="hot", interpolation="nearest")
        ax2.set_xticks(
            np.arange(0, v2 - v1 + 1, v2 - v1),
            [str(x) for x in np.arange(v1, v2 + 1, v2 - v1)],
        )
        ax2.set_yticks(
            np.arange(0, h2 - h1 + 1, h2 - h1),
            [str(x) for x in np.arange(h1, h2 + 1, h2 - h1)],
        )
        ax2.set_title(f"Interaction region PAE Matrix ({p1_name}, {p2_name})")
        ax2.set_xlabel("Chain B")
        ax2.set_ylabel("Chain A")

        ax3.imshow(sub_interaction_mat, cmap="hot", interpolation="nearest")
        ax3.set_xticks(
            np.arange(0, v2 - v1 + 1, v2 - v1),
            [str(x) for x in np.arange(v1, v2 + 1, v2 - v1)],
        )
        ax3.set_yticks(
            np.arange(0, h2 - h1 + 1, h2 - h1),
            [str(x) for x in np.arange(h1, h2 + 1, h2 - h1)],
        )
        ax3.set_title(f"Confident Contacts ({p1_name}, {p2_name})")
        os.makedirs(self.output_dir, exist_ok=True)
        plt.savefig(os.path.join(self.output_dir, f"{p1_name}_{p2_name}.png"))
        plt.close()


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument(
        "--pae_cutoff",
        type=int,
        default=10,
        help="The PAE cutoff."
    )
    args.add_argument(
        "--plddt_cutoff",
        type=int,
        default=70
    )
    args.add_argument(
        "--contact_threshold",
        type=int,
        default=8,
        help="The contact threshold (distance in angstroms)."
    )
    args.add_argument(
        "--af_predictions",
        type=str,
        default="../inputs/AF_predictions",
        help="The path to the pairwise AF predictions."
    )
    args.add_argument(
        "--output_dir",
        type=str,
        default="../output/af_contact_maps",
        help="The output directory to save the plots."
    )
    args.add_argument(
        "--af_pipeline",
        type=str,
        help="The path to the af pipeline.",
        required=True
    )

    args = args.parse_args()
    wc = WeAreConfident()

    wc.contact_threshold = args.contact_threshold
    wc.pae_cutoff = args.pae_cutoff
    wc.plddt_cutoff = args.plddt_cutoff
    wc.af_predictions_pairwsie = args.af_predictions
    wc.output_dir = args.output_dir
    wc.af_pipeline_path = args.af_pipeline

    wc.get_confident_contacts()