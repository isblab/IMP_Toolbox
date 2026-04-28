import os
from pathlib import Path
import time
import tqdm
import getpass
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IMP_Toolbox.analysis.rmf_to_xyzr import XYZRParser
from IMP_Toolbox.analysis.interaction.coarse_grained.interaction_map import (
    PairwiseMaps,
    MoleculePairs,
)
from IMP_Toolbox.constants.analysis_constants import (
    F_DTYPES,
    I_DTYPES,
)

_user = getpass.getuser()

class BindingData:

    def __init__(
        self,
        input: str,
        xyzr_parser: XYZRParser,
        merge_copies: bool=False,
        self_interaction: str="allow_copies",
        f_dtype: np.dtype=np.float64,
        i_dtype: np.dtype=np.int32,
    ):

        if not xyzr_parser.is_setup:
            xyzr_parser.parse_xyzr_h5_file()

        self.xyzr_parser = xyzr_parser
        self.merge_copies = merge_copies
        self.f_dtype = f_dtype
        self.i_dtype = i_dtype
        self.num_frames = self.xyzr_parser.xyzr_mat.shape[1]

        sel_mol_pair_handler = MoleculePairs(
            xyzr_parser=xyzr_parser,
            input=input,
            self_interaction=self_interaction,
        )
        mol_pair_handler = MoleculePairs(
            xyzr_parser=xyzr_parser,
            input=None,
            self_interaction=self_interaction,
        )
        sel_mol_pair_handler.process_molecule_pairs()
        self.sel_mol_pairs = sel_mol_pair_handler.mol_pairs
        mol_pair_handler.process_molecule_pairs(
            filter_by=self.sel_mol_pairs,
        )
        self.mol_pairs = mol_pair_handler.mol_pairs
        self.pairwise_maps_handler = PairwiseMaps(
            mol_pairs=self.sel_mol_pairs,
            xyzr_parser=xyzr_parser,
            cutoff=10,
            self_interaction=self_interaction,
            f_dtype=f_dtype,
            i_dtype=i_dtype,
        )

    def compute_pairwise_distance(
        self,
        binding_data_dir: str,
        nproc: int,
        operation: str="min",
        subtract_radii: bool=True,
        overwrite: bool=False,
    ):

        pairwise_distances = self.pairwise_maps_handler.fetch_pairwise_distances(
            output_dir=binding_data_dir,
            nproc=nproc,
            operation=operation,
            subtract_radii=subtract_radii,
            overwrite=overwrite,
        )

        if self.merge_copies:
            self.xyzr_parser.merge_residue_selection_by_copies()
            self.molwise_residues = self.xyzr_parser.molwise_residues

            pairwise_distances = self.pairwise_maps_handler.merge_maps_by_copies(
                pairwise_maps=pairwise_distances,
                map_type="dmap",
                operation="min",
                binarize_map=False,
            )

        self.pairwise_distances = dict(sorted(pairwise_distances.items(), key=lambda x: x[0]))

    def plot_pairwise_distances(
        self,
        output_dir: str,
        outname: str="fit_to_binding_data",
    ):
        """ Plot pairwise distance maps.

        ## Arguments:

        - **output_dir (str)**:<br />
            Directory to save output plots.

        - **pairwise_distances (dict)**:<br />
            Dictionary where keys are pair names formatted as "Molecule1_Molecule2" and
            values are arrays of shape (num_frames,) containing the minimum distance
            across all bead pairs for each frame.
        """

        assert hasattr(self, "pairwise_distances"), (
            "Pairwise distances not computed yet. Please run compute_pairwise_distance() first."
        )

        all_data = list(self.pairwise_distances.values())
        label_names = list(self.pairwise_distances.keys())
        print("Pairwise distance data calculated for all pairs. Now plotting...\n")

        fig, ax = plt.subplots(figsize=(2*len(label_names), 5))
        violinplot = ax.violinplot(
            all_data,
            showmeans=False,
            showextrema=False,
            showmedians=False,
            # quantiles=[[0.9]]*len(all_data),
        )
        for pc in violinplot['bodies']:
            pc.set_facecolor("#B0ABFF")
            pc.set_edgecolor('black')
            pc.set_alpha(1)
        quartile1, medians, quartile3 = np.percentile(
            all_data, [25, 50, 75], axis=1
        )
        def adjacent_values(vals, q1, q3):
            upper_adjacent_value = q3 + (q3 - q1) * 1.5
            upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

            lower_adjacent_value = q1 - (q3 - q1) * 1.5
            lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
            return lower_adjacent_value, upper_adjacent_value
        whiskers = np.array([
            adjacent_values(sorted_array, q1, q3)
            for sorted_array, q1, q3 in zip(all_data, quartile1, quartile3)])
        whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

        inds = np.arange(1, len(medians) + 1)
        ax.scatter(inds, medians, marker='o', color="#3329BF", s=30, zorder=3)
        ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
        ax.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

        ax.axhline(y=5, color='g', linestyle='--', label='5 Å Threshold')

        ax.set_xticks(np.arange(1, len(label_names) + 1))
        ax.set_xticklabels(label_names, rotation=45, ha='right', fontsize=7)

        ax.set_xlabel('Protein Pair')
        ax.set_ylabel('Minimum Distance across copies (Å)')
        ax.set_title('Minimum Distance between protein pair')

        plt.tight_layout()
        # plt.show()
        plt.savefig(Path(output_dir) / f"{outname}.png", dpi=600)
        plt.close()

def main(
    input: str,
    xyzr_file: str,
    nproc: int,
    float_dtype: np.dtype,
    binding_data_dir: str,
    outname: str="fit_to_binding_data",
    overwrite: bool=False,
    merge_copies: bool=False,
):

    os.makedirs(binding_data_dir, exist_ok=True)

    binding_data_ = BindingData(
        input=input,
        xyzr_parser=XYZRParser(xyzr_file=xyzr_file),
        merge_copies=merge_copies,
        self_interaction="allow_copies",
        f_dtype=float_dtype,
        i_dtype=np.int32,
    )

    binding_data_.compute_pairwise_distance(
        binding_data_dir=binding_data_dir,
        nproc=nproc,
        overwrite=overwrite,
    )

    binding_data_.plot_pairwise_distances(
        output_dir=binding_data_dir,
        outname=outname,
    )

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the JSON file specifying binding data.",
    )
    parser.add_argument(
        "--xyzr_file",
        # default=f"/home/{_user}/Projects/cardiac_desmosome/analysis/prod_runs/output_trans_dsc_5716/set5/rmf_xyzr_data/sampcon_extracted_frames_xyzr.h5",
        default=f"/data/{_user}/Projects/cardiac_desmosome/analysis/test_runs/analysis_cis_dsc_5176/set2/rmf_xyzr_data/sampcon_extracted_frames_xyzr.h5",
        type=str,
        help="Path to the input HDF5 file containing XYZR data.",
    )
    parser.add_argument(
        "--nproc",
        default=16,
        type=int,
        help="Number of processes for parallel execution.",
    )
    parser.add_argument(
        "--output_dir",
        # default=f"/home/{_user}/Projects/cardiac_desmosome/analysis/prod_runs/output_trans_dsc_5716/set5/fit_to_binding_data",
        default=f"/data/{_user}/Projects/cardiac_desmosome/analysis/test_runs/analysis_cis_dsc_5176/set2/fit_to_binding_data",
        type=str,
        help="Directory to save output files.",
    )
    # parser.add_argument(
    #     "--plotting",
    #     type=str,
    #     default="matplotlib",
    #     help="Plotting options: 'plotly' or 'matplotlib'.",
    # )
    parser.add_argument(
        "--merge_copies",
        action="store_true",
        default=False,
        help="Whether to merge maps across copies of the same protein pair.",
    )
    parser.add_argument(
        "--float_dtype",
        type=int,
        default=64,
        choices=[16, 32, 64],
        help="Float dtype for distance map calculations.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Whether to overwrite existing distance map files."
    )

    args = parser.parse_args()

    ############################################################################

    main(
        input=args.input,
        xyzr_file=args.xyzr_file,
        nproc=args.nproc,
        float_dtype=F_DTYPES.get(args.float_dtype, np.float64),
        binding_data_dir=args.output_dir,
        outname="fit_to_binding_data",
        overwrite=bool(args.overwrite),
        merge_copies=bool(args.merge_copies),
    )