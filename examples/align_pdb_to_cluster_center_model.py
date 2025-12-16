import os
import yaml
import argparse
import IMP
import IMP.atom
from IMP_Toolbox.analysis.align_pdb_to_ccm import (
    extract_coordinates,
    find_transform_and_save,
)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Align a PDB model to the cluster center model (CCM)."
    )
    parser.add_argument(
        "--ccm_file",
        type=str,
        required=True,
        help="Path to the RMF file containing the cluster center model.",
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the YAML configuration file specifying chain and protein details.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Path to save the aligned PDB files.",
    )

    args = parser.parse_args()

    fits_to_transfrom = yaml.load(open(args.input), Loader=yaml.FullLoader)
    # print(fits_to_transfrom)
    os.makedirs(args.output_dir, exist_ok=True)
    ccm_file = args.ccm_file

    for fit in fits_to_transfrom:

        pdb_file = fit["pdb_file"]
        pdb_chain_ids = fit["pdb_chain_ids"]
        molecules = fit["molecules"]
        selection_ranges = fit["selection_ranges"]
        copy_indices_meta = fit["copy_indices"]

        for copy_indices in copy_indices_meta:

            assert len(copy_indices) == len(molecules), (
                f"""Copy indices length must match molecules length. Instead got:
                molecules ({len(molecules)}): {molecules}
                copy_indices ({len(copy_indices)}): {copy_indices}"""
            )

            coords_ccm_all = {}
            coords_pdb_all = {}

            for i, molecule in enumerate(molecules):

                pdb_chain_id = pdb_chain_ids[i]
                selection_range = selection_ranges[i]
                coords_ccm, coords_pdb = extract_coordinates(
                    ccm_file,
                    pdb_file,
                    pdb_chain_id,
                    molecule,
                    copy_indices[i],
                    selection_range,
                )
                coords_ccm_all[molecule] = coords_ccm
                coords_pdb_all[molecule] = coords_pdb

            _zipped_names = list(zip(molecules, copy_indices))
            outname = "_".join([f"{name}_{ci}" for name, ci in _zipped_names])
            output_name = f"aligned_{outname}.pdb"
            output_path = os.path.join(args.output_dir, output_name)

            chain_selector = IMP.atom.ChainPDBSelector(pdb_chain_ids)

            find_transform_and_save(
                query=coords_pdb_all,
                template=coords_ccm_all,
                output_name=output_path,
                pdb_file=pdb_file,
                pdb_chain_selector=chain_selector,
            )