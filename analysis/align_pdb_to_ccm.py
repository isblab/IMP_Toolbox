import os
import yaml
import argparse
import IMP
import RMF
import IMP.core
import IMP.rmf
import IMP.atom
import IMP.pmi.analysis
from IMP_Toolbox.utils_imp_toolbox.obj_helpers import (
    get_key_from_res_range,
    get_res_range_from_key,
)

def extract_coordinates(
    ccm_file: str,
    pdb_file: str,
    pdb_chain_id: str,
    protein_name: str,
    copy_index: int,
    selection_range: str,
):
    """ Extract coordinates from CCM and PDB for given selection.

    Args:

        ccm_file (str):
            Path to the cluster center model RMF file.

        pdb_file (str):
            Path to the reference PDB file.

        pdb_chain_id (str):
            Chain ID in the PDB file to select.

        protein_name (str):
            Protein name in the cluster center model to select.

        copy_index (int):
            Copy index of the protein in the cluster center model to select.

        selection_range (str):
            Residue range to select (e.g., "10-50" or "10-20,30-40").

    Returns:
        tuple:
        - coords_ccm (list):
            List of coordinates from the CCM model.
        - coords_pdb (list):
            List of coordinates from the PDB model.
    """

    residue_indices = get_res_range_from_key(selection_range)

    ccm_mdl = IMP.Model()
    ccm = RMF.open_rmf_file_read_only(ccm_file)
    hier = IMP.rmf.create_hierarchies(ccm, ccm_mdl)[0]
    IMP.rmf.load_frame(ccm, 0)
    ccm_mdl.update()

    pdb_ca_mdl = IMP.Model()
    pdb_ca = IMP.atom.read_pdb(
        pdb_file,
        pdb_ca_mdl,
        IMP.atom.CAlphaPDBSelector()
    )
    pdb_ca_mdl.update()

    sel_ca_pdb = IMP.atom.Selection(
        pdb_ca,
        resolution=1,
        chain_id=pdb_chain_id,
        residue_indexes=residue_indices,
    ).get_selected_particles()

    # monkey patch to remove standalone coarse grained beads from CCM selection
    res_nums_in_pdb = [
        int(IMP.atom.Residue(leaf).get_name().split()[-1]) for leaf in sel_ca_pdb
    ]
    mising_res_nums = [
        i for i in residue_indices
        if i not in res_nums_in_pdb
    ]
    standalone_beads = [
        str(b) for b in get_key_from_res_range(mising_res_nums, as_list=True)
        if '-' not in b
    ]

    sel_ccm = IMP.atom.Selection(
        hier,
        resolution=1,
        molecule=protein_name,
        copy_index=copy_index,
        residue_indexes=residue_indices,
    ).get_selected_particles()

    new_sel_ccm = []
    for leaf in sel_ccm:
        if (
            not IMP.atom.Fragment.get_is_setup(leaf)
            and leaf.get_name() not in standalone_beads
        ):
            new_sel_ccm.append(leaf)
        #     print("Kept atomistic bead:", leaf)
        # else:
        #     print("Removed coarse grained bead:", leaf)

    assert len(sel_ca_pdb) == len(new_sel_ccm), (
        f"""Length mismatch: PDB {len(sel_ca_pdb)} vs CCM {len(new_sel_ccm)}\n
        new_sel_ccm: {new_sel_ccm}\n
        sel_ca_pdb: {sel_ca_pdb}
        """
    )

    coords_pdb = [IMP.core.XYZ(i).get_coordinates() for i in sel_ca_pdb]
    coords_ccm = [IMP.core.XYZ(i).get_coordinates() for i in new_sel_ccm]

    return coords_ccm, coords_pdb

def find_transform_and_save(
    query: dict,
    template: dict,
    output_name: str,
    pdb_file: str,
    pdb_chain_selector: IMP.atom.PDBSelector,
):
    """ Find transformation between query and template and save the transformed PDB.

    Args:

        query (dict):
            Query coordinates dictionary with chain IDs as keys.

        template (dict):
            Template coordinates dictionary with chain IDs as keys.

        output_name (str):
            Output PDB file name.

        pdb_file (str):
            Path to the PDB file to be transformed.

        pdb_chain_selector (IMP.atom.PDBSelector):
            PDB selector for selecting chains in the PDB file.
    """

    assert isinstance(query, dict), "Query must be a dictionary"
    assert isinstance(template, dict), "Template must be a dictionary"

    new_mdl = IMP.Model()
    reload = IMP.atom.read_pdb(pdb_file, new_mdl, pdb_chain_selector)
    _tra, trb = IMP.pmi.analysis.Alignment(query=query, template=template).align()

    IMP.atom.transform(reload, trb)
    IMP.atom.write_pdb(reload, output_name)

def fit_pdb_to_ccm(
    ccm_file: str,
    fits_to_perform: list,
    output_dir: str,
):
    """ Fit PDB models to the cluster center model (CCM).

    Args:

        ccm_file (str):
            Path to the RMF file containing the cluster center model.

        fits_to_perform (list):
            List of fit configurations specifying PDB files, chain IDs, molecules,
            selection ranges, and copy indices.

        output_dir (str):
            Path to save the aligned PDB files.
    """

    for fit in fits_to_perform:

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
            output_path = os.path.join(output_dir, output_name)

            chain_selector = IMP.atom.ChainPDBSelector(pdb_chain_ids)

            find_transform_and_save(
                query=coords_pdb_all,
                template=coords_ccm_all,
                output_name=output_path,
                pdb_file=pdb_file,
                pdb_chain_selector=chain_selector,
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

    input_config = yaml.load(open(args.input), Loader=yaml.FullLoader)
    fits_to_perform = input_config["fits_to_perform"]
    # print(fits_to_perform)
    os.makedirs(args.output_dir, exist_ok=True)
    ccm_file = args.ccm_file

    fit_pdb_to_ccm(
        ccm_file=ccm_file,
        fits_to_perform=fits_to_perform,
        output_dir=args.output_dir,
    )