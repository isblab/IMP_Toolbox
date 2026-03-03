import IMP
import IMP.rmf
import IMP.atom
import RMF
import tqdm
import argparse
import os
import getpass
_user = getpass.getuser()

def add_frames(
    out_rmf: str,
    out_hier: IMP.atom.Hierarchy,
    in_rmf: str,
    frame_list: str,
    offset: int = 0,
):
    """ Load frames from an RMF file and save them to another RMF file.

    ## Arguments:

    - **out_rmf (str)**:<br />
        Output RMF file to save the extracted frames.

    - **out_hier (IMP.atom.Hierarchy)**:<br />
        Hierarchy to which the frames belong in the output RMF file.

    - **in_rmf (str)**:<br />
        Input RMF file from which frames are loaded.

    - **frame_list (str)**:<br />
        List of frame IDs to be loaded from input RMF file.

    - **offset (int, optional):**:<br />
        Offset to be applied to frame IDs in the list. Default is 0.
    """

    with open(frame_list) as f:
        rd = f.read().strip().split('\n')

    rd = list(map(lambda x: int(x) - offset, rd))
    f = RMF.open_rmf_file_read_only(in_rmf)
    IMP.rmf.link_hierarchies(f, [out_hier])

    for i in tqdm.tqdm(rd, smoothing=0, desc='Loading RMF'):
        IMP.rmf.load_frame(f, RMF.FrameID(i))
        IMP.rmf.save_frame(out_rmf, str(i))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract frames from RMF files based on a list of frame IDs."
    )
    parser.add_argument(
        "--rmf_out",
        default=f"/home/{_user}/imp_toolbox_test/sampcon_extract/extracted_frames.rmf",
        type=str,
        help="Output RMF file to save the extracted frames."
    )
    parser.add_argument(
        "--rmf1",
        default=f"/home/{_user}/imp_toolbox_test/analysis/pmi_analysis/A_models_clust2.rmf3",
        type=str,
        help="Input RMF file 1."
    )
    parser.add_argument(
        "--list1",
        default=f"/home/{_user}/imp_toolbox_test/analysis/sampcon_output/analysis/cluster.0.sample_A.txt",
        type=str,
        help="List of frame IDs for RMF file 1."
    )
    parser.add_argument(
        "--rmf2",
        default=f"/home/{_user}/imp_toolbox_test/analysis/pmi_analysis/B_models_clust2.rmf3",
        type=str,
        help="Input RMF file 2."
    )
    parser.add_argument(
        "--list2",
        default=f"/home/{_user}/imp_toolbox_test/analysis/sampcon_output/analysis/cluster.0.sample_B.txt",
        type=str,
        help="List of frame IDs for RMF file 2."
    )
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.rmf_out), exist_ok=True)

    m = IMP.Model()
    temp_f = RMF.open_rmf_file_read_only(args.rmf1)
    num_frames_rmf1 = temp_f.get_number_of_frames()
    h0 = IMP.rmf.create_hierarchies(temp_f, m)[0]
    fh_out = RMF.create_rmf_file(args.rmf_out)
    IMP.rmf.add_hierarchy(fh_out, h0)
    del temp_f

    add_frames(fh_out, h0, args.rmf1, args.list1, 0)
    add_frames(fh_out, h0, args.rmf2, args.list2, num_frames_rmf1)