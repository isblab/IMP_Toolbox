import sys
import IMP
import IMP.rmf
import RMF
import tqdm
import argparse
import os

def add_frames(
    out_rmf,
    out_hier,
    in_rmf,
    frame_list,
    offset=0
):
    """ Load frames from an RMF file and save them to another RMF file.

    Args:
        out_rmf : The output RMF file to save the frames.
        out_hier : The hierarchy to which the frames belong.
        in_rmf (_type_): The input RMF file from which frames are loaded.
        frame_list (list): A list of frame IDs to be loaded from the input RMF file.
        offset (int, optional): An offset to be applied to the frame IDs in the list.
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
        default=f'/home/omkar/imp_toolbox_test/sampcon_extract/extracted_frames.rmf',
        type=str,
        help="Output RMF file to save the extracted frames."
    )
    parser.add_argument(
        "--rmf1",
        default=f'/home/omkar/imp_toolbox_test/analysis/pmi_analysis/A_models_clust2.rmf3',
        type=str,
        help="Input RMF file 1."
    )
    parser.add_argument(
        "--list1",
        default=f'/home/omkar/imp_toolbox_test/analysis/sampcon_output/analysis/cluster.0.sample_A.txt',
        type=str,
        help="List of frame IDs for RMF file 1."
    )
    parser.add_argument(
        "--rmf2",
        default=f'/home/omkar/imp_toolbox_test/analysis/pmi_analysis/B_models_clust2.rmf3',
        type=str,
        help="Input RMF file 2."
    )
    parser.add_argument(
        "--list2",
        default=f'/home/omkar/imp_toolbox_test/analysis/sampcon_output/analysis/cluster.0.sample_B.txt',
        type=str,
        help="List of frame IDs for RMF file 2."
    )
    args = parser.parse_args()

    # sys.argv -> final output, file1 rmf, file1 list, file2 rmf, file2 list
    # file2 list is assumed to have an offset = size(file1 rmf)

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