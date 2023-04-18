import sys
import IMP
import IMP.rmf
import RMF
import tqdm

# sys.argv -> final output, file1 rmf, file1 list, file2 rmf, file2 list
# file2 list is assumed to have an offset = size(file1 rmf)
rmf_out, rmf1, list1, rmf2, list2 = sys.argv[1:]


def add_frames(out_rmf, out_hier, in_rmf, frame_list, offset=0):
    with open(frame_list) as f:
        rd = f.read().strip().split('\n')
    rd = list(map(lambda x: int(x) - offset, rd))
    f = RMF.open_rmf_file_read_only(in_rmf)
    IMP.rmf.link_hierarchies(f, [out_hier])
    for i in tqdm.tqdm(rd, smoothing=0, desc='Loading RMF'):
        IMP.rmf.load_frame(f, RMF.FrameID(i))
        IMP.rmf.save_frame(out_rmf, str(i))


m = IMP.Model()
temp_f = RMF.open_rmf_file_read_only(rmf1)
n_rmf1 = temp_f.get_number_of_frames()
h0 = IMP.rmf.create_hierarchies(temp_f, m)[0]
fh_out = RMF.create_rmf_file(rmf_out)
IMP.rmf.add_hierarchy(fh_out, h0)
del temp_f
add_frames(fh_out, h0, rmf1, list1, 0)
add_frames(fh_out, h0, rmf2, list2, n_rmf1)


