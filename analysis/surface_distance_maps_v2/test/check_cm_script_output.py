import IMP
import RMF
import IMP.rmf
import IMP.atom
import IMP.core
import IMP.algebra
import random
import numpy as np

mdl = IMP.Model()
rmf_fname = "./test.rmf3"
rmf_fh = RMF.open_rmf_file_read_only(rmf_fname)
hier = IMP.rmf.create_hierarchies(rmf_fh, mdl)

prot_name = "WDR76"
cp_id = 0
prot_size = 625
r1 = 1
r2 = 626

distances = {}
for frame in range(rmf_fh.get_number_of_frames()):
    IMP.rmf.load_frame(rmf_fh, frame)
    s1 = IMP.atom.Selection(
        hierarchy=hier[0],
        molecule=prot_name,
        copy_index=cp_id,
        residue_index=r1,
        resolution=1,
    ).get_selected_particles()[0]
    s2 = IMP.atom.Selection(
        hierarchy=hier[0],
        molecule=prot_name,
        copy_index=cp_id,
        residue_index=r2,
        resolution=1,
    ).get_selected_particles()[0]

    distance = IMP.core.get_distance(IMP.core.XYZR(s1), IMP.core.XYZR(s2))
    if (r1, r2) not in distances:
        distances[(r1, r2)] = [distance]
    else:
        distances[(r1, r2)].append(distance)


print(IMP.core.XYZR(s1), IMP.core.XYZR(s2))
print(distances)

for k, val in distances.items():
    print(k, np.mean(val))
dmat = np.loadtxt("distance_matrices/WDR76.0-WDR76.0_mean_distances.csv")
print(dmat[r1 - 1, r2 - 1])
