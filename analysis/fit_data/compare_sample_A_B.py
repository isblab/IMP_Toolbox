import numpy as np
import os
import sys
import mrcfile
from scipy.interpolate import RegularGridInterpolator as rgi
from scipy.stats import pearsonr, spearmanr


def calculate_metrics(v1, v2, filter_zero=False):
    v1, v2 = v1.flatten(), v2.flatten()
    if filter_zero:
        v1 = v1[v2 != 0]
        v2 = v2[v2 != 0]
    overlap = np.sum(v1 * v2)
    overlap_mean_sub = np.sum((v1 - v1.mean()) * (v2 - v2.mean()))
    mag1 = np.sqrt(np.sum(v1 ** 2))
    mag2 = np.sqrt(np.sum(v2 ** 2))
    mag_meaned1 = np.sqrt(np.sum((v1 - v1.mean()) ** 2))
    mag_meaned2 = np.sqrt(np.sum((v2 - v2.mean()) ** 2))
    corr_over_zero = overlap / mag1 / mag2
    corr_over_mean = overlap_mean_sub / mag_meaned1 / mag_meaned2
    corr_pearson = pearsonr(v1, v2).statistic
    corr_spearman, pval = spearmanr(v1, v2)
    return overlap, corr_over_zero, corr_over_mean, (corr_pearson, corr_spearman)


def calculate_with_external_grid_with_addition_two_lists(list_of_g, list_of_v1, list_of_v2, name, voxel_spacing=5):
    # the two "v" lists are added and compared
    interpolators = []
    for g, v in zip(list_of_g, list_of_v1 + list_of_v2):
        interpolators.append(rgi(g, v, bounds_error=False, fill_value=0))
    xs, ys, zs = [], [], []
    for g in list_of_g:
        xs.append(g[0])
        ys.append(g[1])
        zs.append(g[2])
    xfinal = np.arange(np.min(np.hstack(xs)), np.max(np.hstack(xs)), voxel_spacing)
    yfinal = np.arange(np.min(np.hstack(ys)), np.max(np.hstack(ys)), voxel_spacing)
    zfinal = np.arange(np.min(np.hstack(zs)), np.max(np.hstack(zs)), voxel_spacing)
    xfinal, yfinal, zfinal = np.meshgrid(xfinal, yfinal, zfinal, indexing='ij')
    desired_points = np.array([(i, j, k) for i, j, k in zip(xfinal.flatten(), yfinal.flatten(), zfinal.flatten())])
    values = []
    for i in interpolators:
        values.append(i(desired_points).reshape(xfinal.shape))
    sum1 = np.sum(values[:len(list_of_v1)], axis=0)
    sum2 = np.sum(values[len(list_of_v1):], axis=0)
    o, c0, cm, cp = calculate_metrics(sum1, sum2)
    print(f'{name} (all points)')
    print(f'\tOverlap: {o:.4f}, C0: {c0:.4f}, Cam: {cm:.4f} ({cp[0]:.4f}, {cp[1]:.4f})')
    o, c0, cm, cp = calculate_metrics(sum1, sum2, False)
    print(f'{name} (all points, filter=False)')
    print(f'\tOverlap: {o:.4f}, C0: {c0:.4f}, Cam: {cm:.4f} ({cp[0]:.4f}, {cp[1]:.4f})')


path = sys.argv[1]
pa = f'{path}/Sample_A'
pb = f'{path}/Sample_B'
fa = [f'{pa}/{f}' for f in os.listdir(pa) if 'LPD' in f]
fb = [f'{pb}/{f}' for f in os.listdir(pa) if 'LPD' in f]

mrc_a = [mrcfile.open(x, 'r', permissive=True) for x in fa if 'GPKP' not in x]
mrc_b = [mrcfile.open(x, 'r', permissive=True) for x in fb if 'GPKP' not in x]

grid_points = []
value_points = []
for mrc in mrc_a + mrc_b:
    xvals = mrc.voxel_size.x * np.arange(mrc.data.shape[2]) + mrc.header.origin.x
    yvals = mrc.voxel_size.y * np.arange(mrc.data.shape[1]) + mrc.header.origin.y
    zvals = mrc.voxel_size.z * np.arange(mrc.data.shape[0]) + mrc.header.origin.z
    grid_points.append((xvals, yvals, zvals))
    value_points.append(mrc.data.transpose(2, 1, 0).copy())

calculate_with_external_grid_with_addition_two_lists(grid_points, value_points[:len(mrc_a)],
                                                     value_points[len(mrc_a):], 'Sample A <-> Sample B')
