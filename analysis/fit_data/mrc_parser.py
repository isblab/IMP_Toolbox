import numpy as np
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


def calculate_with_grid_of_map1(g0, v0, g1, v1, name):
    interpolator = rgi(g0, v0, bounds_error=False, fill_value=0)
    x, y, z = g1
    x, y, z = np.meshgrid(x, y, z, indexing='ij')
    desired_points = np.array([(i, j, k) for i, j, k in zip(x.flatten(), y.flatten(), z.flatten())])
    values = interpolator(desired_points).reshape(x.shape)
    o, c0, cm, cp = calculate_metrics(values, v1)
    print(f'{name} (all points)')
    print(f'\tOverlap: {o:.4f}, C0: {c0:.4f}, Cam: {cm:.4f} ({cp[0]:.4f}, {cp[1]:.4f})')
    o, c0, cm, cp = calculate_metrics(values, v1, False)
    print(f'{name} (all points)')
    print(f'\tOverlap: {o:.4f}, C0: {c0:.4f}, Cam: {cm:.4f} ({cp[0]:.4f}, {cp[1]:.4f})')


def calculate_with_external_grid_with_addition(list_of_g, list_of_v, name, voxel_spacing=5):
    # First member of the list is the reference, rest are added
    interpolators = []
    for g, v in zip(list_of_g, list_of_v):
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
    if len(values) > 2:
        summed_vals = np.sum(values[1:], axis=0)
    else:
        summed_vals = values[1]
    o, c0, cm, cp = calculate_metrics(values[0], summed_vals)
    print(f'{name} (all points)')
    print(f'\tOverlap: {o:.4f}, C0: {c0:.4f}, Cam: {cm:.4f} ({cp[0]:.4f}, {cp[1]:.4f})')
    o, c0, cm, cp = calculate_metrics(values[0], summed_vals, False)
    print(f'{name} (all points)')
    print(f'\tOverlap: {o:.4f}, C0: {c0:.4f}, Cam: {cm:.4f} ({cp[0]:.4f}, {cp[1]:.4f})')
