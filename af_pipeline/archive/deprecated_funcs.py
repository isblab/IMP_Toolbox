
def get_interacting_patches(
    self,
    contact_map: np.array,
    region_of_interest: dict,
    bandwidth: float = 0, # mod this to 1
):
    """Get the interacting patches and the corrresponding segmented map

    Args:
        contact_map (np.array): a binary contact map
        region_of_interest (dict): interacting region in the format (after af_offset) \n
                                    {chain1: (start1, end1), chain2: (start2, end2)}
        bandwidth (float, optional): bandwidth for the MeanShift clustering. Defaults to None. \n
                                    If not given, the bandwidth is estimated using sklearn.cluster.estimate_bandwidth

    Returns:
        segmented_map (np.array): segmented map of interacting patches
        patches (dict): interacting patches in the format \n
                        {patch_id: {chain_id1: "1-50,375", chain_id2: "2-10"}}
    """

    patches = {}
    interacting_res_pairs = np.argwhere(contact_map == 1)

    if len(interacting_res_pairs) == 0:  # no interacting residues
        return contact_map, patches

    chain1, chain2 = region_of_interest.keys()

    ms = MeanShift(bandwidth=bandwidth)
    ms.fit(interacting_res_pairs)
    labels = ms.labels_
    # cluster_centers = ms.cluster_centers_

    for res_idcs, label in zip(interacting_res_pairs, labels):
        res1_idx, res2_idx = res_idcs
        res1_pos = res1_idx + region_of_interest[chain1][0]
        res2_pos = res2_idx + region_of_interest[chain2][0]

        if label + 1 not in patches:
            patches[int(label + 1)] = {
                chain1: [res1_pos],
                chain2: [res2_pos],
            }

        else:
            patches[int(label + 1)][chain1].append(res1_pos)
            patches[int(label + 1)][chain2].append(res2_pos)

    for label, patch in patches.items():
        for chain in patch.keys():
            patch_residues = list(set(patch[chain]))
            patch[chain] = get_key_from_res_range(patch_residues)

    # segmented map of interacting patches
    segmented_map = np.zeros_like(contact_map)
    for i, res_idx in enumerate(interacting_res_pairs):
        segmented_map[res_idx[0], res_idx[1]] = labels[i] + 1

    return segmented_map, patches
