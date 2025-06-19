```python
def save_map(
    contact_map: np.ndarray,
    avg_contact_probs_mat: np.ndarray | None,
    patches: dict,
    chain1: str,
    chain2: str,
    p1_region: tuple,
    p2_region: tuple,
    out_file: str,
    save_plot=False,
    plot_type="static",
    p1_name: str | None = None,
    p2_name: str | None = None,
    concat_residues: bool = True,
    contact_probability: bool = False,
    num_to_idx: dict = None,
    idx_to_num: dict = None,
):
    """Save the interacting patches and the contact map to a file.

    Args:
        contact_map (np.ndarray): binary contact map or contact map
        avg_contact_probs_mat (np.ndarray): average contact_probs_mat map
        patches (dict): interacting patches from the map
        interacting_region (dict): interacting region specified by the user
        out_file (str): path to save the output file
        save_plot (bool, optional): save the plot. Defaults to False.
    """

    if contact_probability:
        assert avg_contact_probs_mat is not None; "avg_contact_probs_mat must be provided if contact_probability is True"

    out_dir = os.path.dirname(out_file)
    file_name = os.path.basename(out_file).split(".")[0]

    # txt_outfile = os.path.join(out_dir, f"{file_name}.txt")

    # print(f"Writing interacting patches to {txt_outfile}")

    # # with open(txt_outfile, "w") as f:
    # #     for patch_id, patch in patches.items():
    # #             f.write(f"Patch {patch_id}\n")
    # #             for chain, res_range in patch.items():
    # #                 f.write(f"{chain}: {res_range}\n")

    csv_outfile = os.path.join(out_dir, f"{file_name}.csv")

    print(f"Writing interacting patches to {csv_outfile}")

    import pandas as pd
    df_rows = []
    for _, patch in patches.items():
        ch1_res_range = patch[chain1].tolist()
        ch2_res_range = patch[chain2].tolist()

        if concat_residues:
            if contact_probability:

                # res1_idxs = patch[chain1] - p1_region[0]
                # res2_idxs = patch[chain2] - p2_region[0] + p1_region[1] - p1_region[0] + 1

                res1_idxs = np.array([num_to_idx[chain1][res_num] for res_num in ch1_res_range])
                res2_idxs = np.array([num_to_idx[chain2][res_num] for res_num in ch2_res_range])

                contact_probs_mat_res_range = avg_contact_probs_mat[np.ix_(res1_idxs, res2_idxs)]
                avg_contact_prob = np.round(np.mean(contact_probs_mat_res_range), 2)

            ch1_res_range = get_key_from_res_range(ch1_res_range)
            ch2_res_range = get_key_from_res_range(ch2_res_range)

            if contact_probability:
                df_rows.append([ch1_res_range, ch2_res_range, avg_contact_prob])

            else:
                df_rows.append([ch1_res_range, ch2_res_range])

        else:
            for res1, res2 in product(ch1_res_range, ch2_res_range):
                if contact_probability:

                    res1_idx = res1 - p1_region[0]
                    res2_idx = res2 - p2_region[0] + p1_region[1] - p1_region[0] + 1

                    contact_probs_mat_res_range = avg_contact_probs_mat[res1_idx, res2_idx]
                    avg_contact_prob = np.mean(contact_probs_mat_res_range)

                    df_rows.append([res1, res2, avg_contact_prob])

                else:
                    df_rows.append([res1, res2])

    column_names = [f"{chain1}", f"{chain2}"]

    if contact_probability:
        column_names.append("avg_contact_probability")

    if p1_name and p2_name:
        column_names[0] = f"{p1_name}_{chain1}"
        column_names[1] = f"{p2_name}_{chain2}"

    df = pd.DataFrame(df_rows, columns=column_names)

    df.to_csv(csv_outfile, index=False, header=True, sep=",")

    if save_plot:
        if plot_type == "interactive" or plot_type == "both":
            from utils import plot_map
            fig = plot_map(
                contact_map=contact_map,
                chain1=chain1 if p1_name is None else f"{p1_name}_{chain1}",
                chain2=chain2 if p2_name is None else f"{p2_name}_{chain2}",
                p1_region=p1_region,
                p2_region=p2_region,
                plot_type="interactive",
            )
            out_file = os.path.join(out_dir, f"{file_name}.html")
            fig.write_html(
                out_file,
                full_html=False,
            )

            if plot_type == "both":
                plot_type = "static"

        if plot_type == "static":

            from utils import plot_map
            fig = plot_map(
                contact_map=contact_map,
                chain1=chain1 if p1_name is None else f"{p1_name}_{chain1}",
                chain2=chain2 if p2_name is None else f"{p2_name}_{chain2}",
                p1_region=p1_region,
                p2_region=p2_region,
                plot_type="static",
            )

            out_file = os.path.join(out_dir, f"{file_name}.png")
            fig.figure.savefig(out_file)
```