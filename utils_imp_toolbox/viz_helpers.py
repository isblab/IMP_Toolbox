import os
import numpy as np
from itertools import product
import sys
IMP_TOOLBOX = "/home/omkar/Omkar/IMP_TOOLBOX"
sys.path.append(IMP_TOOLBOX)
from IMP_Toolbox.utils_imp_toolbox.obj_helpers import get_key_from_res_range

def make_protein(
    protein_name,
    seq_start,
    seq_end,
    *kwargs
):
    template_dict = {
            "mode": "protein",
            "config": {
            "canvasWidth": seq_end - seq_start + 1,
            "canvasHeight": 600,
            "viewBox": "10 -230 1000 600",
            "scaleLength": 3000,
            "zoomRatio": 1
            },
            "data":[
                {
                "type": "protein",
                "length": seq_end - seq_start + 1,
                "option": {
                    "displayName": protein_name,
                    "position": {
                        "start": {
                            "site": seq_start,
                            "display": "top"
                        },
                        "end": {
                            "site": seq_end,
                            "display": "top"
                        }
                    },
                    "coordinate": {
                        "vertical": {
                            "start": 1,
                            "isLocked": False
                        },
                        "horizontal": {
                            "start": 0,
                            "end": seq_end - seq_start,
                            "isLocked": False
                        }
                    },
                    "id": "Protein_" + protein_name,
                    "style":{
                        "align": "custom",
                        "height": 25,
                        "fontSize": 12,
                        "color": "#babdb6",
                        "gradient": "none",
                        "texture": {
                        "type": "none",
                        "color": "#333333"
                        }
                    },
                    "borderStyle": {
                        "color": "#000000",
                        "size": 1,
                        "isDash": False
                    }
                },
                "children": []
            }
            ]
     }

    return template_dict


def generate_cmap(n, scheme="soft-warm"):
    """ Generate a list of n colors
    (modfied from chatgpt)

    Args:
        n (int): number of colors
        scheme (str, optional): Defaults to "soft-warm".

    Returns:
        colors (list): list of n colors
    """

    import random
    import time

    colors = set()
    start = time.time()

    assert n > 0; "Number of colors must be greater than 0"
    if n != 2:
        assert scheme != "binary"; "Binary scheme is only valid for 2 colors"

    if n == 2 and scheme == "binary":
        return ["black", "green"]

    while len(colors) < n:
        if scheme == "standard":
            r = random.randint(0, 255)
            g = random.randint(0, 255)
            b = random.randint(0, 255)

        elif scheme == "non-bright":
            r = random.randint(50, 180)
            g = random.randint(50, 180)
            b = random.randint(50, 180)

        elif scheme == "earth-tone":
            r = random.randint(100, 180)
            g = random.randint(60, 140)
            b = random.randint(40, 120)

        elif scheme == "cool-tone":
            r = random.randint(50, 120)
            g = random.randint(100, 180)
            b = random.randint(120, 255)

        elif scheme == "soft-warm":
            r = random.randint(180, 255)
            g = random.randint(130, 200)
            b = random.randint(90, 160)

        elif scheme == "contrasting-non-bright":
            if len(colors) % 2 == 0:  # Alternate between darker and lighter muted tones
                r = random.randint(40, 120)
                g = random.randint(40, 120)
                b = random.randint(40, 120)
            else:
                r = random.randint(140, 200)
                g = random.randint(140, 200)
                b = random.randint(140, 200)

        else:
            raise ValueError(
                "Invalid scheme. Choose from 'non-bright', 'earth-tone', 'cool-tone', or 'soft-warm'."
            )

        color = "#{:02x}{:02x}{:02x}".format(r, g, b)
        colors.add(color)

        if time.time() - start > 10:
            break

    return list(colors)


def plot_map(
    contact_map: np.ndarray,
    chain1: str,
    chain2: str,
    p1_region: tuple,
    p2_region: tuple,
    plot_type: str
):
    """Plot the contact map

    Args:
        contact_map (np.ndarray): binary contact map or segmented map with labels
    """

    xtick_vals = np.arange(
        0, p2_region[1] - p2_region[0] + 1
    )
    xtick_labels = [str(x+p2_region[0]) for x in xtick_vals]

    ytick_vals = np.arange(
        0, p1_region[1] - p1_region[0] + 1
    )
    ytick_labels = [str(x+p1_region[0]) for x in ytick_vals]

    num_unique_patches = len(np.unique(contact_map))

    colorscale = generate_cmap(
        n=num_unique_patches,
        scheme="binary" if num_unique_patches == 2 else "soft-warm",
    )

    if plot_type == "interactive":

        import plotly.graph_objects as go


        fig = go.Figure(
            data=go.Heatmap(
                z=contact_map,
                colorscale=colorscale,
                xgap=0.2,
                ygap=0.2,
            )
        )

        fig.update_layout(
            title="Contact Map",
            yaxis_title=f"Residue number of {chain1}",
            xaxis_title=f"Residue number of {chain2}",
            xaxis=dict(
                tickmode="array",
                tickformat=".0f",
                tickvals=xtick_vals,
                ticktext=xtick_labels,
            ),
            yaxis=dict(
                tickmode="array",
                tickformat=".0f",
                tickvals=ytick_vals,
                ticktext=ytick_labels,
            ),
        )

    elif plot_type == "static":

            import matplotlib.pyplot as plt
            import matplotlib.colors as mcolors

            cmap = mcolors.ListedColormap(colorscale)
            fig, ax = plt.subplots()

            fig = ax.imshow(
                contact_map,
                cmap=cmap,
                interpolation="nearest",
            )

            xtick_labels = [xtick_labels[0], xtick_labels[-1]]
            ytick_labels = [ytick_labels[0], ytick_labels[-1]]
            xtick_vals = [xtick_vals[0], xtick_vals[-1]]
            ytick_vals = [ytick_vals[0], ytick_vals[-1]]

            ax.set_xticks(xtick_vals)
            ax.set_xticklabels(xtick_labels)
            ax.set_yticks(ytick_vals)
            ax.set_yticklabels(ytick_labels)

            ax.set_xlabel(f"Residue number of {chain2}")
            ax.set_ylabel(f"Residue number of {chain1}")
            ax.set_title("Contact Map")

    return fig

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

