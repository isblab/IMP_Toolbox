```python
def plot_map(contact_map: np.ndarray, chain1: str, chain2: str, p1_region: tuple, p2_region: tuple, plot_type: str):
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
```