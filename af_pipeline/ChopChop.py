import numpy as np
from af_pipeline._Initialize import _Initialize
from utils import get_key_from_res_range


class FreeBirds(_Initialize):

    def __init__(
        self,
        data_path: str,
        monomer_pred_dir: str,
        structure_path: str,
        af_offset: dict | None = None,
        monomer_to_chain_map: dict | None = None,
    ):

        super().__init__(
            data_file_path=data_path,
            struct_file_path=structure_path,
            af_offset=af_offset,
        )

        self.data_path = data_path
        self.monomer_pred_dir = monomer_pred_dir
        self.structure_path = structure_path
        self.af_offset = af_offset
        self.pae_cutoff = 5
        self.monomer_to_chain_map = monomer_to_chain_map
        self.pae

    def get_non_confident_residues(self, chain_id: str, num_confident_pairs: int=0):
        """ Get the non-confident residues in a chain based on the PAE matrix. \n
        These residues have PAE < pae_cutoff with at most `num_confident_pairs` \n
        residues in the prediction excluding the residues in the chain

        Args:
            chain_id (str): _description_
            num_confident_pairs (int): _description_
        """

        avg_pae_mat = self.avg_pae
        chain_idx = [
            idx for idx, chain in enumerate(self.token_chain_ids) if chain == chain_id
        ]
        chain_pae = [avg_pae_mat[idx] for idx in chain_idx]
        chain_pae = chain_pae[0 : self.lengths_dict[chain_id]]
        chain_pae = np.array(chain_pae)

        chain_pae[chain_pae < self.pae_cutoff] = 1
        chain_pae[chain_pae >= self.pae_cutoff] = 0

        c_start, c_end = self.af_offset[chain_id]
        res_num = list(np.arange(c_start, c_end + 1))
        # print(res_num)
        res_idx = [self.num_to_idx[chain_id][res] for res in res_num]
        # print(res_idx)
        chain_pae_row_sums = np.sum(chain_pae, axis=1).astype(int)
        # print(chain_pae_row_sums)
        self_chain_pae_row_sums = np.sum(chain_pae[:, res_idx], axis=1).astype(int)
        # print(self_chain_pae_row_sums)
        chain_pae_row_sums -= self_chain_pae_row_sums
        # print(chain_pae_row_sums)
        non_confident_residues = [
            res_num[idx]
            for idx, row_sum in enumerate(chain_pae_row_sums)
            if row_sum <= num_confident_pairs
        ]

        non_confident_residues = get_key_from_res_range(
            non_confident_residues, as_list=True
        )

        print("You should remove these residues from the complex prediction:")
        print(non_confident_residues)

        from plotly import graph_objects as go

        fig = go.Figure()
        fig.add_trace(go.Heatmap(z=chain_pae, colorscale="hot"))
        for patch_to_remove in non_confident_residues:
            if "-" in patch_to_remove:
                fig.add_shape(
                    type="rect",
                    x0=0,
                    y0=int(patch_to_remove.split("-")[0]) - c_start,
                    x1=chain_pae.shape[1],
                    y1=int(patch_to_remove.split("-")[1]) - c_start,
                    line=dict(color="red", width=2),
                    fillcolor="red",
                    opacity=0.3,
                )
        fig.add_shape(
            type="rect",
            x0=self.num_to_idx[chain_id][c_start],
            y0=0,
            x1=self.num_to_idx[chain_id][c_end],
            y1=chain_pae.shape[0],
            line=dict(color="green", width=2),
            fillcolor="green",
            opacity=0.3,
        )
        fig.show()

        return non_confident_residues

    # def get_chain_monomer_plddt(self, chain_id: str):
    #     monomer_pae_path, monomer_structure_path = self.get_chain_monomer_paths(
    #         chain_id
    #     )
    #     chain_instance = AfParser(
    #         data_file_path=monomer_pae_path, struct_file_path=monomer_structure_path
    #     )
    #     plddt_list = chain_instance.structureparser.get_ca_plddt()
    #     return plddt_list

    # def get_chain_monomer_paths(self, chain_id: str):
    #     """ Get the monomer prediction paths for a chain.

    #     Args:
    #         chain_id (str): _description_

    #     Returns:
    #         _type_: _description_
    #     """        
    #     chain_pred_dir = os.path.join(
    #         self.monomer_pred_dir, self.monomer_to_chain_map[chain_id]
    #     )

    #     monomer_pae_path = [
    #         os.path.join(chain_pred_dir, f)
    #         for f in os.listdir(chain_pred_dir)
    #         if f.endswith("full_data_0.json")
    #     ][0]

    #     monomer_structure_path = [
    #         os.path.join(chain_pred_dir, f)
    #         for f in os.listdir(chain_pred_dir)
    #         if f.endswith("model_0.cif")
    #     ][0]

    #     return monomer_pae_path, monomer_structure_path

    # def compare_plddt(self):
    #     unique_chain_ids = list(set(self.token_chain_ids))
    #     for chain_id in unique_chain_ids:
    #         monomer_chain_plddt = self.get_chain_monomer_plddt(chain_id)
    #         region_of_interest_monomer_chain_plddt = monomer_chain_plddt[
    #             self.af_offset[chain_id][0] : self.af_offset[chain_id][1] + 1
    #         ]
    #         complex_chain_res_idx = [
    #             idx
    #             for idx, chain in enumerate(self.token_chain_ids)
    #             if chain == chain_id
    #         ]
    #         complex_chain_plddt = [
    #             self.plddt_list[idx] for idx in complex_chain_res_idx
    #         ]

    #         print(f"Chain ID: {chain_id}")
    #         # print(f"Monomer: {monomer_chain_plddt}")
    #         # print(f"Complex: {complex_chain_plddt}")
    #         figure = plt.figure()
    #         plt.plot(region_of_interest_monomer_chain_plddt, label="Monomer")
    #         plt.plot(complex_chain_plddt, label="Complex")
    #         plt.legend()
    #         plt.title(
    #             f"Chain ID: {chain_id} ({self.monomer_to_chain_map[chain_id].split("_")[0]}:{str(self.af_offset[chain_id][0])}-{str(self.af_offset[chain_id][1])}) pLDDT"
    #         )
    #         plt.show()
    #         # plt.close()

    # def plot_pae_histogram(self):
    #     # pae_mat = self.pae
    #     avg_pae_mat = self.avg_pae
    #     # pae_mat = pae_mat[pae_mat < self.pae_cutoff]
    #     import numpy as np

    #     # plot chain-wise pae histogram
    #     for chain_id in list(set(self.token_chain_ids)):
    #         chain_idx = [
    #             idx
    #             for idx, chain in enumerate(self.token_chain_ids)
    #             if chain == chain_id
    #         ]
    #         chain_pae = [avg_pae_mat[idx] for idx in chain_idx]
    #         chain_pae = chain_pae[0 : self.lengths_dict[chain_id]]
    #         chain_pae = np.array(chain_pae)
    #         # fig = plt.figure()
    #         # plt.hist(chain_pae, bins=30)
    #         # plt.title(f"Chain {chain_id} PAE Histogram")
    #         # plt.show()

    #         fig = plt.figure()
    #         chain_pae[chain_pae < self.pae_cutoff] = 0
    #         chain_pae[chain_pae >= self.pae_cutoff] = 1
    #         plt.imshow(chain_pae, cmap="hot", interpolation="nearest")
    #         plt.title(
    #             f"Chain {chain_id} PAE Map {self.monomer_to_chain_map[chain_id].split('_')[0]}"
    #         )

        # fig = plt.figure()
        # plt.hist(avg_pae_mat.flatten(), bins=30)
        # plt.title("PAE Histogram")
        # plt.show()

        # fig = plt.figure()
        # # convert pae_map to binary
        # avg_pae_mat[avg_pae_mat < self.pae_cutoff] = 0
        # avg_pae_mat[avg_pae_mat >= self.pae_cutoff] = 1
        # plt.imshow(avg_pae_mat, cmap="hot", interpolation="nearest")
