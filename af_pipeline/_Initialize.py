from collections import defaultdict
from af_pipeline.Parser import AfParser


class _Initialize(AfParser):
    """ Initialize the AF2/3 data and structure files.
    """

    def __init__(
        self,
        data_file_path: str,
        struct_file_path: str | None = None,
        af_offset: dict | None = None,
    ):

        super().__init__(
            struct_file_path=struct_file_path,
            data_file_path=data_file_path,
            af_offset=af_offset,
        )

        # AF2/3 structure file path.
        self.struct_file_path = struct_file_path
        # AF2/3 structure data path.
        self.data_file_path = data_file_path
        self.af_offset = af_offset
        self.get_attributes()


    def get_attributes(self):
        """
        Extract the following from the input files:
            0. Residue positions of all residues for each chain.
            1. Ca or representative atom coordinates.
            2. Ca or representative atom pLDDT.
            3. Tokenized chain IDs.
            4. Tokenized residue IDs.
            5. Chain lengths.
            6. PAE matrix.
            7. Average PAE matrix.
            8. Min PAE for each residue.
        """

        data = self.dataparser.get_data_dict()
        self.pae = self.dataparser.get_pae(data=data)
        self.avg_pae = self.dataparser.get_avg_pae(pae=self.pae)

        if self.struct_file_path:
            # Residue positions of all residues for each chain.
            self.res_dict = self.structureparser.get_residue_positions()
            self.token_chain_ids = self.structureparser.get_token_chain_ids(
                self.res_dict
            )
            # self.token_res_ids = self.structureparser.get_token_res_ids(self.res_dict)
            self.lengths_dict = self.structureparser.get_chain_lengths(self.res_dict)
            # Ca-coords of all residues for each chain.
            self.coords_dict = self.structureparser.get_ca_coordinates()
            # Ca-plddt of all residues for each chain.
            self.plddt_dict = self.structureparser.get_ca_plddt()
            # Get minPAE for each residue.
            self.min_pae_dict = self.get_min_pae(
                avg_pae=self.avg_pae,
                lengths_dict=self.lengths_dict,
                mask_intrachain=True,
                return_dict=True,
            )

        else:
            self.token_chain_ids = self.dataparser.get_token_chain_ids(data=data)
            # self.token_res_ids = self.dataparser.get_token_res_ids(data=data)
            self.res_dict = self.dataparser.get_residue_positions(data=data)
            self.lengths_dict = self.dataparser.get_chain_lengths(data=data)

        self.sanity_check()
        self.idx_to_num, self.num_to_idx = self.residue_map(af_offset=self.af_offset)


    def sanity_check(self):
        """
        Perform sanity checks on the input data. \n
        If structure file path is not provided, the input data file should be in AF3 format.
        """

        error_statement = "Input data file needs to be in AF3 format if structure path is not provided."

        if not self.token_chain_ids:
            raise Exception(f"No chain IDs found. {error_statement}")

        # if not self.token_res_ids:
        #     raise Exception(f"No residue IDs found. {error_statement}")

        if not self.lengths_dict:
            raise Exception(f"No chain lengths found. {error_statement}")


    def residue_map(self, af_offset: dict | None = None):
        """
        Create a mapping of residue indices to residue numbers and vice-versa. \n
        res_idx is essentially token index. \n
        res_num is the residue number. \n
        res_num = res_idx + 1 if af_offset is not provided. \n
        res_num = res_idx + af_offset if af_offset is provided. \n
        af_offset informs what is the starting residue number for each chain.
        """

        # res_num = 1

        idx_to_num = defaultdict(dict)
        num_to_idx = defaultdict(dict)

        # for res_idx, chain_id in enumerate(self.token_chain_ids):

        #     res_num = self.token_res_ids[res_idx]

        #     if af_offset and chain_id in af_offset:
        #         res_num += af_offset[chain_id][0] - 1

        #     idx_to_num[chain_id][res_idx] = res_num
        #     num_to_idx[chain_id][res_num] = res_idx

        res_idx = 0
        for chain_id, res_pos_list in self.res_dict.items():
            for res_pos in res_pos_list:
                res_num = res_pos[0]
                if af_offset and chain_id in af_offset:
                    res_num += af_offset[chain_id][0] - 1
                idx_to_num[chain_id][res_idx] = res_num
                num_to_idx[chain_id][res_num] = res_idx
                res_idx += 1

        return idx_to_num, num_to_idx