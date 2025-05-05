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
        self.pae = self.dataparser.get_modified_pae(data=data)
        self.avg_pae = self.dataparser.get_avg_pae(pae=self.pae)
        self.contact_probs_mat = self.dataparser.get_contact_probs_mat(data=data)
        if self.contact_probs_mat is not None:
            self.avg_contact_probs_mat = self.dataparser.get_avg_contact_probs_mat(contact_probs_mat=self.contact_probs_mat)

        if self.struct_file_path:
            self.token_chain_ids, self.token_res_ids = self.structureparser.get_token_chain_res_ids()
            # Chain lengths for each chain.
            self.lengths_dict = self.structureparser.get_chain_lengths(self.token_chain_ids)
            # Ca-coords of all residues for each chain.
            self.coords_list = self.structureparser.get_ca_coordinates()
            # Ca-plddt of all residues for each chain.
            self.plddt_list = self.structureparser.get_ca_plddt()

        else:
            self.token_chain_ids = self.dataparser.get_token_chain_ids(data=data)
            self.token_res_ids = self.dataparser.get_token_res_ids(data=data)
            self.lengths_dict = self.dataparser.get_chain_lengths(data=data)

        self.sanity_check()
        self.idx_to_num, self.num_to_idx = self.renumber.residue_map(
            self.token_chain_ids, self.token_res_ids
        )


    def sanity_check(self):
        """
        Perform sanity checks on the input data. \n
        If structure file path is not provided, the input data file should be in AF3 format.
        """

        error_statement = "Input data file needs to be in AF3 format if structure path is not provided."

        if not self.lengths_dict:
            raise Exception(f"No chain lengths found. {error_statement}")