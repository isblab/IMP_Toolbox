from af_pipeline.Parser import AfParser


class _Initialize(AfParser):
    """ Initialize the AF2/3 data and structure files.
    """

    def __init__(
        self,
        data_file_path: str,
        struct_file_path: str,
        af_offset: dict | None = None,
        **kwargs,
    ):

        super().__init__(
            struct_file_path=struct_file_path,
            data_file_path=data_file_path,
            af_offset=af_offset,
            **kwargs,
        )

        # AF2/3 structure file path.
        self.struct_file_path = struct_file_path
        # AF2/3 structure data path.
        self.data_file_path = data_file_path
        self.af_offset = af_offset
        self.average_atom_pae = kwargs.get("average_atom_pae", False)
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
        self.contact_probs_mat = self.dataparser.get_contact_probs_mat(data=data)

        if self.contact_probs_mat is not None:
            self.avg_contact_probs_mat = self.dataparser.get_avg_contact_probs_mat(
                contact_probs_mat=self.contact_probs_mat
            )

        self.token_chain_ids = self.dataparser.get_token_chain_ids(data=data)
        self.token_res_ids = self.dataparser.get_token_res_ids(data=data)

        if self.token_chain_ids is None or self.token_res_ids is None:
            self.token_chain_ids, self.token_res_ids = self.structureparser.get_token_chain_res_ids()

        # Ca-coords of all residues for each chain.
        self.coords_list = self.structureparser.get_ca_coordinates()
        # Ca-plddt of all residues for each chain.
        self.plddt_list = self.structureparser.get_ca_plddt()

        self.pae = self.update_pae(
            pae=self.pae,
            token_res_ids=self.token_res_ids,
            token_chain_ids=self.token_chain_ids,
            average_atom_pae=self.average_atom_pae,
        )
        self.avg_pae = self.update_pae(
            pae=self.avg_pae,
            token_res_ids=self.token_res_ids,
            token_chain_ids=self.token_chain_ids,
            average_atom_pae=self.average_atom_pae,
        )

        if self.contact_probs_mat is not None:
            self.contact_probs_mat = self.update_contact_probs(
                contact_probs_mat=self.contact_probs_mat,
                token_chain_ids=self.token_chain_ids,
                token_res_ids=self.token_res_ids,
                average_atom_pae=self.average_atom_pae,
            )
            self.avg_contact_probs_mat = self.update_contact_probs(
                contact_probs_mat=self.avg_contact_probs_mat,
                token_chain_ids=self.token_chain_ids,
                token_res_ids=self.token_res_ids,
                average_atom_pae=self.average_atom_pae,
            )

        self.token_chain_ids, self.token_res_ids = self.update_token_ids(
            token_chain_ids=self.token_chain_ids,
            token_res_ids=self.token_res_ids,
            average_atom_pae=self.average_atom_pae,
        )

        self.lengths_dict = self.get_chain_lengths(
            token_chain_ids=self.token_chain_ids,
        )
        self.sanity_check()
        self.idx_to_num, self.num_to_idx = self.renumber.residue_map(
            token_chain_ids=self.token_chain_ids,
            token_res_ids=self.token_res_ids,
        )


    def sanity_check(self):
        """
        Perform sanity checks on the input data. \n
        If structure file path is not provided, the input data file should be in AF3 format.
        """

        error_statement = "Input data file needs to be in AF3 format if structure path is not provided."

        if not self.lengths_dict:
            raise Exception(f"No chain lengths found. {error_statement}")