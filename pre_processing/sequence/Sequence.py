from utils import request_session

class FetchSequences:
    """Fetch sequences for given uniprot ids and save them in fasta format
    """

    def __init__(self, uniprot_ids):
        self.uniprot_ids = uniprot_ids

    def uniprot_to_sequences(self, max_retries=3):
        """Get sequences for given uniprot ids in fasta format

        Args:
            uniprot_ids (list): list of uniprot ids
            max_retries (int, optional): Defaults to 3.
        """

        UniProt_API = "https://rest.uniprot.org/"
        joined = ",".join(self.uniprot_ids)
        req_sess = request_session(max_retries=max_retries)
        response = req_sess.get(
            f"{UniProt_API}/uniprotkb/accessions?accessions={joined}&format=fasta"
        )

        if response.status_code == 200:
            print("Successfully fetched sequences for given uniprot ids")
            return response.text

        else:
            raise Exception("Error while requesting sequences for given uniprot ids")


    def only_uniprot_id_as_name(self, fasta=None):
        """Keep only uniprot id as name in fasta file

        Args:
            fasta (str): fasta file

        Returns:
            str: fasta file with only uniprot id as name
        """

        if fasta is None:
            fasta = self.uniprot_to_sequences()

        sequences_lines = fasta.split("\n")
        for i, line in enumerate(sequences_lines):
            if line.startswith(">"):
                sequences_lines[i] = f">{line.split("|")[1]}"

        fasta = "\n".join(sequences_lines)

        return fasta