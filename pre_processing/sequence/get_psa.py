############### Align sequences given Uniprot IDs ################
########## ------>"May the Force serve u well..." <------#########
##################################################################

############# One above all #############
##-------------------------------------##
import time
import numpy as np
import pandas as pd
import sys
import requests

import Bio
from Bio import SeqIO
from io import StringIO
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO

uni_id1 = sys.argv[1]    # Uniprot ID for query sequence.
uni_id2 = sys.argv[2]    # Uniprot ID for target sequence
res_range = sys.argv[3]  # Residue range to be mapped. Input as res1,res2 (None if not required).

class PSA():
	def __init__( self, uni_id1, uni_id2, res_range ):
		self.uni_id1 = uni_id1
		self.uni_id2 = uni_id2
		self.res_range = list( map( int, res_range.split( "," ) ) )

	def forward( self ):
		print( "Downloading Uniprot sequences..." )
		seqa = self.get_uniprot_seq( self.uni_id1 )
		seqb = self.get_uniprot_seq( self.uni_id2 )
		print( "Sequence download complete..." )

		# Using NeedleCommandLine for PSA.
		aligned_seqa, aligned_seqb = self.Needle_alignment( seqa, seqb )

		if self.res_range[0] > len( seqa ) or self.res_range[1] > len( seqa ):
			raise Exception( "Residue range exceeds the length of Protein A..." )
		
		print( "Obtaining aligned positions..." )
		if self.res_range == None:
			exit()
		elif self.res_range[0] == self.res_range[1]:
			aligned_pos = self.get_aligned_pos( aligned_seqb, aligned_seqa, self.res_range[0] )
			print( f"Residue {self.res_range[0]} in protein A aligned to residue {aligned_pos} in protein B." )
		else:
			aligned_pos = { f"{self.uni_id1}":[], f"{self.uni_id2}":[] }
			for res in np.arange( self.res_range[0], self.res_range[1] + 1, 1 ):
				aligned_pos[self.uni_id1].append( res )
				aligned_pos[self.uni_id2].append( self.get_aligned_pos( aligned_seqb, aligned_seqa, res ) )
			df = pd.DataFrame()
			df[self.uni_id1] = aligned_pos[self.uni_id1]
			df[self.uni_id2] = aligned_pos[self.uni_id2]
			df.to_csv( f"Aligned_{self.uni_id1}-{self.uni_id2}.csv", index = False )
		print( "May the Force be with you..." )


	def send_request( self, url, max_trials = 5 ):
		# Send a request to the server to fetch the data.
		# url --> URL of the server from where to fetch the data.
		for i in range( 1, 10 ):
			try:
				response = requests.get( url)
			
			# If there is an error (say Connection Timed Out).
			except Exception as e:
				# Wait for 10 seconds.
				print("There there...\nTry again fella...")
				time.sleep(i*2)
				if i == max_trials-1:
					print("Aah...just leave it man...\n")
					exit()
				# Try mapping them again...
				continue
			else:
				return response.text
				break


	def get_uniprot_seq( self, uni_id ):
		# Obtain Uniprot seq for the specified Uniprot ID.
		# uni_id --> Uniprot accession.
		url = f"http://www.uniprot.org/uniprot/{uni_id}.fasta"
		data = self.send_request( url )
		seq_record = [str( record.seq ) for record in SeqIO.parse( StringIO( data ), 'fasta' )]
		# Some Uniprot IDs might be obsolete.
		if seq_record == []:
			return []
		else:
			return seq_record[0]


	def get_aligned_pos(self, target, query, position ):
		# Return the residue position from the aligned seq for protein2.
		# target --> sequence for which the aligned position is required.
		# query --> sequence to which target is mapped.
		# position --> residue position of query seq for which the aligned target residue position is required.

		iter1 = 0
		idx = 0
		for i, pos in enumerate( query ):
			if pos != "-":
				if iter1 == position:
					idx = i
					break
				iter1 += 1
		

		iter2 = 0
		for j in range( len( target[:idx] ) ):
			# For residues mapped to a gap, 0 will be returned.
			if target[idx] == "-":
				return 0
				break
			elif target[j] != "-":
				iter2 += 1
			else:
				continue
		return iter2


	def Needle_alignment(self, seqa, seqb ):
		# Sequence alignment using EMBOSS Needle.
		# seqa --> seq for protein A.
		# seqb --> seq for protein B.
		outputfile = f"./Aligned_{self.uni_id1}-{self.uni_id2}.txt"
		
		# Define parameters for the alignment.
		needle_cline = NeedleCommandline()
		needle_cline.asequence = "asis:" + seqa
		needle_cline.bsequence = "asis:" + seqb
		needle_cline.gapopen = 10
		needle_cline.gapextend = 0.5
		needle_cline.endopen = 10
		needle_cline.endextend = 0.5
		needle_cline.outfile = outputfile

		stdout, stderr = needle_cline()
		alignments = AlignIO.read( outputfile, "emboss" )

		return "".join( [a for a in alignments[0]] ), "".join([b for b in alignments[1]] )

####################################################
####################################################
PSA( uni_id1, uni_id2, res_range ).forward()