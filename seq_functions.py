import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC123

def codon_GC(fasta,output=()):
	df = pd.DataFrame(columns=["Gene","Length","GC","GC1","GC2","GC3"])
	for sequence in SeqIO.parse(fasta, "fasta"):
		GC = GC123(sequence)
		length = len(sequence.seq)
		df = df.append({"Gene":sequence.id,"Length":length,"GC":GC[0],"GC1":GC[1],"GC2":GC[2],"GC3":GC[3]},ignore_index=True)
	if output:
        	df.to_csv(output, sep='\t', index=False)
	return df
	
	
	
