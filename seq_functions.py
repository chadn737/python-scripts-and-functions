import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC123

#interpret sequence context, taken from methylpy.utils
def expand_nucleotide_code(mc_type=["C"]):
    iub_dict = {"N":["A","C","G","T"],"H":["A","C","T"],"C":["C"],"G":["G"],"T":["T"],"A":["A"]}
    for type in mc_type[:]:
        type += "N" * (3 - len(type))
        mc_type.extend(["".join(i) for i in itertools.product(*[iub_dict[nuc] for nuc in type])])
    if "C" in mc_type:
        mc_type.extend(["CG", "CHG", "CHH","CNN"])
    if "CG" in mc_type:
        mc_type.extend(["CGN"])
    return mc_type

#
def codon_GC(fasta,output=()):
	df = pd.DataFrame(columns=["Gene","Length","GC","GC1","GC2","GC3"])
	for sequence in SeqIO.parse(fasta, "fasta"):
		GC = GC123(sequence)
		length = len(sequence.seq)
		df = df.append({"Gene":sequence.id,"Length":length,"GC":GC[0],"GC1":GC[1],"GC2":GC[2],"GC3":GC[3]},ignore_index=True)
	if output:
        	df.to_csv(output, sep='\t', index=False)
	return df

# count the subcontexts in fasta
def count_subcontext_fasta(fasta,context=['CG','CHG','CHH'],output=(),filter_chr=[]):
    df = pd.DataFrame(columns=['context','total_bases'])
    for c in context:
        count = 0
        for i in expand_nucleotide_code([c]):
            for sequence in SeqIO.parse(fasta, "fasta"):
                if sequence.name not in filter_chr:
                    count = count + sequence.seq.count(i) + sequence.seq.reverse_complement().count(i)
        df = df.append({'context': c, 'total_bases': count}, ignore_index=True)
    if output:
        df.to_csv(output, sep='\t', index=False)
    else:
        return df

#
def something(fasta,context,output=(),filter_chr=[]):
	
