import sys
import pandas as pd
import pybedtools as pbt
import math
from scipy import stats
import itertools
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import Counter 

# Print iterations progress
def printProgress (iteration, total, prefix = '', suffix = '', decimals = 1, barLength = 100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    formatStr       = "{0:." + str(decimals) + "f}"
    percents        = formatStr.format(100 * (iteration / float(total)))
    filledLength    = int(round(barLength * iteration / float(total)))
    bar             = '█' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()

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

#filter allc file based on sequence context
def filter_context(allc,context=["C"]):
    a = pd.read_table(allc)
    a = a[a.mc_class.isin(expand_nucleotide_code(context))]
    return a

#converts allc file to special bed file format for use with pbt
#use: allc2bed("allc file",bed=True)
#if bed=False, returns a pandas dataframe, not a bed file
def allc2bed(allc,context=["C"],bed=True):
    a = filter_context(allc,context)
    a['pos2'] = a.pos
    a['name'] = a.index 
    a['score'] = "."
    a = a[['chr','pos','pos2','name','score','strand','mc_class','mc_count','total','methylated']]
    if bed is True:
        a = pbt.BedTool.from_dataframe(a)
    return a

#map cytosines from allc to features in bed file
#use: map_allc2feature( "allc file", "bed file of features", accend=(sort by ascending or decending, default is accend/True))
#outputs tab-delimited file: column 1 = feature name, column 2 = sequence context, column 3 = methylated reads, 
#column 4 = total reads, column 5 = is site methylated or unmethylated (determined by methylpy binomial test, 1 = yes, 0 = no)
def map_allc2feature(allc,features_bed,ascend=True):
    mC_bed = allc2bed(allc)
    f_bed = pbt.BedTool(features_bed)
    mapping = pbt.bedtool.BedTool.intersect(mC_bed,f_bed,wa=True,wb=True)
    m = pd.read_table(mapping.fn, header=None, usecols = [13,6,7,8,9])
    m = m[[13,6,7,8,9]]
    m = m.sort_values(by = 13,ascending=ascend)
    return m

#calculate the genome-wide weighted methylation using plant specific contexts
#use: weighted_mC_plant("input allc", "output file", <cutoff for coverage needed to include site, default=0>)
#outputs tab-delimited file: column 1 = sequence context, column 2 = total number reads mapping, 
#column 3 = total methylated reads, column 4 = weighted methylation (decimal), column 5 = weighted methylation (as %)
def weighted_mC_plant(allc, mC_results, cutoff=0):
    CG = mCG = CHG = mCHG = CHH = mCHH = CNN = mCNN = 0
    with open(allc) as f:
        next(f)
        for l in f:
            c = l.split('\t')
            if int(c[5]) >= int(cutoff):
                if c[3].startswith("CG"):
                    CG = CG + int(c[5])
                    mCG = mCG + int(c[4])
                elif c[3].endswith("G"):
                    CHG = CHG + int(c[5])
                    mCHG = mCHG + int(c[4])
                elif c[3].startswith("CN") or c[3].endswith("N"):
                    CNN = CNN + int(c[5])
                    mCNN = mCNN + int(c[4])
                else:
                    CHH = CHH + int(c[5])
                    mCHH = mCHH + int(c[4])     
    with open(mC_results, "w") as out:
        out.write('\t'.join(["Context","Total","Methylated","Weighted_mC","%Weighted_mC"]) + '\n')
        out.write('\t'.join(["CG",str(CG),str(mCG),str(np.float64(mCG)/np.float64(CG)),str(pCG*100)]) + '\n')
        out.write('\t'.join(["CHG",str(CHG),str(mCHG),str(np.float64(mCHG)/np.float64(CHG)),str(pCHG*100)]) + '\n')
        out.write('\t'.join(["CHH",str(CHH),str(mCHH),str(np.float64(mCHH)/np.float64(CHH)),str(pCHH*100)]) + '\n')
        out.write('\t'.join(["CNN",str(CNN),str(mCNN),str(np.float64(mCNN)/np.float64(CNN)),str(pCNN*100)]) + '\n')

#calculate the genome-wide weighted methylation using non-plant specific contexts
#use: weighted_mC_nonplant("input allc", "output file", <cutoff for coverage needed to include site, default=0>)
#outputs tab-delimited file: column 1 = sequence context, column 2 = total number reads mapping, 
#column 3 = total methylated reads, column 4 = weighted methylation (decimal), column 5 = weighted methylation (as %)
def weighted_mC_nonplant(allc, mC_results, cutoff=0):    
    CG = mCG = CH = mCH = CNN = mCNN = 0    
    with open(allc) as f:
        next(f)
        for l in f:
            c = l.split('\t')
            if int(c[5]) >= int(cutoff):
                if c[3].startswith("CG"):
                    CG = CG + int(c[5])
                    mCG = mCG + int(c[4])
                elif c[3].startswith("CN") or c[3].endswith("N"):
                    CNN = CNN + int(c[5])
                    mCNN = mCNN + int(c[4])
                else:
                    CH = CH + int(c[5])
                    mCH = mCH + int(c[4])      
    with open(mC_results, "w") as out:
        out.write('\t'.join(["Context","Total","Methylated","Weighted_mC","%Weighted_mC"]) + '\n')
        out.write('\t'.join(["CG",str(CG),str(mCG),str(np.float64(mCG)/np.float64(CG)),str(pCG*100)]) + '\n')
        out.write('\t'.join(["CH",str(CH),str(mCH),str(np.float64(mCH)/np.float64(CH)),str(pCH*100)]) + '\n')
        out.write('\t'.join(["CN",str(CNN),str(mCNN),str(np.float64(mCNN)/np.float64(CNN)),str(pCNN*100)]) + '\n')

#calculate the genome-wide percent methylated sites using plant specific contexts. 
#This relies on the binomial test conducted by methylpy to determine if a site is methylated or unmethylated
#use: weighted_mC_plant("input allc", "output file", <cutoff for coverage needed to include site, default=0>)
#outputs tab-delimited file: column 1 = sequence context, column 2 = total number reads mapping, 
#column 3 = total methylated reads, column 4 = percent methylation (decimal), column 5 = percent methylation (as %)
def percent_mC_plant(allc, pC_results, cutoff=0):
    CG = mCG = CHG = mCHG = CHH = mCHH = CNN = mCNN = 0
    with open(allc) as f:
        next(f)
        for l in f:
            c = l.split('\t')
            if int(c[5]) >= int(cutoff):
                if c[3].startswith("CG"):
                    CG = CG + 1
                    mCG = mCG + int(c[6])
                elif c[3].endswith("G"):
                    CHG = CHG + 1
                    mCHG = mCHG + int(c[6])
                elif c[3].startswith("CN") or c[3].endswith("N"):
                    CNN = CNN + 1
                    mCNN = mCNN + int(c[6])
                else:
                    CHH = CHH + 1
                    mCHH = mCHH + int(c[6])       
    with open(pC_results, "w") as out:
        out.write('\t'.join(["Context","Total","Methylated","Percent_mC","%Percent_mC"]) + '\n')
        out.write('\t'.join(["CG",str(CG),str(mCG),str(np.float64(mCG)/np.float64(CG)),str(pCG*100)]) + '\n')
        out.write('\t'.join(["CHG",str(CHG),str(mCHG),str(np.float64(mCHG)/np.float64(CHG)),str(pCHG*100)]) + '\n')
        out.write('\t'.join(["CHH",str(CHH),str(mCHH),str(np.float64(mCHH)/np.float64(CHH)),str(pCHH*100)]) + '\n')
        out.write('\t'.join(["CNN",str(CNN),str(mCNN),str(np.float64(mCNN)/np.float64(CNN)),str(pCNN*100)]) + '\n')

#calculate the genome-wide percent methylated sites using nonplant specific contexts. 
#This relies on the binomial test conducted by methylpy to determine if a site is methylated or unmethylated
#use: weighted_mC_plant("input allc", "output file", <cutoff for coverage needed to include site, default=0>)
#outputs tab-delimited file: column 1 = sequence context, column 2 = total number reads mapping, 
#column 3 = total methylated reads, column 4 = percent methylation (decimal), column 5 = percent methylation (as %)        
def percent_mC_nonplant(allc, pC_results, cutoff=0):
    CG = mCG = CH = mCH = CNN = mCNN = 0    
    with open(allc) as f:
        next(f)
        for l in f:
            c = l.split('\t')
            if int(c[5]) >= int(cutoff):
                if c[3].startswith("CG"):
                    CG = CG + 1
                    mCG = mCG + int(c[6])
                elif c[3].startswith("CN") or c[3].endswith("N"):
                    CNN = CNN + 1
                    mCNN = mCNN + int(c[6])
                else:
                    CH = CH + 1
                    mCH = mCH + int(c[6])         
    with open(mC_results, "w") as out:
        out.write('\t'.join(["Context","Total","Methylated","Percent_mC","%Percent_mC"]) + '\n')
        out.write('\t'.join(["CG",str(CG),str(mCG),str(np.float64(mCG)/np.float64(CG)),str(pCG*100)]) + '\n')
        out.write('\t'.join(["CH",str(CH),str(mCH),str(np.float64(mCH)/np.float64(CH)),str(pCH*100)]) + '\n')
        out.write('\t'.join(["CN",str(CNN),str(mCNN),str(np.float64(mCNN)/np.float64(CNN)),str(pCNN*100)]) + '\n')
        
#calculate methylation levels for a region using plant specific contexts
#use: calc_feature_mC_plant( <input file>, <output file>, <cutoff for coverage needed to include site, default=0>)
#takes as input, output from map_allc2feature()
#outputs tab-delimited file: column 1 = feature name, column 2 = total CG sites, column 3 = methylated CG sites, 
#column 4 = percent methylated CGs, column 5 = total CG reads, column 6 = methylated CG reads, 
#column 7 = weighted CG methylation, column 8 = total CHG sites, column 9 = methylated CHG sites, 
#column 10 = percent methylated CHGs, column 11 = total CHG reads, column 12 = methylated CHG reads, 
#column 13 = weighted CHG methylation, column 14 = total CHH sites, column 15 = methylated CHH sites, 
#column 16 = percent methylated CHHs, column 17 = total CHH reads, column 18 = methylated CHH reads, 
#column 19 = weighted CHH methylation 
def calc_feature_mC_plant( allc2feature, feature_mC, cutoff=0):
    CG = mCG = CGr = mCGr = CHG = mCHG = CHGr = mCHGr = CHH = mCHH = CHHr = mCHHr = 0
    name = "none"
    out = open(feature_mC, "w")
    for c in allc2feature.itertuples():
        if name == "none":
            name = c[1]
            if int(c[4]) >= int(cutoff):
                if c[2].startswith("CN") or c[2].endswith("N"):
                    continue
                elif c[2].startswith("CG"):   
                    CG = CG + 1
                    mCG = mCG + int(c[5])
                    CGr = CGr + int(c[4])
                    mCGr = mCGr + int(c[3])
                elif c[2].endswith("G"):
                    CHG = CHG + 1
                    mCHG = mCHG + int(c[5])
                    CHGr = CHGr + int(c[4])
                    mCHGr = mCHGr + int(c[3])
                else:
                    CHH = CHH + 1
                    mCHH = mCHH + int(c[5])
                    CHHr = CHHr + int(c[4])
                    mCHHr = mCHHr + int(c[3])
            out.write('\t'.join(["Feature","Num_CG","mCG","Percent_mCG","Num_CG_reads","mCG_reads","Weighted_mCG","Num_CHG","mCHG","Percent_mCHG","Num_CHG_reads","mCHG_reads","Weighted_mCHG","Num_CHH","mCHH","Percent_mCHH","Num_CHH_reads","mCHH_reads","Weighted_mCHH"]) + '\n')
        elif c[1] != name:
            out.write('\t'.join([str(name),str(CG),str(mCG),str(np.float64(mCG)/np.float64(CG)),str(CGr),str(mCGr),str(np.float64(mCGr)/np.float64(CGr)),str(CHG),str(mCHG),str(np.float64(mCHG)/np.float64(CHG)),str(CHGr),str(mCHGr),str(np.float64(mCHGr)/np.float64(CHGr)),str(CHH),str(mCHH),str(np.float64(mCHH)/np.float64(CHH)),str(CHHr),str(mCHHr),str(np.float64(mCHHr)/np.float64(CHHr))]) + '\n')
            name = c[1]
            CG = mCG = CGr = mCGr = CHG = mCHG = CHGr = mCHGr = CHH = mCHH = CHHr = mCHHr = 0
            if int(c[4]) >= int(cutoff):
                if c[2].startswith("CN") or c[2].endswith("N"):
                    continue
                elif c[2].startswith("CG"):   
                    CG = CG + 1
                    mCG = mCG + int(c[5])
                    CGr = CGr + int(c[4])
                    mCGr = mCGr + int(c[3])
                elif c[2].endswith("G"):
                    CHG = CHG + 1
                    mCHG = mCHG + int(c[5])
                    CHGr = CHGr + int(c[4])
                    mCHGr = mCHGr + int(c[3])
                else:
                    CHH = CHH + 1
                    mCHH = mCHH + int(c[5])
                    CHHr = CHHr + int(c[4])
                    mCHHr = mCHHr + int(c[3])
        elif c[1] == name:
            if int(c[4]) >= int(cutoff):
                if c[2].startswith("CN") or c[2].endswith("N"):
                    continue
                elif c[2].startswith("CG"):   
                    CG = CG + 1
                    mCG = mCG + int(c[5])
                    CGr = CGr + int(c[4])
                    mCGr = mCGr + int(c[3])
                elif c[2].endswith("G"):
                    CHG = CHG + 1
                    mCHG = mCHG + int(c[5])
                    CHGr = CHGr + int(c[4])
                    mCHGr = mCHGr + int(c[3])
                else:
                    CHH = CHH + 1
                    mCHH = mCHH + int(c[5])
                    CHHr = CHHr + int(c[4])
                    mCHHr = mCHHr + int(c[3])
    out.write('\t'.join([str(name),str(CG),str(mCG),str(np.float64(mCG)/np.float64(CG)),str(CGr),str(mCGr),str(np.float64(mCGr)/np.float64(CGr)),str(CHG),str(mCHG),str(np.float64(mCHG)/np.float64(CHG)),str(CHGr),str(mCHGr),str(np.float64(mCHGr)/np.float64(CHGr)),str(CHH),str(mCHH),str(np.float64(mCHH)/np.float64(CHH)),str(CHHr),str(mCHHr),str(np.float64(mCHHr)/np.float64(CHHr))]) + '\n')
    out.close()
              
#calculate methylation levels for a region using nonplant specific contexts
#use: calc_feature_mC_nonplant( <input file>, <output file>, <cutoff for coverage needed to include site, default=0>)
#takes as input, output from map_allc2feature()
#outputs tab-delimited file: column 1 = feature name, column 2 = total CG sites, column 3 = methylated CG sites, 
#column 4 = percent methylated CGs, column 5 = total CG reads, column 6 = methylated CG reads, 
#column 7 = weighted CG methylation, column 8 = total CH sites, column 9 = methylated CH sites, 
#column 10 = percent methylated CHs, column 11 = total CH reads, column 12 = methylated CH reads, 
#column 13 = weighted CH methylation
def calc_feature_mC_nonplant( allc2feature, feature_mC, cutoff=0):
    CG = mCG = CGr = mCGr = CH = mCH = CHr = mCHr = 0
    name = "none"
    out = open(feature_mC, "w")
    for c in allc2feature.itertuples():
        if name == "none":
            name = c[1]
            if int(c[4]) >= int(cutoff):
                if c[2].startswith("CN") or c[2].endswith("N"):
                    continue
                elif c[2].startswith("CG"):   
                    CG = CG + 1
                    mCG = mCG + int(c[5])
                    CGr = CGr + int(c[4])
                    mCGr = mCGr + int(c[3])
                else:
                    CH = CH + 1
                    mCH = mCH + int(c[5])
                    CHr = CHr + int(c[4])
                    mCHr = mCHr + int(c[3])
            out.write('\t'.join(["Feature","Num_CG","mCG","Percent_mCG","Num_CG_reads","mCG_reads","Weighted_mCG","Num_CH","mCH","Percent_mCH","Num_CH_reads","mCH_reads","Weighted_mCH"]) + '\n')
        elif c[1] != name:
            out.write('\t'.join([str(name),str(CG),str(mCG),str(np.float64(mCG)/np.float64(CG)),str(CGr),str(mCGr),str(np.float64(mCGr)/np.float64(CGr)),str(CH),str(mCH),str(np.float64(mCH)/np.float64(CH)),str(CHr),str(mCHr),str(np.float64(mCHr)/np.float64(CHr))]) + '\n')
            name = c[1]
            CG = mCG = CGr = mCGr = CH = mCH = CHr = mCHHr = 0
            if int(c[4]) >= int(cutoff):
                if c[2].startswith("CN") or c[2].endswith("N"):
                    continue
                elif c[2].startswith("CG"):   
                    CG = CG + 1
                    mCG = mCG + int(c[5])
                    CGr = CGr + int(c[4])
                    mCGr = mCGr + int(c[3])
                else:
                    CH = CH + 1
                    mCH = mCH + int(c[5])
                    CHr = CHr + int(c[4])
                    mCHr = mCHr + int(c[3])
        elif c[1] == name:
            if int(c[4]) >= int(cutoff):
                if c[2].startswith("CN") or c[2].endswith("N"):
                    continue
                elif c[2].startswith("CG"):   
                    CG = CG + 1
                    mCG = mCG + int(c[5])
                    CGr = CGr + int(c[4])
                    mCGr = mCGr + int(c[3])
                else:
                    CH = CH + 1
                    mCH = mCH + int(c[5])
                    CHr = CHr + int(c[4])
                    mCHr = mCHr + int(c[3])
    out.write('\t'.join([str(name),str(CG),str(mCG),str(np.float64(mCG)/np.float64(CG)),str(CGr),str(mCGr),str(np.float64(mCGr)/np.float64(CGr)),str(CH),str(mCH),str(np.float64(mCH)/np.float64(CH)),str(CHr),str(mCHr),str(np.float64(mCHr)/np.float64(CHr))]) + '\n')
    out.close()

#Calculates average weighted methylation in plant specific contexts across bins for all input features for use in making
#metaplots
#Use: feature_metaplot_plant( "input allc file", "input bed file of features", "input required genome file for bedtools",
#"output file", ignore standedness (default = False), number of windows (default 60: 20 upstream, 20 downstream, 20 across
#feature), basepairs upstream and downstream (default is 2000 bp), cutoff on min number reads per site (default is 0)
#outputs tab-delimited file: column 1 = Bin number, column 2 = CG methylation level, column 3 = CHG methylation level, 
#column 4 = CHH methylation level
def feature_metaplot_plant(allc,features_bed,genome_file,metaplot_out,ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0):
    windows2 = int(windows/3)
    counter = 1
    bed_list = ["u","f","d"]
    mC_bed = allc2bed(allc)
    a = pd.read_table(features_bed,header=None)
    if ignoreStrand is True:
        p_bed = pbt.BedTool.from_dataframe(a)
    else:
        p_bed = pbt.BedTool.from_dataframe(a[a[5]=='+'])
        n_bed = pbt.BedTool.from_dataframe(a[a[5]=='-'])
    out = open(metaplot_out, "w")
    CG = mCG = CHG = mCHG = CHH = mCHH = 0
    out.write('\t'.join(["Bin","mCG","mCHG","mCHH"]) + '\n')
    for y in bed_list:
        if y == "u":
            if ignoreStrand is True:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=updown_stream,r=0,s=True)
            else:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=updown_stream,r=0,s=True)
                nf_bed = pbt.bedtool.BedTool.flank(n_bed,g=genome_file,l=updown_stream,r=0,s=True)
                pw_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=windows2,i="srcwinnum")
                nw_bed = pbt.bedtool.BedTool.window_maker(nf_bed,b=nf_bed,n=windows2,i="srcwinnum",reverse=True)
        elif y == "f":
            if ignoreStrand is True:
                w_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=windows2,i="srcwinnum")
            else:
                pw_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=windows2,i="srcwinnum")
        elif y == "d":
            if ignoreStrand is True:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=0,r=updown_stream,s=True)
            else:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=0,r=updown_stream,s=True)
                nf_bed = pbt.bedtool.BedTool.flank(n_bed,g=genome_file,l=0,r=updown_stream,s=True)
                pw_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=windows2,i="srcwinnum")
                nw_bed = pbt.bedtool.BedTool.window_maker(nf_bed,b=nf_bed,n=windows2,i="srcwinnum",reverse=True)
        if ignoreStrand is True:
            w_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=windows2,i="srcwinnum")
        else:
            w_bed = pw_bed.cat(nw_bed, postmerge=False)
        mapping = pbt.bedtool.BedTool.intersect(mC_bed,w_bed,wa=True,wb=True)
        m = pd.read_table(mapping.fn, header=None, usecols = [13,6,7,8])
        for x in list(range(1,windows2+1)):
            for c in m.itertuples():
                if c[4].endswith("_"+str(x)):
                    if int(c[3]) >= int(cutoff):
                        if c[1].startswith("CN") or c[1].endswith("N"):
                            continue
                        elif c[1].startswith("CG"):
                            CG = CG + int(c[3])
                            mCG = mCG + int(c[2])
                        elif c[1].endswith("G"):
                            CHG = CHG + int(c[3])
                            mCHG = mCHG + int(c[2])
                        else:
                            CHH = CHH + int(c[3])
                            mCHH = mCHH + int(c[2])
            out.write('\t'.join([str(counter),str(np.float64(mCG)/np.float64(CG)),str(np.float64(mCHG)/np.float64(CHG)),str(np.float64(mCHH)/np.float64(CHH))]) + '\n')                
            counter = counter + 1
            CG = mCG = CHG = mCHG = CHH = mCHH = 0
    out.close()
    
#Calculates average weighted methylation in nonplant specific contexts across bins for all input features for use in making
#metaplots
#Use: feature_metaplot_nonplant( "input allc file", "input bed file of features", "input required genome file for bedtools",
#"output file", ignore standedness (default = False), number of windows (default 60: 20 upstream, 20 downstream, 20 across
#feature), basepairs upstream and downstream (default is 2000 bp), cutoff on min number reads per site (default is 0)
#outputs tab-delimited file: column 1 = Bin number, column 2 = CG methylation level, column 3 = CH methylation level
def feature_metaplot_nonplant(allc,features_bed,genome_file,metaplot_out,ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0):
    windows2 = int(windows/3)
    counter = 1
    bed_list = ["u","f","d"]
    mC_bed = allc2bed(allc)
    a = pd.read_table(features_bed,header=None)
    if ignoreStrand is True:
        p_bed = pbt.BedTool.from_dataframe(a)
    else:
        p_bed = pbt.BedTool.from_dataframe(a[a[5]=='+'])
        n_bed = pbt.BedTool.from_dataframe(a[a[5]=='-'])
    out = open(metaplot_out, "w")
    CG = mCG = CH = mCH = 0
    out.write('\t'.join(["Bin","mCG","mCHH"]) + '\n')
    for y in bed_list:
        if y == "u":
            if ignoreStrand is True:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=updown_stream,r=0,s=True)
            else:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=updown_stream,r=0,s=True)
                nf_bed = pbt.bedtool.BedTool.flank(n_bed,g=genome_file,l=updown_stream,r=0,s=True)
                pw_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=windows2,i="srcwinnum")
                nw_bed = pbt.bedtool.BedTool.window_maker(nf_bed,b=nf_bed,n=windows2,i="srcwinnum",reverse=True)
        elif y == "f":
            if ignoreStrand is True:
                w_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=windows2,i="srcwinnum")
            else:
                pw_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=windows2,i="srcwinnum")
        elif y == "d":
            if ignoreStrand is True:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=0,r=updown_stream,s=True)
            else:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=0,r=updown_stream,s=True)
                nf_bed = pbt.bedtool.BedTool.flank(n_bed,g=genome_file,l=0,r=updown_stream,s=True)
                pw_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=windows2,i="srcwinnum")
                nw_bed = pbt.bedtool.BedTool.window_maker(nf_bed,b=nf_bed,n=windows2,i="srcwinnum",reverse=True)
        if ignoreStrand is True:
            w_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=windows2,i="srcwinnum")
        else:
            w_bed = pw_bed.cat(nw_bed, postmerge=False)
        mapping = pbt.bedtool.BedTool.intersect(mC_bed,w_bed,wa=True,wb=True)
        m = pd.read_table(mapping.fn, header=None, usecols = [13,6,7,8])
        for x in list(range(1,windows2+1)):
            for c in m.itertuples():
                if c[4].endswith("_"+str(x)):
                    if int(c[3]) >= int(cutoff):
                        if c[1].startswith("CN") or c[1].endswith("N"):
                            continue
                        elif c[1].startswith("CG"):
                            CG = CG + int(c[3])
                            mCG = mCG + int(c[2])
                        else:
                            CHH = CHH + int(c[3])
                            mCHH = mCHH + int(c[2])
            out.write('\t'.join([str(counter),str(np.float64(mCG)/np.float64(CG)),str(np.float64(mCH)/np.float64(CH))]) + '\n')                
            counter = counter + 1
            CG = mCG = CH = mCH = 0
    out.close()
    
#create chromosome/scaffold level metaplots for plant specific contexts
#use genome_metaplot_plant( input allc, input .genome file required by bedtools, window size (default is 100000bp), 
#stepsize for creating sliding windows (default is 100000bp, equal to window size, so no slide, must be lower than window size),
#cutoff for minimum number reads (default is 0), save window methylation data for each chrom/scaffold, will create corresponding
#tsv file for each (default is False))
#outputs pdf for each sequence in .genome file, if save_data=True, will output also text files with windows methylation data
#to customize plots, choose save_data=True and import in python, R, etc to make own plot.
def genome_metaplot_plant(allc,genome_file,windows=100000,stepsize=100000,cutoff=0,save_data=False,no_figs=False):
    mC_bed = allc2bed(allc)
    genome = pbt.BedTool(genome_file)                                                                                  
    w_bed = pbt.bedtool.BedTool.window_maker(genome,g=genome_file,w=windows,s=stepsize,i="winnum")
    mapping = pbt.bedtool.BedTool.intersect(mC_bed,w_bed,wa=True,wb=True)
    m = pd.read_table(mapping.fn, header=None, usecols = [10,13,6,7,8])
    m = m.sort_values(by = 13,ascending=True)
    genome = pd.read_table(genome_file,header=None)
    for i in genome.itertuples():
        chrn=i[1]
        chr=m[m[10]==chrn]
        df = pd.DataFrame(columns=['mCG','mCHG','mCHH'])
        name = "none"
        CG = mCG = CHG = mCHG = CHH = mCHH = 0
        for c in chr.itertuples():
            if name == "none":
                name = c[5]
                if int(c[3]) >= int(cutoff):
                    if c[1].startswith("CN") or c[1].endswith("N"):
                        continue
                    elif c[1].startswith("CG"):   
                        CG = CG + int(c[3])
                        mCG = mCG + int(c[2])                 
                    elif c[1].endswith("G"):
                        CHG = CHG + int(c[3])
                        mCHG = mCHG + int(c[2])
                    else:
                        CHH = CHH + int(c[3])
                        mCHH = mCHH + int(c[2])
            elif c[5] != name:
                df = df.append({'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH':(np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
                name = c[5]
                CG = mCG = CHG = mCHG = CHH = mCHH = 0
                if int(c[3]) >= int(cutoff):
                    if c[1].startswith("CN") or c[1].endswith("N"):
                        continue
                    elif c[1].startswith("CG"):   
                        CG = CG + int(c[3])
                        mCG = mCG + int(c[2])
                    elif c[1].endswith("G"):
                        CHG = CHG + int(c[3])
                        mCHG = mCHG + int(c[2])
                    else:
                        CHH = CHH + int(c[3])
                        mCHH = mCHH + int(c[2])
            elif c[5] == name:
                if int(c[3]) >= int(cutoff):
                    if c[1].startswith("CN") or c[1].endswith("N"):
                        continue
                    elif c[1].startswith("CG"):   
                        CG = CG + int(c[3])
                        mCG = mCG + int(c[2])                 
                    elif c[1].endswith("G"):
                        CHG = CHG + int(c[3])
                        mCHG = mCHG + int(c[2])
                    else:
                        CHH = CHH + int(c[3])
                        mCHH = mCHH + int(c[2])
        df = df.append({'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH':(np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
        if save_data is True:
            df.to_csv(str(chrn)+'.tsv', sep='\t')
        if no_figs is True:
            return df
        else:
            fig = df.plot()
            fig2 = fig.get_figure()
            fig2.savefig(str(chrn)+"_window_methylation.pdf",format="pdf")
        
#create chromosome/scaffold level metaplots for nonplant specific contexts
#use genome_metaplot_nonplant( input allc, input .genome file required by bedtools, window size (default is 100000bp), 
#stepsize for creating sliding windows (default is 100000bp, equal to window size, so no slide, must be lower than window size),
#cutoff for minimum number reads (default is 0), save window methylation data for each chrom/scaffold, will create corresponding
#tsv file for each (default is False))
#outputs pdf for each sequence in .genome file, if save_data=True, will output also text files with windows methylation data
#to customize plots, choose save_data=True and import in python, R, etc to make own plot.
def genome_metaplot_nonplant(allc,genome_file,windows=100000,stepsize=100000,cutoff=0,save_data=False):
    mC_bed = allc2bed(allc)
    genome = pbt.BedTool(genome_file)                                                                                  
    w_bed = pbt.bedtool.BedTool.window_maker(genome,g=genome_file,w=windows,s=stepsize,i="winnum")
    mapping = pbt.bedtool.BedTool.intersect(mC_bed,w_bed,wa=True,wb=True)
    m = pd.read_table(mapping.fn, header=None, usecols = [10,13,6,7,8])
    m = m.sort_values(by = 13,ascending=True)
    genome = pd.read_table(genome_file,header=None)
    for i in genome.itertuples():
        chrn=i[1]
        chr=m[m[10]==chrn]
        df = pd.DataFrame(columns=['mCG','mCH'])
        name = "none"
        CG = mCG = CH = mCH = 0
        for c in chr.itertuples():
            if name == "none":
                name = c[5]
                if int(c[3]) >= int(cutoff):
                    if c[1].startswith("CN") or c[1].endswith("N"):
                        continue
                    elif c[1].startswith("CG"):   
                        CG = CG + int(c[3])
                        mCG = mCG + int(c[2])                 
                    else:
                        CH = CH + int(c[3])
                        mCH = mCH + int(c[2])
            elif c[5] != name:
                df = df.append({'mCG':(np.float64(mCG)/np.float64(CG)), 'mCH':(np.float64(mCH)/np.float64(CH))}, ignore_index=True)
                name = c[5]
                CG = mCG = CH = mCH = 0
                if int(c[3]) >= int(cutoff):
                    if c[1].startswith("CN") or c[1].endswith("N"):
                        continue
                    elif c[1].startswith("CG"):   
                        CG = CG + int(c[3])
                        mCG = mCG + int(c[2])
                    else:
                        CH = CH + int(c[3])
                        mCH = mCH + int(c[2])
            elif c[5] == name:
                if int(c[3]) >= int(cutoff):
                    if c[1].startswith("CN") or c[1].endswith("N"):
                        continue
                    elif c[1].startswith("CG"):   
                        CG = CG + int(c[3])
                        mCG = mCG + int(c[2])                 
                    else:
                        CH = CH + int(c[3])
                        mCH = mCH + int(c[2])
        df = df.append({'mCG':(np.float64(mCG)/np.float64(CG)), 'mCH':(np.float64(mCH)/np.float64(CH))}, ignore_index=True)
        if save_data is True:
            df.to_csv(str(chrn)+'.tsv', sep='\t')
        fig = df.plot()
        fig2 = fig.get_figure()
        fig2.savefig(str(chrn)+"_window_methylation.pdf",format="pdf")
        
#get distribution of methylation levels per site for plant specific contexts
#use site_mC_level_plant( allc file, output kernal density plot for each context (default is True), save txt file with per-site
#methylation levels for each context (use for input in other applications/plotting, default is False), cutoff for minimum number
#reads mapping to site (default is 3, this is typical cutoff used in methylpy, should not set lower than cutoff in methylpy)
#outputs three kernal density plots if kernal_density = True and/or three txt files with per-site methylation level in single 
#column, for each context if save_raw = True
def site_mC_level_plant(allc,kernal_density=True,save_raw=False,cutoff=3):
    a = pd.read_table(allc)
    CG = pd.DataFrame(columns=['mCG'])
    CHG = pd.DataFrame(columns=['mCHG'])
    CHH = pd.DataFrame(columns=['mCHH'])
    for c in a.itertuples():
        if int(c[7]) == 1 and int(c[6]) >= int(cutoff):
            if c[4].startswith("CN") or c[4].endswith("N"):
                continue
            elif c[4].startswith("CG"):
                CG = CG.append({'mCG':(np.float64(c[5])/np.float64(c[6]))}, ignore_index=True)
            elif c[4].endswith("G"):
                CHG = CHG.append({'mCHG':(np.float64(c[5])/np.float64(c[6]))}, ignore_index=True)
            else:
                CHH = CHH.append({'mCHH':(np.float64(c[5])/np.float64(c[6]))}, ignore_index=True)
    if save_raw is True:
        CG.to_csv("CG_per_site_mC.tsv", sep='\t', index=False)
        CHG.to_csv("CHG_per_site_mC.tsv", sep='\t', index=False)
        CHH.to_csv("CHH_per_site_mC.tsv", sep='\t', index=False)
    if kernal_density is True:
        fig = CG.plot(kind='density')
        fig2 = fig.get_figure()
        fig2.savefig("CG_site_mC_dist.pdf",format="pdf")
        fig = CHG.plot(kind='density')
        fig2 = fig.get_figure()
        fig2.savefig("CHG_site_mC_dist.pdf",format="pdf") 
        fig = CHH.plot(kind='density')
        fig2 = fig.get_figure()
        fig2.savefig("CHH_site_mC_dist.pdf",format="pdf")
        
#get distribution of methylation levels per site for nonplant specific contexts
#use site_mC_level_plant( allc file, output kernal density plot for each context (default is True), save txt file with per-site
#methylation levels for each context (use for input in other applications/plotting, default is False), cutoff for minimum number
#reads mapping to site (default is 3, this is typical cutoff used in methylpy, should not set lower than cutoff in methylpy)
#outputs three kernal density plots if kernal_density = True and/or three txt files with per-site methylation level in single 
#column, for each context if save_raw = True
def site_mC_level_nonplant(allc,kernal_density=True,save_raw=False,cutoff=3):
    a = pd.read_table(allc)
    CG = pd.DataFrame(columns=['mCG'])
    CH = pd.DataFrame(columns=['mCH'])
    for c in a.itertuples():
        if int(c[7]) == 1 and int(c[6]) >= int(cutoff):
            if c[4].startswith("CN") or c[4].endswith("N"):
                continue
            elif c[4].startswith("CG"):
                CG = CG.append({'mCG':(np.float64(c[5])/np.float64(c[6]))}, ignore_index=True)
            else:
                CH = CH.append({'mCH':(np.float64(c[5])/np.float64(c[6]))}, ignore_index=True)
    if save_raw is True:
        CG.to_csv("CG_per_site_mC.tsv", sep='\t', index=False)
        CH.to_csv("CH_per_site_mC.tsv", sep='\t', index=False)
    if kernal_density is True:
        fig = CG.plot(kind='density')
        fig2 = fig.get_figure()
        fig2.savefig("CG_site_mC_dist.pdf",format="pdf")
        fig = CH.plot(kind='density')
        fig2 = fig.get_figure()
        fig2.savefig("CH_site_mC_dist.pdf",format="pdf")

#get fraction methylated sites per each trinucleotide context for both allc and fasta files
#TO DO: add in plotting functions
def subcontext_mC_hist(allc,fasta,output=(),cutoff=3):
    allc_df = count_subcontext_allc(allc,cutoff=cutoff)
    fasta_df = count_subcontext_fasta(fasta)
    df = pd.concat([fasta_df,allc_df])
    df.loc['fraction_fasta'] = df.loc['mC'] / df.loc['total_fasta']
    if output:
        df.to_csv(output, sep='\t', index=False)
        
#get fraction methylated sites per each trinucleotide context for allc files only
def count_subcontext_allc(allc,output=(),cutoff=3):
    a = pd.read_table(allc)
    df = pd.DataFrame(columns=["CAA","CAT","CAC","CAG","CTA","CTT","CTC","CTG","CCA","CCT","CCC","CCG","CGA","CGT","CGC","CGG"])
    for i in ["CAA","CAT","CAC","CAG","CTA","CTT","CTC","CTG","CCA","CCT","CCC","CCG","CGA","CGT","CGC","CGG"]:
        i_table = a[(a['mc_class']==i) & (a['total']>=cutoff)]
        total = len(i_table.index)
        mC = len(i_table[i_table['methylated']==1].index)
        df.set_value('mC', i, mC)
        df.set_value('total_allc', i, total) 
        df.set_value('fraction_allc', i, np.float64(mC)/np.float64(total)) 
    if output:
        df.to_csv(output, sep='\t', index=False)
    else:
        return df   
    
#get number sites per each trinucleotide context for fasta files
def count_subcontext_fasta(fasta,output=(),per_seq=False):
    CAA = CAT = CAC = CAG = CTA = CTT = CTC = CTG = CCA = CCT = CCC = CCG = CGA = CGT = CGC = CGG = 0
    df = pd.DataFrame(columns=["Seq_name","CAA","CAT","CAC","CAG","CTA","CTT","CTC","CTG","CCA","CCT","CCC","CCG","CGA","CGT","CGC","CGG"])
    for sequence in SeqIO.parse(fasta, "fasta"):
        if sequence.name not in ["C","L","M","chrC","chrL","chrM","CHLOROPLAST","MITOCHONDRIA"]:
            CAA_count = sequence.seq.count('CAA') + sequence.seq.reverse_complement().count('CAA')
            CAA = CAA + CAA_count
            CAT_count = sequence.seq.count('CAT') + sequence.seq.reverse_complement().count('CAT')
            CAT = CAT + CAT_count
            CAC_count = sequence.seq.count('CAC') + sequence.seq.reverse_complement().count('CAC')
            CAC = CAC + CAC_count
            CAG_count = sequence.seq.count('CAG') + sequence.seq.reverse_complement().count('CAG')
            CAG = CAG + CAG_count
            CTA_count = sequence.seq.count('CTA') + sequence.seq.reverse_complement().count('CTA')
            CTA = CTA + CTA_count
            CTT_count = sequence.seq.count('CTT') + sequence.seq.reverse_complement().count('CTT')
            CTT = CTT + CTT_count
            CTC_count = sequence.seq.count('CTC') + sequence.seq.reverse_complement().count('CTC')
            CTC = CTC + CTC_count
            CTG_count = sequence.seq.count('CTG') + sequence.seq.reverse_complement().count('CTG')
            CTG = CTG + CTG_count
            CCA_count = sequence.seq.count('CCA') + sequence.seq.reverse_complement().count('CCA')
            CCA = CCA + CCA_count
            CCT_count = sequence.seq.count('CCT') + sequence.seq.reverse_complement().count('CCT')
            CCT = CCT + CCT_count
            CCC_count = sequence.seq.count('CCC') + sequence.seq.reverse_complement().count('CCC')
            CCC = CCC + CCC_count
            CCG_count = sequence.seq.count('CCG') + sequence.seq.reverse_complement().count('CCG')
            CCG = CCG + CCG_count
            CGA_count = sequence.seq.count('CGA') + sequence.seq.reverse_complement().count('CGA')
            CGA = CGA + CGA_count
            CGT_count = sequence.seq.count('CGT') + sequence.seq.reverse_complement().count('CGT')
            CGT = CGT + CGT_count
            CGC_count = sequence.seq.count('CGC') + sequence.seq.reverse_complement().count('CGC')
            CGC = CGC + CGC_count
            CGG_count = sequence.seq.count('CGG') + sequence.seq.reverse_complement().count('CGG')
            CGG = CGG + CGG_count
            if per_seq is True:
                df = df.append({'Seq_name':sequence.name,'CAA':CAA_count,'CAT':CAT_count,'CAC':CAC_count,'CAG':CAG_count,'CTA':CTA_count,'CTT':CTT_count,'CTC':CTC_count,'CTG':CTG_count,'CCA':CCA_count,'CCT':CCT_count,'CCC':CCC_count,'CCG':CCG_count,'CGA':CGA_count,'CGT':CGT_count,'CGC':CGC_count,'CGG':CGG_count}, ignore_index=True)
    if per_seq is False:
        df = df.append({'Seq_name':"All",'CAA':CAA,'CAT':CAT,'CAC':CAC,'CAG':CAG,'CTA':CTA,'CTT':CTT,'CTC':CTC,'CTG':CTG,'CCA':CCA,'CCT':CCT,'CCC':CCC,'CCG':CCG,'CGA':CGA,'CGT':CGT,'CGC':CGC,'CGG':CGG}, ignore_index=True)
    if output:
        df.to_csv(output, sep='\t', index=False)
    else:
        df = df.rename(index={0: 'total_fasta'})
        df = df.drop('Seq_name', 1)
        return df

#filter DMR output from methylpy
#use filter_dmr( collapsed DMR file from methylpy, output file for filtered data, minimum number DMS (default = 5)
#minimum difference in methylation between highest and lowest methylated sample (default = 0.1 or 10% difference)
def filter_dmr(dmr_file,output=(),min_dms=5,min_mC_diff=0.1):
    a=pd.read_table(dmr_file)
    list=[]
    for c in a.itertuples():
        if c[4] >= min_dms:
            if max(c[7:])-min(c[7:]) >= min_mC_diff:
                list.append(c[0])
    a=a.ix[list]
    if output:
        a.to_csv(output, sep='\t', index=False)
    else:
        return a
    
#This script creates a 2x2 table showing methylation status when two cytosines next to each other.
def mCmC_association(allc,cutoff=3,context=["C"]):
    a = filter_context(allc,context)
    t = np.zeros([2, 2])
    i=0
    l=len(a)
    printProgress(i, l, prefix = 'Progress:', suffix = 'Complete', barLength = 50)
    for c in a.itertuples():
        if "CC" in c[4] and c[6] >= cutoff:
            chr2 = c[1]
            if c[3] == "+":
                ind = c[0] + 1
                pos2 = c[2] + 1
            elif c[3] == "-":
                ind = c[0] - 1
                pos2 = c[2] - 1
            n = a.loc[ind]
            if not n.empty:
                if n.chr == chr2 and n.pos == pos2 and int(n.total) >= cutoff:
                    if c[7] == 1 and int(n.methylated) == 1:
                        t[0,0] = t[0,0] + 1
                    elif c[7] == 0 and int(n.methylated) == 1:
                        t[1,0] = t[1,0] + 1
                    elif c[7] == 1 and int(n.methylated) == 0: 
                        t[0,1] = t[0,1] + 1
                    elif c[7] == 0 and int(n.methylated) == 0: 
                        t[1,1] = t[1,1] + 1  
        i += 1
        printProgress(i, l, prefix = 'Progress:', suffix = 'Complete', barLength = 50)
    return t

#Get total number, number of methylated, and percent methylated sites for expanded (more or less than 3 bases) seq contexts
def mC_multibase_context(allc,fasta,genome_file,output=(),up=1,down=2,cutoff=3,context=["C"]):
    a = allc2bed(allc,context,bed=False)
    a = a[a.total >= cutoff]
    p_bed = pbt.bedtool.BedTool.slop(pbt.BedTool.from_dataframe(a[a.strand == '+']),g=genome_file,l=up+1,r=down,s=True)
    n_bed = pbt.bedtool.BedTool.slop(pbt.BedTool.from_dataframe(a[a.strand == '-']),g=genome_file,l=up,r=down+1,s=True)
    m = pd.read_table(p_bed.cat(n_bed,postmerge=False).nucleotide_content(fi=fasta,s=True,seq=True).fn)
    C = Counter(m["20_seq"])
    mC = Counter(m[m["10_usercol"] == 1]["20_seq"])
    df = pd.DataFrame(columns=["Seq","Total","Methylated","Percent_Methylated"])
    for i in C:
        df = df.append({'Seq':i,'Total':C[i],'Methylated':mC[i],'Percent_Methylated':np.float64(mC[i])/np.float64(C[i])},ignore_index=True)
    if output:
        df.to_csv(output, sep='\t', index=False)
    return df

def coverage_stats(allc,fasta,output=(),plot=False,context=['C']):
    a=filter_context(allc,context)
    a=pd.DataFrame.from_dict(Counter(a.total),orient='index')
    a.columns=['percent']
    C_count=0
    for sequence in SeqIO.parse(fasta, "fasta"):
        if sequence.name not in ["C","L","M","chrC","chrL","chrM","CHLOROPLAST","MITOCHONDRIA"]:
            if "C" in context:
                C_count = sequence.seq.count('C') + sequence.seq.reverse_complement().count('C') + C_count
            else:
                for i in expand_nucleotide_code(context):
                    C_count = sequence.seq.count(i) + sequence.seq.reverse_complement().count(i) + C_count
    a=a[[0]]/C_count
    a.columns=['percent']
    if output:
        a.to_csv(output, sep='\t', index=False)
    if plot:
        density = stats.gaussian_kde(list(a.percent))
        dist_space = linspace( min(a.percent), max(a.percent), 100 )
        #density.covariance_factor = lambda : .25
        #density._compute_covariance()
        plt.plot(dist_space,density(dist_space))
        plt.savefig("cytosine_coverage_dist.pdf",format='pdf')
    return a
        
if __name__ == "__main__":
    pass