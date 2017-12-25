import os
import sys
import pybedtools as pbt
import pandas as pd
import itertools
import numpy as np
from scipy import stats
from scipy.stats.stats import pearsonr
from collections import Counter
from Bio import SeqIO

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
	if open(allc).readline().rstrip() == 'chr\tpos\tstrand\tmc_class\tmc_count\ttotal\tmethylated':
		a = pd.read_table(allc,dtype={'chr':str,'pos':int,'strand':str,'mc_class':str,
			'mc_count':int,'total':int,'methylated':int})
	else:
		a = pd.read_table(allc,names=['chr','pos','strand','mc_class','mc_count',
			'total','methylated'],dtype={'chr':str,'pos':int,'strand':str,
			'mc_class':str,'mc_count':int,'total':int,'methylated':int})
	a = a[a.mc_class.isin(expand_nucleotide_code(context))]
	return a

def allc2bed(allc,context=["C"],bed=True):
    a = filter_context(allc,context)
    a['pos2'] = a.pos
    a['name'] = a.index
    a['score'] = "."
    a = a[['chr','pos','pos2','name','score','strand','mc_class','mc_count','total','methylated']]
    if bed is True:
        a = pbt.BedTool.from_dataframe(a)
    return a

#simple function for filtering gff files based on feature (gene, exon, mRNA, etc)
def feat_filter(x,feature):
    if feature:
        return x[2] == feature
    else:
        return x

#simple function for filtering gff files based on strand
def strand_filter(x,strand):
    return x.strand == strand

#simple function for filtering gff files based on chromosome
def chr_filter(x,chr):
    return x.chrom not in chr

#make an intron gff file
def intron_gff(features,output=(),remove_db=True):
    db=gffutils.create_db(features, 'gff.db')
    gff_out=gffutils.gffwriter.GFFWriter(output, in_place=True)
    name="none"
    count=1
    for rec in db.create_introns():
    	if rec['Parent'] == name:
        	count=count+1
        	rec.attributes['ID'] = [str(rec['Parent'])[2:-2] + '.intron' + str(count)]
       		gff_out.write_rec(rec)
    	else:
        	name = rec['Parent']
        	count = 1
        	rec.attributes['ID'] = [str(rec['Parent'])[2:-2] + '.intron' + str(count)]
        	gff_out.write_rec(rec)
    gff_out.close()
    if remove_db:
        os.remove('gff.db')

#generic function for calculating methylation levels in windows
def window_methylation_levels(m,cutoff=0,nuc_bed=(),output_mC_counts=False):
    if nuc_bed:
        if output_mC_counts:
            a = pd.DataFrame(columns=['window','mCG_reads','CG_reads','mCG',
                                      'mCHG_reads','CHG_reads','mCHG','mCHH_reads',
                                       'CHH_reads','mCHH','coverage'])
        else:
             a = pd.DataFrame(columns=['window','mCG','mCHG','mCHH','coverage'])
    else:
        if output_mC_counts:
            a = pd.DataFrame(columns=['window','mCG_reads','CG_reads','mCG',
                                      'mCHG_reads','CHG_reads','mCHG','mCHH_reads',
                                       'CHH_reads','mCHH'])
        else:
             a = pd.DataFrame(columns=['window','mCG','mCHG','mCHH'])
    name = "none"
    C = CG = mCG = CHG = mCHG = CHH = mCHH = GC = 0
    if nuc_bed:
        nuc = pd.read_table(nuc_bed.fn, usecols = [3,7,8])
        m = pd.merge(m,nuc,left_on=13,right_on='4_usercol')
    for c in m.itertuples():
      if name == "none":
        name = c[5]
        if nuc_bed:
            GC = int(c[7]) + int(c[8])
        if int(c[3]) >= int(cutoff):
          C = C + 1
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
        if nuc_bed:
            if output_mC_counts:
                a = a.append({'window':str(name), 'mCG_reads':int(mCG), 'CG_reads':int(CG), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG_reads':int(mCHG), 'CHG_reads':int(CHG), 'mCHG':(np.float64(mCHG)/np.float64(CHG)),  'mCHH_reads':int(mCHH), 'CHH_reads':int(CHH), 'mCHH':(np.float64(mCHH)/np.float64(CHH)), 'coverage':(np.float64(C)/np.float(GC))}, ignore_index=True)
            else:
                a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH':(np.float64(mCHH)/np.float64(CHH)), 'coverage':(np.float64(C)/np.float(GC))}, ignore_index=True)
        else:
            if output_mC_counts:
                a = a.append({'window':str(name), 'mCG_reads':int(mCG), 'CG_reads':int(CG), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG_reads':int(mCHG), 'CHG_reads':int(CHG), 'mCHG':(np.float64(mCHG)/np.float64(CHG)),  'mCHH_reads':int(mCHH), 'CHH_reads':int(CHH), 'mCHH':(np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
            else:
                a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH':(np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
        name = c[5]
        if nuc_bed:
            GC = int(c[7]) + int(c[8])
        C = CG = mCG = CHG = mCHG = CHH = mCHH = 0
        if int(c[3]) >= int(cutoff):
            C = C + 1
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
        C = C + 1
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
    if nuc_bed:
        if output_mC_counts:
            a = a.append({'window':str(name), 'mCG_reads':int(mCG), 'CG_reads':int(CG), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG_reads':int(mCHG), 'CHG_reads':int(CHG), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH_reads':int(mCHH), 'CHH_reads':int(CHH), 'mCHH':(np.float64(mCHH)/np.float64(CHH)), 'coverage':(np.float64(C)/np.float(GC))}, ignore_index=True)
        else:
            a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH':(np.float64(mCHH)/np.float64(CHH)), 'coverage':(np.float64(C)/np.float(GC))}, ignore_index=True)
    else:
        if output_mC_counts:
            a = a.append({'window':str(name), 'mCG_reads':int(mCG), 'CG_reads':int(CG), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG_reads':int(mCHG), 'CHG_reads':int(CHG), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH_reads':int(mCHH), 'CHH_reads':int(CHH), 'mCHH':(np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
        else:
            a = a.append({'window':str(name), 'mCG':(np.float64(mCG)/np.float64(CG)), 'mCHG':(np.float64(mCHG)/np.float64(CHG)), 'mCHH':(np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
    return a

#get correlations of mC with genes
def mC_distribution(genome_file,allc,gff_file,fasta,windows=100000,stepsize=50000,cutoff=0,filter=0.5,feature='gene',output_path=(),save_windows=False,return_windows=False):
    w_bed = pbt.bedtool.BedTool.window_maker(pbt.BedTool(genome_file),g=genome_file,w=windows,s=stepsize,i='srcwinnum')
    nuc = pbt.BedTool.nucleotide_content(w_bed,fasta)
    mC_bed = allc2bed(allc)
    gene_gff = pbt.BedTool(gff_file).filter(feat_filter,feature=feature)
    gene_mapping = pbt.bedtool.BedTool.intersect(w_bed,gene_gff,wa=True)
    g = pd.DataFrame.from_dict(Counter([(f[3]) for f in gene_mapping]), orient='index').reset_index()
    g.columns=['window','genes']
    allc_mapping = pbt.bedtool.BedTool.intersect(mC_bed,w_bed,wa=True,wb=True)
    m = pd.read_table(allc_mapping.fn,header=None,usecols=[10,13,6,7,8])
    m = m.sort_values(by = 13,ascending=True)
    a = window_methylation_levels(m,cutoff,nuc)
    df = pd.merge(a,g,on='window')
    df2 = df[df['coverage'] >= filter]
    correlations = pd.DataFrame(columns=['Context','r','p-value'])
    corCG = pearsonr(df2['mCG'],df2['genes'])
    correlations = correlations.append({'Context': 'mCG', 'r': corCG[0], 'p-value': corCG[1]}, ignore_index=True)
    corCHG = pearsonr(df2['mCHG'],df2['genes'])
    correlations = correlations.append({'Context': 'mCHG', 'r': corCHG[0], 'p-value': corCHG[1]}, ignore_index=True)
    corCHH = pearsonr(df2['mCHH'],df2['genes'])
    correlations = correlations.append({'Context': 'mCHH', 'r': corCHH[0], 'p-value': corCHH[1]}, ignore_index=True)
    if output_path:
        correlations.to_csv(output_path + "correlations.tsv", sep='\t', index=False)
        df.to_csv(output_path + "unfiltered_genome_windows_data.tsv", sep='\t', index=False)
        df2.to_csv(output_path + "filtered_genome_windows_data.tsv", sep='\t', index=False)
    else:
        return correlations

#get total weighted methylation
def weighted_mC(allc, output=(), cutoff=0):
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
    a = pd.DataFrame(columns=['Context','Total','Methylated','Weighted_mC'])
    a = a.append({'Context': 'mCG', 'Total': CG, 'Methylated': mCG, 'Weighted_mC': (np.float64(mCG)/np.float64(CG))}, ignore_index=True)
    a = a.append({'Context': 'mCHG', 'Total': CHG, 'Methylated': mCHG, 'Weighted_mC': (np.float64(mCHG)/np.float64(CHG))}, ignore_index=True)
    a = a.append({'Context': 'mCHH', 'Total': CHH, 'Methylated': mCHH, 'Weighted_mC': (np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
    a = a.append({'Context': 'mCNN', 'Total': CNN, 'Methylated': mCNN, 'Weighted_mC': (np.float64(mCNN)/np.float64(CNN))}, ignore_index=True)
    if output:
        a.to_csv(output, sep='\t', index=False)
    else:
        return a

#plot methylation level for features
def feature_metaplot(allc,features,genome_file,output=(),ignoreStrand=False,windows=60,updown_stream=2000,cutoff=0,filter_features=(),filter_chr=[]):
    counter = 1
    mC_bed = allc2bed(allc)
    a = pbt.BedTool(features)
    if ignoreStrand:
        p_bed = a.filter(feat_filter,filter_features).filter(chr_filter,filter_chr).saveas('p_tmp')
    else:
        p_bed = a.filter(strand_filter,strand='+').filter(feat_filter,filter_features).filter(chr_filter,filter_chr).saveas('p_tmp')
        n_bed = a.filter(strand_filter,strand='-').filter(feat_filter,filter_features).filter(chr_filter,filter_chr).saveas('n_tmp')
    CG = mCG = CHG = mCHG = CHH = mCHH = 0
    metaplot = pd.DataFrame(columns=['Bin','mCG','mCHG','mCHH'])
    for y in ['u','f','d']:
        if y == 'u':
            if ignoreStrand:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=updown_stream,r=0,s=True)
            else:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=updown_stream,r=0,s=True)
                nf_bed = pbt.bedtool.BedTool.flank(n_bed,g=genome_file,l=updown_stream,r=0,s=True)
                pw_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=int(windows/3),i='srcwinnum')
                nw_bed = pbt.bedtool.BedTool.window_maker(nf_bed,b=nf_bed,n=int(windows/3),i='srcwinnum',reverse=True)
        elif y == 'f':
            if ignoreStrand:
                w_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=int(windows/3),i='srcwinnum')
            else:
                pw_bed = pbt.bedtool.BedTool.window_maker(p_bed,b=p_bed,n=int(windows/3),i='srcwinnum')
                nw_bed = pbt.bedtool.BedTool.window_maker(n_bed,b=n_bed,n=int(windows/3),i='srcwinnum',reverse=True)
        elif y == 'd':
            if ignoreStrand:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=0,r=updown_stream,s=True)
            else:
                pf_bed = pbt.bedtool.BedTool.flank(p_bed,g=genome_file,l=0,r=updown_stream,s=True)
                nf_bed = pbt.bedtool.BedTool.flank(n_bed,g=genome_file,l=0,r=updown_stream,s=True)
                pw_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=int(windows/3),i='srcwinnum')
                nw_bed = pbt.bedtool.BedTool.window_maker(nf_bed,b=nf_bed,n=int(windows/3),i='srcwinnum',reverse=True)
        if ignoreStrand:
            w_bed = pbt.bedtool.BedTool.window_maker(pf_bed,b=pf_bed,n=int(windows/3),i='srcwinnum')
        else:
            w_bed = pw_bed.cat(nw_bed, postmerge=False)
        mapping = pbt.bedtool.BedTool.intersect(mC_bed,w_bed,wa=True,wb=True)
        m = pd.read_table(mapping.fn, header=None, usecols = [13,6,7,8])
        for x in list(range(1,int(windows/3)+1)):
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
            metaplot = metaplot.append({'Bin': counter,'mCG': (np.float64(mCG)/np.float64(CG)),
                                        'mCHG': (np.float64(mCHG)/np.float64(CHG)),
                                        'mCHH': (np.float64(mCHH)/np.float64(CHH))}, ignore_index=True)
            counter = counter + 1
            CG = mCG = CHG = mCHG = CHH = mCHH = 0
    os.remove('p_tmp')
    os.remove('n_tmp')
    if output:
        metaplot.to_csv(output, sep='\t', index=False)
    else:
        return metaplot

#map methylation to features
def map2features(allc,features,genome_file,updown_stream=2000,first_feature=(),second_feature=(),filter_chr=[]):
    bed = pbt.BedTool(features).filter(feat_filter,first_feature).filter(chr_filter,filter_chr)
    flank_bed = pbt.bedtool.BedTool.flank(bed,g=genome_file,l=updown_stream,r=updown_stream,s=True).saveas('f_tmp')
    cds_bed = pbt.BedTool(features).filter(feat_filter,second_feature).filter(chr_filter,filter_chr).saveas('c_tmp')
    bed = cds_bed.cat(flank_bed, postmerge=False)
    mC_bed = allc2bed(allc)
    mapping = pbt.bedtool.BedTool.intersect(mC_bed,bed,wa=True)
    m = pd.read_table(mapping.fn, header=None, usecols = [0,1,5,6,7,8,9])
    m.columns = ['chr','pos','strand','mc_class','mc_count','total','methylated']
    m = m.drop_duplicates()
    m.to_csv('CDS_allc.tmp', sep='\t', index=False)
    os.remove('f_tmp')
    os.remove('c_tmp')

#plot methylation levels for genes
def gene_metaplot(allc,features,genome_file,output=(),ignoreStrand=False,windows=60,
                  updown_stream=2000,cutoff=0,first_feature=(),second_feature=(),
                  filter_chr=[],remove_tmp=True):
    map2features(allc,features,genome_file,updown_stream,first_feature,second_feature,filter_chr)
    feature_metaplot('CDS_allc.tmp',features,genome_file,output,ignoreStrand,
                     windows,updown_stream,cutoff,first_feature,filter_chr)
    if remove_tmp:
        os.remove('CDS_allc.tmp')

#output per-site methylation levels for mCs in each specified context
def per_site_mC(allc,output_path,context=['CG','CHG','CHH']):
    for i in context:
        a = filter_context(allc,[i])
        a = a[a['methylated'] == 1]
        a['mc_level'] = a['mc_count']/a['total']
        a.to_csv(output_path+i+'_site_methylation.txt', sep='\t', index=False)

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

#Count up subcontext methylation
def count_subcontext_allc(allc,context=['CG','CHG','CHH'],output=(),cutoff=3,filter_chr=[]):
    a = pd.read_table(allc)
    df = pd.DataFrame(columns=['context','total_passed','methylated','weighted_mC'])
    for c in context:
        tp = mC = t = m = 0
        for i in expand_nucleotide_code([c]):
            i_table = a[(a['mc_class']==i) & (~a['chr'].isin(filter_chr))]
            t = t + i_table['total'].values.sum()
            m = m + i_table['mc_count'].values.sum()
            i_table = i_table[i_table['total']>=cutoff]
            tp = tp + len(i_table.index)
            mC = mC + len(i_table[i_table['methylated']==1].index)
        df = df.append({'context': c, 'total_passed': tp, 'methylated': mC,
                        "weighted_mC": np.float64(m)/np.float64(t)}, ignore_index=True)
    if output:
        df.to_csv(output, sep='\t', index=False)
    else:
        return df

#Get subcontexts of genome and methylation
def subcontext_methylation(allc,fasta,context=['CG','CHG','CHH'],output=(),cutoff=3,filter_chr=[]):
    a = count_subcontext_fasta(fasta,context,output=(),filter_chr=filter_chr)
    b = count_subcontext_allc(allc,context,output=(),cutoff=3,filter_chr=filter_chr)
    df = pd.merge(a,b,on='context')
    if output:
        df.to_csv(output, sep='\t', index=False)
    else:
        return df

#Get feature methylation data
def feature_mC_levels(allc,features,output=(),cutoff=0,filter_features=(),filter_chr=[]):
    bed = pbt.BedTool(features).filter(feat_filter,filter_features).filter(chr_filter,filter_chr)
    mC_bed = allc2bed(allc)
    mapping = pbt.bedtool.BedTool.intersect(mC_bed,bed,wa=True,wb=True)
    m = pd.read_table(mapping.fn, header=None, usecols = [18,6,7,8,9])
    m = m.sort_values(by = 18,ascending=True)
    a = pd.DataFrame(columns=['Gene','CG_sites','mCG_sites','CG_reads','mCG_reads',
                              'CG_methylation_level','CHG_sites','mCHG_sites','CHG_reads',
                              'mCHG_reads','CHG_methylation_level','CHH_sites','mCHH_sites',
                              'CHH_reads','mCHH_reads','CHH_methylation_level','C_sites',
                              'mC_sites','C_reads','mC_reads','C_methylation_level'])
    name = "none"
    rCG = mrCG = CG = mCG = rCHG = mrCHG = CHG = mCHG = rCHH = mrCHH = CHH = mCHH = 0
    for c in m.itertuples():
        if name == "none":
            name = c[5]
            if int(c[3]) >= int(cutoff):
                if c[1].startswith("CN"):
                    continue
                elif c[1].startswith("CG"):
                    rCG = rCG + int(c[3])
                    mrCG = mrCG + int(c[2])
                    CG = CG + 1
                    mCG = mCG + int(c[4])
                elif c[1].endswith("G"):
                    rCHG = rCHG + int(c[3])
                    mrCHG = mrCHG + int(c[2])
                    CHG = CHG + 1
                    mCHG = mCHG + int(c[4])
                else:
                    rCHH = rCHH + int(c[3])
                    mrCHH = mrCHH + int(c[2])
                    CHH = CHH + 1
                    mCHH = mCHH + int(c[4])
        elif c[5] != name:
            a = a.append({'Gene':name, 'CG_sites':CG, 'mCG_sites':mCG,
                          'CG_reads':rCG, 'mCG_reads':mrCG,
                          'CG_methylation_level':str(np.float64(mrCG)/np.float64(rCG)),
                          'CHG_sites':CHG, 'mCHG_sites':mCHG, 'CHG_reads':rCHG,
                          'mCHG_reads':mrCHG, 'CHG_methylation_level':str(np.float64(mrCHG)/np.float64(rCHG)),
                          'CHH_sites':CHH, 'mCHH_sites':mCHH, 'CHH_reads':rCHH, 'mCHH_reads':mrCHH,
                          'CHH_methylation_level':str(np.float64(mrCHH)/np.float64(rCHH)),
                          'C_sites':(CG+CHG+CHH), 'mC_sites':(mCG+mCHG+mCHH),
                          'C_reads':(rCG+rCHG+rCHH), 'mC_reads':(mrCG+mrCHG+mrCHH),
                          'C_methylation_level':str(np.float64(mrCG+mrCHG+mrCHH)/np.float64(rCG+rCHG+rCHH))
                          },ignore_index=True)
            name = c[5]
            rCG = mrCG = CG = mCG = rCHG = mrCHG = CHG = mCHG = rCHH = mrCHH = CHH = mCHH = 0
            if int(c[3]) >= int(cutoff):
                if c[1].startswith("CN"):
                    continue
                elif c[1].startswith("CG"):
                    rCG = rCG + int(c[3])
                    mrCG = mrCG + int(c[2])
                    CG = CG + 1
                    mCG = mCG + int(c[4])
                elif c[1].endswith("G"):
                    rCHG = rCHG + int(c[3])
                    mrCHG = mrCHG + int(c[2])
                    CHG = CHG + 1
                    mCHG = mCHG + int(c[4])
                else:
                    rCHH = rCHH + int(c[3])
                    mrCHH = mrCHH + int(c[2])
                    CHH = CHH + 1
                    mCHH = mCHH + int(c[4])
        elif c[5] == name:
            if int(c[3]) >= int(cutoff):
                if c[1].startswith("CN"):
                    continue
                elif c[1].startswith("CG"):
                    rCG = rCG + int(c[3])
                    mrCG = mrCG + int(c[2])
                    CG = CG + 1
                    mCG = mCG + int(c[4])
                elif c[1].endswith("G"):
                    rCHG = rCHG + int(c[3])
                    mrCHG = mrCHG + int(c[2])
                    CHG = CHG + 1
                    mCHG = mCHG + int(c[4])
                else:
                    rCHH = rCHH + int(c[3])
                    mrCHH = mrCHH + int(c[2])
                    CHH = CHH + 1
                    mCHH = mCHH + int(c[4])
    a = a.append({'Gene':name, 'CG_sites':CG, 'mCG_sites':mCG,
                  'CG_reads':rCG, 'mCG_reads':mrCG,
                  'CG_methylation_level':str(np.float64(mrCG)/np.float64(rCG)),
                  'CHG_sites':CHG, 'mCHG_sites':mCHG, 'CHG_reads':rCHG,
                  'mCHG_reads':mrCHG, 'CHG_methylation_level':str(np.float64(mrCHG)/np.float64(rCHG)),
                  'CHH_sites':CHH, 'mCHH_sites':mCHH, 'CHH_reads':rCHH, 'mCHH_reads':mrCHH,
                  'CHH_methylation_level':str(np.float64(mrCHH)/np.float64(rCHH)),
                  'C_sites':(CG+CHG+CHH), 'mC_sites':(mCG+mCHG+mCHH),
                  'C_reads':(rCG+rCHG+rCHH), 'mC_reads':(mrCG+mrCHG+mrCHH),
                  'C_methylation_level':str(np.float64(mrCG+mrCHG+mrCHH)/np.float64(rCG+rCHG+rCHH))
                  },ignore_index=True)
    if output:
        a.to_csv(output, sep='\t', index=False)
    else:
        return a

#
def FDR(a,column,new_col):
    b=a[a[column] != "NaN"][column].sort_values(ascending=False)
    lp=len(b)
    x=list(reversed(range(1,lp+1)))
    y=np.minimum.accumulate([lp/i*j for i,j in zip(x,b)])
    z=[i if i < 1.0 else 1.0 for i in y]
    c=pd.DataFrame(z,columns=[new_col])
    c.index=b.index
    a=pd.concat([a,c],axis=1)
    return a

#
def gene_binom_test(df,cutoff=10,calc_baseline=True,mCG=(),mCHG=(),mCHH=(),output=()):
    a=pd.read_table(df,sep="\t")
    if calc_baseline:
        mCG=pd.DataFrame.sum(a['mCG_sites'])/pd.DataFrame.sum(a['CG_sites'])
        mCHG=pd.DataFrame.sum(a['mCHG_sites'])/pd.DataFrame.sum(a['CHG_sites'])
        mCHH=pd.DataFrame.sum(a['mCHH_sites'])/pd.DataFrame.sum(a['CHH_sites'])
    elif not mCG or not mCHG or not mCHH:
        print('Use must specify baseline mCG, mCHG, and mCHH levels')
        return
    a['CG_pvalue']=stats.binom.sf(a.mCG_sites-1,a.CG_sites,mCG)
    a['CG_pvalue']=[i if s >= cutoff else "NaN" for i,s in zip(a['CG_pvalue'],a['CG_sites'])]
    a=FDR(a,column="CG_pvalue",new_col="CG_qvalue")
    a['CHG_pvalue']=stats.binom.sf(a.mCHG_sites-1,a.CHG_sites,mCHG)
    a['CHG_pvalue']=[i if s >= cutoff else "NaN" for i,s in zip(a['CHG_pvalue'],a['CHG_sites'])]
    a=FDR(a,column="CHG_pvalue",new_col="CHG_qvalue")
    a['CHH_pvalue']=stats.binom.sf(a.mCHH_sites-1,a.CHH_sites,mCHH)
    a['CHH_pvalue']=[i if s >= cutoff else "NaN" for i,s in zip(a['CHH_pvalue'],a['CHH_sites'])]
    a=FDR(a,column="CHH_pvalue",new_col="CHH_qvalue")
    if output:
        a.to_csv(output, sep='\t', index=False)
    else:
        return a

#
def classify(row,min_sites=0,qvalue=0.05):
    if row['CHH_qvalue'] <= qvalue and row['CHH_sites'] >= min_sites:
        return 'CHH'
    elif row['CHG_qvalue'] <= qvalue and row['CHG_sites'] >= min_sites:
        return 'CHG'
    elif row['CG_qvalue'] <= qvalue and row['CG_sites'] >= min_sites:
        return 'CG'
    elif (row['mCG_sites']+row['mCHG_sites']+row['mCHH_sites']) == 0 and \
         (row['CG_sites']+row['CHG_sites']+row['CHH_sites']) >= min_sites:
        return 'UM'
    else:
        return 'NA'

#
def classify_genes(df,output=(),min_sites=0,qvalue=0.05):
    a = pd.read_table(df,sep="\t")
    a['Classification'] = a.apply(classify,axis=1,min_sites=min_sites,qvalue=qvalue)
    if output:
        a.to_csv(output, sep='\t', index=False)
    else:
        return a

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

#
def genome_window_methylation(allc,genome_file,output=(),window_size=100,cutoff=0,filter_chr=[],output_mC_counts=True):
	w_bed = pbt.bedtool.BedTool.window_maker(pbt.BedTool(genome_file),g=genome_file,w=window_size,i='srcwinnum')
	mC_bed = allc2bed(allc)
	allc_mapping = pbt.bedtool.BedTool.intersect(mC_bed,w_bed,wa=True,wb=True)
	m = pd.read_table(allc_mapping.fn,header=None,usecols=[10,13,6,7,8])
	#m = m.sort_values(by = 13,ascending=True)
	b = window_methylation_levels(m,cutoff=cutoff,nuc_bed=(),output_mC_counts=True)
	if output:
		b.to_csv(output, sep='\t', index=False)
	else:
		return b
	for i in ['w_bed','mC_bed','allc_mapping','m','b']:
		del(i)
