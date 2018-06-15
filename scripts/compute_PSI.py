import argparse
import pandas as pd
import numpy as np
import math

def compute_gene_PSI(PSI_df,exon_gene_df,junction_gene_df,annot_gene_df,read_len,norm):
    'Compute Percent Spliced In (PSI)'
    exons = list(exon_gene_df.index)
    min_reads = 10

    for exon in exons:

        exon_chrom = annot_gene_df.loc[exon]['probe.chr'].split('chr')[1]
        exon_start = annot_gene_df.loc[exon]['probe.start']
        exon_stop = annot_gene_df.loc[exon]['probe.stop']
        exon_len = math.fabs(exon_start - exon_stop) + 1

        B_intervals = []
        C_intervals = []
        for interval in junction_gene_df.index:
            chrom,start,stop = interval.split('_')
            if chrom != exon_chrom:
                print(chrom)
                print(exon_chrom)
                print('Chromosomes aren\'t matching, something is wrong!!')

            if int(start) < exon_start and int(stop) > exon_stop:#C
                C_intervals.append(interval)
            if int(start) < exon_start and (int(stop) >= exon_start and int(stop) <= exon_stop):#B1
                B_intervals.append(interval)
            if (int(start) >= exon_start and int(start) <= exon_stop) and int(stop) > exon_stop:#B2
                B_intervals.append(interval)

        A_reads = exon_gene_df.loc[exon]
        B_reads = junction_gene_df.loc[B_intervals]
        B_reads = B_reads.iloc[:,2:].sum()

        C_reads = junction_gene_df.loc[C_intervals]
        C_reads = C_reads.iloc[:,2:].sum()

        total_reads = A_reads + B_reads + C_reads
        if len(total_reads[total_reads < min_reads])/len(total_reads) > 0.2:
            continue
        else:
            A_B_norm = (A_reads + B_reads)/(read_len + exon_len-1)
            C_norm = C_reads/(read_len-1)

            mask = (A_B_norm == 0) & (C_norm == 0)
            PSI_norm = 100*((A_B_norm)/(A_B_norm + C_norm)) if (A_B_norm.sum() + C_norm.sum()) > 0.00 else 0
            PSI_norm = PSI_norm.where(~mask,0)
            if norm:
                PSI_norm = ZNorm(PSI_norm)
                if PSI_norm == None:
                    print('{} had no variance'.format(exon))
                    continue
            PSI_df[exon] = PSI_norm

    return PSI_df

def ZNorm(vals):
    m = np.mean(vals)
    sd = math.sqrt(np.var(vals))
    if sd == 0: return None
    return [(item-m)/sd for item in vals]

def max_coverage_samples(exon_df):
    'Given exon expression dataframe'
    'Take first two elements of sample name'
    'If there is duplicate sample'
    'Return list of max coverage non-dup samples'
    
    med_exon_df = pd.DataFrame(exon_df.median(axis=0))
    med_exon_df['samples'] = ["-".join(i.split('-')[:2]) for i in med_exon_df.index]
    med_exon_df['original_samples'] = med_exon_df.index
    highest_cov_samples = list(med_exon_df.groupby('samples').idxmax()[0]) 
    
    return highest_cov_samples

def filter_exons(PSI_df):
    'Given PSI dataframe'
    'Filter out exons with little variation:'
    'PSI = 100% or PSI = 0% in >90% patients'
    'Return filtered PSI'
    
    PSI_filtered_df = PSI_df.dropna(axis=1)
    exons = PSI_filtered_df.columns
    exons_to_remove = []
    num_samples = len(PSI_df.index)
    for exon in exons:
        exon_PSI = PSI_filtered_df[exon]
        num_zero_or_hunded = len(exon_PSI[exon_PSI == 0]) + len(exon_PSI[exon_PSI == 100])
        if num_zero_or_hunded/num_samples > 0.9:
            exons_to_remove.append(exon)
    PSI_filtered_df = PSI_filtered_df.drop(exons_to_remove,axis=1).dropna()
    new_index = ["-".join(i.split('-')[:2]) for i in PSI_filtered_df.index]
    PSI_filtered_df.index = new_index

    return PSI_filtered_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute Percent Spliced In (PSI) for any exon")
    parser.add_argument("--exon_expr", help="File path to exon expression dataset", type=str, required=True)
    parser.add_argument("--junction_expr", help="File path to junction expression dataset", type=str, required=True)
    parser.add_argument("--exprannot", help="Exon expression annotation file", type=str, required=True)
    parser.add_argument("--chrom", help="Chromsome for PSI calculations", type=str, required=True)
    parser.add_argument("--min_samples", help="Require data for this many samples", type=int, default=0)
    parser.add_argument("--read_len", help="Paired end length, default is GTEx read length: 75bp ", type=int, default=75)
    parser.add_argument("--norm", help="Normalize PSI values", required=False)
    parser.add_argument("--out", help="Write data files to this file", type=str, required=True)

    args = parser.parse_args()
    EXONFILE = args.exon_expr
    JUNCTIONFILE = args.junction_expr
    ANNOTFILE = args.exprannot
    CHROM = args.chrom
    MINSAMPLES = args.min_samples
    READLEN = args.read_len
    OUTFILE = args.out
    NORM = True if args.norm else False        

    exon_df = pd.read_table(EXONFILE, index_col=0)
    annot_df = pd.read_table(ANNOTFILE, index_col=0,sep=',')
    junction_df = pd.read_table(JUNCTIONFILE, index_col=0,
                                          sep='\t',
                                          skipinitialspace=True)
    if CHROM:
        if CHROM.startswith('chr'):
            CHROM = CHROM.split('chr')[1]
        annot_df = annot_df[annot_df['gene.chr'] == 'chr' + CHROM]
        junction_df =junction_df[junction_df['Gene_Symbol'].isin(annot_df['gene.id'])]
        exon_df = exon_df.loc[list(annot_df['probe.id'])]
    #keeping only samples which are in both exon and junction files:
    common_samples = list(junction_df.columns.intersection(exon_df.columns))
    exon_df = exon_df[common_samples]
    common_samples[:0] = ['Gene_Symbol','Chr']
    junction_df = junction_df[common_samples]

    #removing duplicated samples by taking highest coverage (in exon expression) sample

    highest_cov_samples = max_coverage_samples(exon_df) 
    exon_df = exon_df[highest_cov_samples]
    highest_cov_samples[:0] = ['Gene_Symbol','Chr']
    junction_df = junction_df[highest_cov_samples]
    PSI_df = pd.DataFrame(index=exon_df.columns, columns=exon_df.index)
    genes = set(annot_df['gene.id'])
    print('Chromosome: {}'.format(CHROM))
    print('Exons total: {}'.format(len(exon_df.index)))
    
    for gene in genes:
        exon_gene_df = exon_df.loc[exon_df.index.str.startswith(gene)]
        junction_gene_df = junction_df[junction_df['Gene_Symbol'] == gene]
        annot_gene_df = annot_df[annot_df['gene.id'] == gene]
        annot_gene_df.index = annot_gene_df['probe.id']
        PSI_df = compute_gene_PSI(PSI_df,
                                  exon_gene_df,
                                  junction_gene_df,
                                  annot_gene_df,
                                  READLEN,
                                  NORM)
    PSI_filtered_df = filter_exons(PSI_df.dropna(axis=1))
    print('Exons with PSI: {}\n'.format(len(PSI_filtered_df.columns)))
    PSI_filtered_df.to_csv(OUTFILE)