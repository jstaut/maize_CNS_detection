import numpy as np
import pandas as pd
import pybedtools
from pybedtools import BedTool
import gzip
import subprocess

def get_gff(promoter_gff, output_name):
    
#     fasta = BedTool(read_zipped_file(fasta_file))
    promoter_gff.saveas().moveto(output_name + '.gff3')
#     promoter_gff.sequence(fi = fasta, **{"name+": True}, s = True).save_seqs(output_name + '.fasta')

def get_output_name(out_name, up, down):
    
    out = out_name + "_" + str(up) + "up" + "_" + str(down) + "down"
    
    return out

def subtract_genebody(gene_coords_df, proms_gff):
    
    #if rel_chrom != None:
        #gene_coords_df = gene_coords_df.replace({0: rel_chrom})
    #else: pass
    
    gff_cds = get_gff_object(gene_coords_df)
    final_prom = proms_gff.subtract(gff_cds).saveas()    
    
    return final_prom

def subtract_exons(gff_path, proms_gff):
    
    gff_df = pd.read_csv(gff_path, sep = '\t', header = None, comment = '#')
    
    #if rel_chrom != None:
        #gff_df = gff_df.replace({0: rel_chrom})
    #else: pass
    
    exon_df = gff_df[gff_df.loc[:, 2] == 'exon']
    UTRs_df = gff_df[(gff_df.loc[:, 2] == 'five_prime_UTR') | (gff_df.loc[:, 2] == 'three_prime_UTR')]
        
    gff_exons = BedTool().from_dataframe(exon_df)
    gff_utrs = BedTool().from_dataframe(UTRs_df)
    exons_noutr = gff_exons.subtract(gff_utrs).saveas()
    final_prom = proms_gff.subtract(exons_noutr).saveas()
        
    return final_prom

def get_gff_object(promoter_df):
    
    promoter_df.columns = ['chrom', 'start', 'end', 'strand']
    promoter_df.insert(loc = 1, column = 'source', value = 'niman')
    promoter_df.insert(loc = 2, column = 'feature', value = promoter_df.index)
    promoter_df.insert(loc = 5, column = 'score', value = '.')
    promoter_df.insert(loc = 7, column = 'frame', value = '.')
    promoter_df.insert(loc = 8, column = 'attribute', value = promoter_df.index)
    promoter_df.reset_index(inplace = True, drop = True)
    promoter_df.chrom.astype(str)
    promoter_df.start.astype(int)
    promoter_df.end.astype(int)
    prom_coord_gff = BedTool().from_dataframe(promoter_df)
    
    return prom_coord_gff

def concat_others(others_list, gff_path):
    
    other_d = {}
    
    for i in others_list:
        oth = get_other(gff_path, i)[0]
        other_d.update(oth)
    
    return other_d

def adjust_start_end(prom_coords_df, chroms_len_d):
    
    for i in np.array(prom_coords_df[0].unique()):
        index = prom_coords_df[(prom_coords_df.iloc[:, 0] == i) & (prom_coords_df.iloc[:, 2] > chroms_len_d[i])].index
        prom_coords_df.loc[index, 2] = chroms_len_d[i]
    
    prom_coords_df.loc[prom_coords_df[1] <= 0, 1] = 1
    
    return prom_coords_df

def get_sorted_df(prom_coord_dict, chrom_len_dict):
    
    prom_coord_df = pd.DataFrame().from_dict(prom_coord_dict, orient = 'index')
    
    prom_coord_df = prom_coord_df.sort_values(by = [0, 1, 2], axis = 0, ascending = [True, True, True])
    
    return prom_coord_df 


def get_sorted_gene_array(genes_dic):
    
    df = pd.DataFrame().from_dict(genes_dic, orient='index')
    df_sort = df.sort_values(by = [0, 1, 2], axis = 0, ascending = [True, True, False])
    sort_gene_names = np.array(df_sort.index)
    
    return sort_gene_names, df_sort

def read_zipped_file(file_path):
    
    gzf = gzip.open(file_path, 'rt')
    
    return gzf

def get_chrom_lengths(fasta_path):
    
    chrom_lengths_d = {}
    out = str()
    
    cmd =  "samtools faidx " + fasta_path

    p1 = subprocess.run(cmd, shell=True)
    
    with open(fasta_path + '.fai', 'r') as f:
        for line in f:
            record = line.strip().split('\t')
            chrom_lengths_d[record[0]] = int(record[1])
    
    return chrom_lengths_d

def merge_gene_dicts(L):
    genes_dict = {}
    for i in L:
        genes_dict.update(i)
    
    return genes_dict

def filter_isoforms(pcod_dic):
    '''
    Version 2
    '''
    
    coord_iso = pd.DataFrame().from_dict(pcod_dic, orient = 'index')
    coord_iso.reset_index(inplace = True)
    coord_iso = coord_iso.rename(columns = {'index': 'isoform', 0: 'chrom', 1: 'start', 2: 'end', 3: 'strand', 4: 'gene'})
    flg = all(coord_iso.isoform.str.split('.', expand = True)[0] == coord_iso.gene.str.split('.', expand = True)[0])
    
    if flg == True:

        coord_iso['dist'] = coord_iso.end - coord_iso.start + 1
        coord_iso = coord_iso.drop_duplicates(subset = ['gene', 'dist'], keep = 'first')
        coord_iso.reset_index(inplace = True, drop = True)
        coord_iso = coord_iso[coord_iso.groupby(['gene'])['dist'].transform(max) == coord_iso['dist']]
    elif flg == False:
        coord_iso.isoform = coord_iso.gene + ':' + coord_iso.isoform
        coord_iso['dist'] = coord_iso.end - coord_iso.start + 1
        coord_iso = coord_iso.drop_duplicates(subset = ['gene', 'dist'], keep = 'first')
        coord_iso.reset_index(inplace = True, drop = True)
        coord_iso = coord_iso[coord_iso.groupby(['gene'])['dist'].transform(max) == coord_iso['dist']]

    coord_iso = coord_iso[['isoform', 'chrom', 'start', 'end', 'strand']]
    coord_iso = coord_iso.set_index('isoform')
    final_isoforms = np.array(coord_iso.index)
    pcod_dic_filtered = coord_iso.T.to_dict(orient='list')

    return pcod_dic_filtered, final_isoforms

def get_pcg(gff_path):
    '''
    Version 2 - stores also gene name
    '''
    p_cod_dic = {}
    pc_genes = []
    with open(gff_path) as f:
        for line in f:
            if not line.startswith('#'):
                rec = line.strip().split('\t')
                if rec[2] == 'gene':
                    gene_name = rec[8].split(';')[0].split('=')[1].split(':')[-1]
                if rec[2] == 'mRNA':
                    name = rec[8].split(';')[0].split('=')[1].split(':')[-1]
                    flag = 1
                    flag2 = 0
                if rec[2] == 'CDS':
                    chromo = rec[0]
                    start = int(rec[3])
                    end = int(rec[4])
                    strand = rec[6]
                    if flag == 1:
                        p_cod_dic[gene_name] = [chromo, start, end, strand, gene_name]
                        pc_genes.append(gene_name)
                    else:
                        p_cod_dic[gene_name][2] = end
                    flag = 0

    return p_cod_dic, np.array(pc_genes)

def get_npcg(gff_path):
    n_pcod_dic = {}
    nc_genes = []
    flag2 = 0
    with open(gff_path) as f:
        for line in f:
            if not line.startswith('#'):
                rec = line.strip().split('\t')
                if flag2 == 1:
                    if rec[2] != 'mRNA':
                        n_pcod_dic[gene_name] = [chromo, start, end, strand]
                        nc_genes.append(gene_name)
                        flag2 = 0
                if rec[2] == 'gene':
                    flag2 = 1
                    gene_name = rec[8].split(';')[0].split('=')[1].split(':')[-1]
                    chromo = rec[0]
                    start = int(rec[3])
                    end = int(rec[4])
                    strand = rec[6]
                if rec[2] == 'mRNA':
                    flag2 = 0
                    
                    
    nc_genes = np.array(nc_genes)
    
    return n_pcod_dic, nc_genes

def get_trans(gff_path):

    trans_dic = {}
    trans_names = []
    with open(gff_path) as f:
        for line in f:
            if not line.startswith('#'):
                rec = line.strip().split('\t')
                if rec[2] == 'transposable_element':
                    t_name = rec[8].split(';')[0].split('=')[1].split(':')[-1]
                    chromo = rec[0]
                    start = int(rec[3])
                    end = int(rec[4])
                    strand = rec[6]
                    trans_dic[t_name] = [chromo, start, end, strand]
                    trans_names.append(t_name)


    trans_names = np.array(trans_names)
    
    return trans_dic, trans_names

def get_other(gff_path, other):
    other_dic = {}
    other_names = []
    with open(gff_path) as f:
        for line in f:
            if not line.startswith('#'):
                rec = line.strip().split('\t')
                if rec[2] == other:
                    t_name = rec[8].split(';')[0].split('=')[1].split(':')[-1]
                    chromo = rec[0]
                    start = int(rec[3])
                    end = int(rec[4])
                    strand = rec[6]
                    other_dic[t_name] = [chromo, start, end, strand]
                    other_names.append(t_name)


    other_names = np.array(other_names)
    
    return other_dic, other_names

def check_singleGene_chrom(gene_list, coord_dict, j):
    gene = gene_list[j]
    chrom = coord_dict[gene][0]
    try:
        n_gene = gene_list[j+1]
        chrom_next = coord_dict[n_gene][0]
    except:
        chrom_next = None
    try:
        p_gene = gene_list[j-1]
        chrom_prev = coord_dict[p_gene][0]
    except:
        chrom_prev = None
    if chrom_prev == None and chrom != chrom_next:
        return gene
    elif chrom_prev != chrom != chrom_next:
        return gene
    elif chrom != chrom_prev and chrom_next == None:
        return gene
    else:
        return None
    
def up_overlap(n, L, D, leap):
    gene = L[n]
    up_gene = L[n + leap]
    gene_start = D[gene][1]
    up_end = D[up_gene][2]
    while up_end >= gene_start:
        if D[gene][0] != D[up_gene][0]:
            up_gene = None
            break
        leap = leap - 1
        up_gene = L[n + leap]
        up_end = D[up_gene][2]
    return up_gene

def down_overlap(n, L, D, leap):
    gene = L[n]
    down_gene = L[n + leap]
    gene_end = D[gene][2]
    down_start = D[down_gene][1]
    while down_start <= gene_end:
        if D[gene][0] != D[down_gene][0]:
            down_gene = None
            break
        leap = leap + 1
        down_gene = L[n + leap]
        down_start = D[down_gene][1]
    return down_gene

def get_wind_up_down(coord_dict, gene_list, win_up, win_down):
    prom_coord = {}
    counter = 0
    counter2 = 0
    for i in range(len(gene_list)):
        
        is_single = check_singleGene_chrom(gene_list, coord_dict, i)
        
        if is_single != None:
            
            gene = is_single
            chrom = coord_dict[gene][0]
            start = coord_dict[gene][1]
            end = coord_dict[gene][2]
            strand = coord_dict[gene][3]

            if strand == '+':
                start = start - win_up
                end = end + win_down
            elif strand == '-':
                start = start - win_down
                end = end + win_up
            
            prom_coord[gene] = [chrom, start, end, strand]
            
        elif is_single == None:
            
            gene = gene_list[i]
            chrom = coord_dict[gene][0]
            start = coord_dict[gene][1]
            end = coord_dict[gene][2]
            strand = coord_dict[gene][3]
            if counter == 0:
                n_gene = down_overlap(i, gene_list, coord_dict, 1)
                if n_gene == None:
                    n_gene_s = float('inf')
                else:
                    n_gene_s = coord_dict[n_gene][1]
                if strand == '+':
                    start = start - win_up
                    end = end + win_down
                    if end >= n_gene_s:
                        end = n_gene_s - 1
                    else: pass
                elif strand == '-':
                    start = start - win_down
                    end = end + win_up
                    if end >= n_gene_s:
                        end = n_gene_s - 1
                    else: pass
                counter = 1
                counter2 = 0
            elif counter2 == 1:
                try:
                    n_gene = gene_list[i+1]
                    chrom_next = coord_dict[n_gene][0]
                    if chrom != chrom_next:
                        counter = 0
                        p_gene = up_overlap(i, gene_list, coord_dict, -1)
                        if p_gene == None:
                            p_gene_e = float('-inf')
                        else:
                            p_gene_e = coord_dict[p_gene][2]
                        if strand == '+':
                            start = start - win_up
                            end = end + win_down
                            if start <= p_gene_e:
                                start = p_gene_e + 1
                            else: pass
                        elif strand == '-':
                            start = start - win_down
                            end = end + win_up
                            if start <= p_gene_e:
                                start = p_gene_e + 1
                            else: pass
                    elif chrom == chrom_next:
                        p_gene = up_overlap(i, gene_list, coord_dict, -1)
                        if p_gene == None:
                            p_gene_e = float('-inf')
                        else:
                            p_gene_e = coord_dict[p_gene][2]
                        n_gene = down_overlap(i, gene_list, coord_dict, 1)
                        if n_gene == None:
                            n_gene_s = float('inf')
                        else:
                            n_gene_s = coord_dict[n_gene][1]
                        if strand == '+':
                            start = start - win_up
                            end = end + win_down
                            if start <= p_gene_e:
                                start = p_gene_e + 1
                            else: pass
                            if end >= n_gene_s:
                                end = n_gene_s - 1
                            else: pass
                        elif strand == '-':
                            start = start - win_down
                            end = end + win_up
                            if start <= p_gene_e:
                                start = p_gene_e + 1
                            else: pass
                            if end >= n_gene_s:
                                end = n_gene_s -1
                            else: pass
                except IndexError:
                    p_gene = up_overlap(i, gene_list, coord_dict, -1)
                    if p_gene == None:
                        p_gene_e = float('-inf')
                    else:
                        p_gene_e = coord_dict[p_gene][2]
                    if strand == '+':
                        start = start - win_up
                        end = end + win_down
                        if start <= p_gene_e:
                            start = p_gene_e + 1
                        else: pass
                    elif strand == '-':
                        start = start - win_down
                        end = end + win_up
                        if start <= p_gene_e:
                            start = p_gene_e + 1
                        else: pass
            counter2 = 1
            prom_coord[gene] = [chrom, start, end, strand]
    return prom_coord
