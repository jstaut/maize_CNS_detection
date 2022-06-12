from promoter_extraction import *
from arguments import parse_command_line

args = parse_command_line()

genes_coords = []
pcg_dic = get_pcg(args.gff[0])[0]

if args.isoforms == False:
    
    pcg_dic = filter_isoforms(pcg_dic)[0]
    
elif args.isoforms == True:
    
    pass

genes_coords.append(pcg_dic)

if args.non_protein_coding == True:
    
    npcg_dic = get_npcg(args.gff[0])[0]
    
    genes_coords.append(npcg_dic)

elif args.non_protein_coding == False:
    
    pass


if args.other_gene_like_features != None:
    
    others_gene_dic = concat_others(args.other_gene_like_features, args.gff[0])
    
    genes_coords.append(others_gene_dic)
    
elif args.other_gene_like_features == None:
    
    pass

genes_coords_D = merge_gene_dicts(genes_coords)

gene_names = get_sorted_gene_array(genes_coords_D)

gene_names_L = gene_names[0]

gene_names_df = gene_names[1]

prom_coord_D = get_wind_up_down(genes_coords_D, gene_names_L, args.up_window, args.down_window)


if args.transposable_elements == True:

    trans_dic = get_trans(args.gff[0])[0]
    
    prom_coord_D.update(trans_dic)
    
elif args.transposable_elements == False:
    
    pass


if args.other_features != None:
    
    others_dic = concat_others(args.other_features, args.gff[0])
    
    prom_coord_D.update(others_dic)
    
elif args.other_features == None:
    
    pass

chrom_len_D = get_chrom_lengths(args.genome_fasta[0])

prom_coord_df  = get_sorted_df(prom_coord_D, chrom_len_D)

prom_coord_df = adjust_start_end(prom_coord_df, chrom_len_D)

promoters_gff = get_gff_object(prom_coord_df)

if args.introns == True and args.only_adjacent == False:
    
    promoters_gff = subtract_exons(args.gff[0], promoters_gff)

elif args.introns == False and args.only_adjacent == True:
    
    promoters_gff = subtract_genebody(gene_names_df, promoters_gff)
    
elif args.introns == False and args.only_adjacent == False:
    pass

output_file_name = get_output_name(args.outputfile, args.up_window, args.down_window)

get_gff(promoters_gff, output_file_name)
