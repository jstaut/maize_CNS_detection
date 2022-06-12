import argparse

def parse_command_line():

    '''
    Function that parses the arguments of the command line using the
    python module argparse
    '''

    parser = argparse.ArgumentParser(prog = 'Promoter extraction',
                        description = 'Script for extraction of promoter sequences '
                        'from a GFF and a genome fasta files as inputs.',
                        conflict_handler='resolve')

    parser.add_argument('gff', nargs = 1, type = str,
                        help = '',
                        metavar = 'GFF formatted file')

    parser.add_argument('genome_fasta', nargs = 1, type = str,
                        help = '',
                        metavar = 'FASTA formatted genome')

    parser.add_argument('-up', '--up_window', nargs = '?', type = int,
                        default = 5000,
                        metavar = 'Upstream window to extract promoter sequence')
    
    parser.add_argument('-down', '--down_window', nargs = '?', type = int,
                        default = 1000,
                        metavar = 'Downstream window to extract promoter sequence')
    
    parser.add_argument('-o', '--outputfile', nargs = '?', type = str,
                        default = 'promoter',
                        metavar = 'name of the output file without extension')
    
    group = parser.add_mutually_exclusive_group()
    
    group.add_argument('-int', '--introns', action = 'store_true',
                        help = '')
    
    group.add_argument('-adj', '--only_adjacent', action = 'store_true',
                        help = '')
    
    parser.add_argument('-npc', '--non_protein_coding', action = 'store_true',
                        help = '')
 
    parser.add_argument('-iso', '--isoforms', action = 'store_true',
                        help = '')

    parser.add_argument('-trans', '--transposable_elements', action = 'store_true',
                        help = '')
    
    parser.add_argument('-other', '--other_features', action = 'append',
                        nargs = "?", type = str,
                        help = '')
    
    parser.add_argument('-othergene', '--other_gene_like_features', action = 'append',
                        nargs = "?", type = str,
                        help = '')

    parser.add_argument('-stats', '--compute_statistics', action = 'store_true',
                        help = '')
    
    args = parser.parse_args()

    return args