import pandas as pd
from pybedtools import BedTool
import sys

def subtract_exons(gff_path, output_name):

    gff_df = pd.read_csv(gff_path, sep = '\t', header = None, comment = '#')

    exon_df = gff_df[gff_df.loc[:, 2] == 'exon']
    UTRs_df = gff_df[(gff_df.loc[:, 2] == 'five_prime_UTR') | (gff_df.loc[:, 2] == 'three_prime_UTR')]

    if UTRs_df.empty:
        gff_exons = BedTool().from_dataframe(exon_df)
        gff_utrs = BedTool().from_dataframe(UTRs_df)
        exons = gff_exons.sort().saveas().moveto(output_name)
    else:
        try:
            gff_exons = BedTool().from_dataframe(exon_df)
            gff_utrs = BedTool().from_dataframe(UTRs_df)
            exons_noutr = gff_exons.sort().subtract(gff_utrs.sort()).saveas().moveto(output_name)
        except:
            gff_exons = BedTool().from_dataframe(exon_df)
            gff_utrs = BedTool().from_dataframe(UTRs_df)
            exons_noutr = gff_exons.sort().subtract(gff_utrs.sort(), sorted = True).saveas().moveto(output_name)


def main():

    gff = sys.argv[1]

    out = sys.argv[2]

    subtract_exons(gff, out)


if __name__ == "__main__":
    main()
