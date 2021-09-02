# this script was designed based on Periscope and artic
# https://github.com/sheffield-bioinformatics-core/periscope


from artic.align_trim import find_primer
from artic.vcftagprimersites import read_bed_file

import argparse
import os
import glob
import pysam
import re
import pandas as pd

def get_mapped_reads(bam):
    # find out how many mapped reads there are for bam # periscope
    mapped_reads = int(pysam.idxstats(bam).split("\n")[0].split("\t")[2])
    return mapped_reads


# bam="/Volumes/TOSHIBA/artic/animal_sample/periscope//0628MN002-14BC12/0628MN002-14BC12_gRNA.bam"
# primer_bed="/Users/yaosmacbook/Downloads/periscope-master/periscope/resources/artic_primers_V3.bed"

def main(args):

    if not args.bam:
            print("no bam file")
            exit(1)

    if not args.primer_bed:
        print("no primer bed file")
        exit(1)
        
        
    file_dir = os.path.split(args.bam)[0]
    
    sample_name_full = os.path.split(args.bam)[1]
    sample_name = re.search( '([A-Za-z0-9-_]+).bam', sample_name_full).group(1)
    
    
    inbamfile = pysam.AlignmentFile( args.bam, 'rb' )
    primer_bed_object = read_bed_file(args.primer_bed)

    mapped_reads = get_mapped_reads(args.bam)
    print( "mapped reads: %s" %(mapped_reads), file=open(os.path.join( file_dir, sample_name +".txt" ), 'w') )

    pos=[]
    lth=[]
    l_primer=[]
    r_primer=[]

    for read in inbamfile:
        if read.is_unmapped != True:

            left_primer = find_primer(primer_bed_object, read.reference_start, '+')
            right_primer = find_primer(primer_bed_object, read.reference_end, '-')

            left_amplicon = int( left_primer[2]['Primer_ID'].split("_")[1])
            right_amplicon = int(right_primer[2]['Primer_ID'].split("_")[1])

            if read.qlen > args.len_cutoff and right_amplicon >= left_amplicon:
                pos.append( read.pos )
                lth.append( read.qlen )
                l_primer.append( left_amplicon )
                r_primer.append( right_amplicon )
                
    dict_read = {'sample': sample_name_full, 'pos':pos, 'lth': lth, 'l_primer':l_primer, 'r_primer': r_primer}
    dict_read_df = pd.DataFrame( dict_read )
    
    dict_read_df.to_csv( os.path.join( file_dir, sample_name +".csv" ) )


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='sum_bam: to log read info')
    parser.add_argument('-b', '--bam-file', dest='bam', help="The bam file of artic reads")
    parser.add_argument('-p', '--primer-bed', dest='primer_bed', help="The bed file with artic primer positions")
    parser.add_argument('-l', '--length-cutoff',dest='len_cutoff', help='Cut-off for read length', default=500)
    parser.add_argument('-s', '--sample',dest='sample', help='Given sample name', default=None)

    args = parser.parse_args()

    out = main(args)

    print("done")