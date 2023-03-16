import argparse
import sys, os
import random
import pysam
import errno
from collections import defaultdict
import parasail
from collections import defaultdict
import edlib
import argparse
import errno
import re
def get_abundance_aligned_reads(sam_file,outfile):

    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    references = SAM_file.references
    valid_genes = set(["SIRV1","SIRV2","SIRV3","SIRV4", "SIRV5", "SIRV6","SIRV7"])
    #transcript_cov = {}
    transcript_cov = defaultdict(set)
    # amgiguous_primary = defaultdict(set)
    #sam_dict=SAM_file.to_dict()
    #print(sam_dict)
    outfile = open(outfile, "w")
    transcript_dict={}

    for read in SAM_file.fetch(until_eof=True):
        if (read.flag == 0 or read.flag == 16):
            transcript_id = read.reference_name
            if not transcript_id in transcript_dict:
                transcript_dict[transcript_id]=[]
            transcript_dict[transcript_id].append(read.query_name)
            print(transcript_id,":",read.query_name)
    for id,reads in transcript_dict.items():
        outfile.write("{0}:\n {1}\n".format(id, reads))
            #if transcript_id[:5] in valid_genes:
                #gene_id = transcript_id[4]
                #transcript_cov[transcript_id].add(read.query_name)
                #print(read.cigarstring)


def main(args):
    #print(args)
    #dirname = os.path.dirname(args.outfile)
    #print(args.sirv_ref)
    #print(args.fastq)
    #mkdir_p(dirname)
    #file_exists=os.path.exists(os.path.join(dirname,"minimap.txt"))
    #if not file_exists:
    transcript_cov = get_abundance_aligned_reads(args.alignments,args.outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Parses alignments and subsamples a fastq file based on alignments.")
    #parser.add_argument('--fastq', type=str, help='fastq file. ')
    #parser.add_argument('--sirv_ref', type=str, help='fastq/a file. ')
    parser.add_argument('--alignments', type=str, help='fastq file. ')
    parser.add_argument('--outfile', type=str, help='Mapping file. ')
    #parser.add_argument('--depth_dist', type=str, help="type of read distribution: either uniform (20 reads per isoform) or exp (2^(x-2)|x<=10)")
    #parser.add_argument('--nr_isoforms', type=int, help='Number of nr_isoforms.')
    args = parser.parse_args()
    main(args)