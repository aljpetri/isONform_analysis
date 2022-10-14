#!/usr/bin/env python

import argparse
import sys, os
import random
import pysam
import errno
from collections import defaultdict

import argparse
import errno
"""
sirv_subsample_isoform: A script to subsample reads from a SIRV dataset according to a minimap alignment
USAGE:  python sirv_subsample_isoform.py --fastq {inputfile} --alignments {alignmentfile} --outfile {outfile}

"""
'''
    Below code taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''
#TODO: generate fasta file holding the reference. I.e. each SIRV that we choose should be added to the reference file
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].split()[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break


def get_subsamples(transcript_cov, depth_dist, nr_isoforms,reference,new_reference):
    #print("TC",transcript_cov)
    subsamples = {}
    already_used=set()
    depth_list=[]
    #print(reference)

    for i in range(0,nr_isoforms):
        if depth_dist=="uniform":
            depth=random.randrange(10,50)
        elif depth_dist=="exp":
            depth=random.choice([1,2,4,8,16,32,64,128])
        depth_list.append(depth)
    for depth_elem in depth_list:
    #for gene_id in transcript_cov:
        not_set = True
        while not_set:
            tr_acc = random.sample(list(transcript_cov), 1)
            this_acc=tr_acc[0]
            #print("THIS",this_acc)
            if len(transcript_cov[this_acc]) > depth_elem:
                #print("LEN",len(transcript_cov[this_acc]))
                if not this_acc in already_used:
                    subset = random.sample(transcript_cov[this_acc], depth_elem) # this is without replacement
                    subsamples[this_acc] = subset
                    not_set=False
                    already_used.add(this_acc)
                    new_reference[this_acc]=reference[this_acc]
                    print("Sampled ", depth_elem," reads for ",this_acc)
    return subsamples

    # subsamples = {}
    # already_ampled = set()
    # for tr_acc in transcript_cov:
    #     print(tr_acc, len(transcript_cov[tr_acc]))
    #     subset = random.sample(transcript_cov[tr_acc], depth) # this is without replacement
    #     subsamples[tr_acc] = subset
    # return subsamples



def get_abundance_aligned_reads(sam_file):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    references = SAM_file.references
    valid_genes = set(["SIRV1","SIRV2","SIRV3","SIRV4", "SIRV5", "SIRV6","SIRV7"])
    transcript_cov = {}
    # amgiguous_primary = defaultdict(set)

    for read in SAM_file.fetch(until_eof=True):
        if (read.flag == 0 or read.flag == 16):
            transcript_id = read.reference_name
            if transcript_id[:5] in valid_genes:
                #gene_id = transcript_id[4]
                #transcript_cov[gene_id][transcript_id].add(read.query_name)
                if not transcript_id in transcript_cov:
                    e_list=[]
                    e_list.append(read.query_name)
                    transcript_cov[transcript_id]=e_list
                else:
                    lista=transcript_cov[transcript_id]
                    lista.append(read.query_name)
                    up_dict={transcript_id: lista}
                    transcript_cov.update(up_dict)
                #transcript_cov[transcript_id].add(read.query_name)
    # except:
    #     pass

    # print(transcript_cov)

        # elif (read.flag == 0 or read.flag == 16) and read.mapping_quality == 0:
        #     transcript_id = read.reference_name
        #     amgiguous_primary[transcript_id].add(read.query_name)

    return transcript_cov

def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def write_new_reference(dirname,new_reference):
    ref_file_name=os.path.join(dirname,"reference.fasta")
    new_ref_file = open(ref_file_name, "w")
    for read_acc,seq in new_reference.items():
            seq = new_reference[read_acc]
            new_ref_file.write(">{0}\n{1}\n".format(read_acc, seq))
    new_ref_file.close()
def main(args):
    #print(args)
    dirname = os.path.dirname(args.outfile)
    #print(args.sirv_ref)
    #print(args.fastq)
    mkdir_p(dirname)
    transcript_cov = get_abundance_aligned_reads(args.alignments)
    reference = {acc: (seq) for acc, (seq, qual) in readfq(open(args.sirv_ref, 'r'))}
    #print(len(transcript_cov), [len(transcript_cov[g]) for g in transcript_cov])
    new_reference={}
    subsamples = get_subsamples(transcript_cov, args.depth_dist, args.nr_isoforms,reference,new_reference)
    fastq = { acc : (seq,qual) for acc, (seq,qual) in readfq(open(args.fastq, 'r'))}
    write_new_reference(dirname,new_reference)
    #print(reference)
    outfile = open(args.outfile, "w")
    for tr_id, set_of_reads in subsamples.items():
        for read_acc in set_of_reads:
            seq, qual  = fastq[read_acc]
            outfile.write("@{0}\n{1}\n+\n{2}\n".format(read_acc, seq, qual))

    outfile.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Parses alignments and subsamples a fastq file based on alignments.")
    parser.add_argument('--fastq', type=str, help='fastq file. ')
    parser.add_argument('--sirv_ref', type=str, help='fastq/a file. ')
    parser.add_argument('--alignments', type=str, help='fastq file. ')
    parser.add_argument('--outfile', type=str, help='Fastq file. ')
    parser.add_argument('--depth_dist', type=str, help="type of read distribution: either uniform (20 reads per isoform) or exp (2^(x-2)|x<=10)")
    parser.add_argument('--nr_isoforms', type=int, help='Number of nr_isoforms.')
    args = parser.parse_args()
    main(args)
