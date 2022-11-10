#!/usr/bin/env python

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
    print("REF",reference)
    depthsum=0
    multi=0
    for i in range(0,nr_isoforms):
        if depth_dist=="uniform":
            depth=random.randrange(10,50)
        elif depth_dist=="exp":
            depth=random.choice([8,16,32,64,128])
        depth_list.append(depth)
    for depth_elem in depth_list:  
    #for gene_id in transcript_cov:
        not_set = True
        while not_set:
            tr_acc = random.sample(list(transcript_cov), 1)
            this_acc=tr_acc[0]
            #print("THIS",this_acc)
            if len(transcript_cov[this_acc]) > depth_elem:
                #print(this_acc)
                #print("LEN",len(transcript_cov[this_acc]))
                if not this_acc in already_used:
                    subset = random.sample(transcript_cov[this_acc], depth_elem) # this is without replacement
                    #print(subset)
                    subsamples[this_acc] = subset
                    not_set=False
                    already_used.add(this_acc)
                    new_reference[this_acc]=reference[this_acc]
                    print("REFSEQ",reference[this_acc])
                    #print("Sampled ", depth_elem," reads for ",this_acc)
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
    #transcript_cov = {}
    transcript_cov = defaultdict(set)
    # amgiguous_primary = defaultdict(set)

    for read in SAM_file.fetch(until_eof=True):
        if (read.flag == 0 or read.flag == 16):
            transcript_id = read.reference_name
            if transcript_id[:5] in valid_genes:
                #gene_id = transcript_id[4]
                transcript_cov[transcript_id].add(read.query_name)
                #print(read.cigarstring)
                """if not transcript_id in transcript_cov:
                    e_list=[]
                    e_list.append(read.query_name)
                    transcript_cov[transcript_id]=e_list
                else:
                    lista=transcript_cov[transcript_id]
                    lista.append(read.query_name)
                    up_dict={transcript_id: lista}
                    transcript_cov.update(up_dict)"""
    # except:
    #     pass

    # print(transcript_cov)

        # elif (read.flag == 0 or read.flag == 16) and read.mapping_quality == 0:
        #     transcript_id = read.reference_name
        #     amgiguous_primary[transcript_id].add(read.query_name)
    print(len(transcript_cov))
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
    print(new_reference)
    for read_acc,seq in new_reference.items():
            seq = new_reference[read_acc]
            new_ref_file.write(">{0}\n{1}\n".format(read_acc, seq))
    new_ref_file.close()


def cigar_to_seq(cigar, query, ref):
    cigar_tuples = []
    result = re.split(r'[=DXSMI]+', cigar)
    cig_pos = 0
    for length in result[:-1]:
        cig_pos += len(length)
        type_ = cigar[cig_pos]
        cig_pos += 1
        cigar_tuples.append((int(length), type_))

    r_index = 0
    q_index = 0
    q_aln = []
    r_aln = []
    for length_, type_ in cigar_tuples:
        if type_ == "=" or type_ == "X":
            q_aln.append(query[q_index: q_index + length_])
            r_aln.append(ref[r_index: r_index + length_])

            r_index += length_
            q_index += length_

        elif type_ == "I":
            # insertion w.r.t. reference
            r_aln.append('-' * length_)
            q_aln.append(query[q_index: q_index + length_])
            #  only query index change
            q_index += length_

        elif type_ == 'D':
            # deletion w.r.t. reference
            r_aln.append(ref[r_index: r_index + length_])
            q_aln.append('-' * length_)
            #  only ref index change
            r_index += length_

        else:
            # print("error")
            # print(cigar)
            sys.exit()

    return "".join([s for s in q_aln]), "".join([s for s in r_aln]), cigar_tuples


def parse_cigar_diversity(cigar_tuples,delta_perc,delta_len):
    miss_match_length=0
    alignment_len=0
    #print("Now we are parsing....")
    #print(cigar_tuples)
    too_long_indel = False
    for i, elem in enumerate(cigar_tuples):

        cig_len = elem[0]
        cig_type = elem[1]
        alignment_len += cig_len
        if (cig_type != '=') and (cig_type != 'M') and (cig_type!='S'):
            #we want to add up all missmatches to compare to sequence length
            miss_match_length += cig_len
            #we also want to make sure that we do not have too large internal sequence differences
            if cig_len > delta_len:
                # print(cig_len, delta_len, cig_type )
                #if i>1 and i<(len(cigar_tuples)-1):
                    #print("ELE",elem)
                    #print("No pop due to delta_len")
                return False
    return True
def parasail_alignment(s1, s2,new_transcript_cov,transcript_id,read, match_score = 2, mismatch_penalty = -2, opening_penalty = 2, gap_ext = 1, ends_discrepancy_threshold=0):
    #print("s1", s1)
    #print("s2",s2)
    user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
    result = parasail.sg_trace_scan_16(s1, s2, opening_penalty, gap_ext, user_matrix)
    if result.saturated:
        #print("SATURATED!",len(s1), len(s2))
        result = parasail.sg_trace_scan_32(s1, s2, opening_penalty, gap_ext, user_matrix)
        #print("computed 32 bit instead")

    # difference in how to obtain string from parasail between python v2 and v3...
    if sys.version_info[0] < 3:
        cigar_string = str(result.cigar.decode).decode('utf-8')
    else:
        cigar_string = str(result.cigar.decode, 'utf-8')
    #print(cigar_string)
    #print(s1)
    #print(s2)
    delta=100
    delta_len=10
    s1_alignment, s2_alignment, cigar_tuples = cigar_to_seq(cigar_string, s1, s2)
    good_to_pop = parse_cigar_diversity(cigar_tuples, delta, delta_len)
    if good_to_pop:
        if not transcript_id in new_transcript_cov:
            new_transcript_cov[transcript_id]=[]
            new_transcript_cov[transcript_id].append(read)
        else:
            new_transcript_cov[transcript_id].append(read)
    # else:
    #     print(read, cigar_tuples)

"""def parasail_alignment(read, reference, x_acc="", y_acc="", match_score=2, mismatch_penalty=-2, opening_penalty=2,
                       gap_ext=1, ends_discrepancy_threshold=0):
    user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
    result = parasail.sg_trace_scan_16(read, reference, opening_penalty, gap_ext, user_matrix)
    if result.saturated:
        print("SATURATED!")
        result = parasail.sg_trace_scan_32(read, reference, opening_penalty, gap_ext, user_matrix)
    if sys.version_info[0] < 3:
        cigar_string = str(result.cigar.decode).decode('utf-8')
    else:
        cigar_string = str(result.cigar.decode, 'utf-8')
    print(cigar_string)
    read_alignment, ref_alignment = cigar_to_seq(cigar_string, read, reference)
    return read_alignment, ref_alignment"""

def filter_reference(transcript_cov,reference,fastq):
    throw_out={}
    i=0
    iter=0
    overall_len=len(transcript_cov)
    new_transcript_cov = {}
    for transcript_id, read_list in transcript_cov.items():
        iter+=1
        print("Iter", iter, "/", overall_len)
        for read in read_list:
            i+=1

            read_seq=fastq[read][0]
            ref_seq=reference[transcript_id]
            #print("READSEQ",read_seq)
            #print("REFSEQ",ref_seq)
            parasail_alignment(read_seq,ref_seq, new_transcript_cov,transcript_id,read)
            #insertions = ref_alignment.count("-")
            #deletions = read_alignment.count("-")
            #indels = insertions + deletions
            #mismatches = len([1 for n1, n2 in zip(read_alignment, ref_alignment) if n1 != n2 and n1 != "-" and n2 != "-"])
            #edit_distance = deletions + insertions + mismatches
    return new_transcript_cov

def main(args):
    #print(args)
    dirname = os.path.dirname(args.outfile)
    #print(args.sirv_ref)
    #print(args.fastq)
    mkdir_p(dirname)
    transcript_cov = get_abundance_aligned_reads(args.alignments)
    cov_file=open(os.path.join(dirname,"minimap.txt"),'w')

    sum_good = 0
    for key,value in transcript_cov.items():
        sum_good += len(value)
        print(key,len(value))
        cov_file.write("{0}:{1}\n".format(key,value))
    cov_file.close()
    print("Total transcripts minimap2 aligments only:", sum_good)
    #print("transcript_cov ",transcript_cov)
    reference = {acc: (seq) for acc, (seq, qual) in readfq(open(args.sirv_ref, 'r'))}
    #print(args.fastq)
    fastq = {acc: (seq, qual) for acc, (seq, qual) in readfq(open(args.fastq, 'r'))}
    #print(len(transcript_cov), [len(transcript_cov[g]) for g in transcript_cov])
    new_reference={}
    transcript_cov=filter_reference(transcript_cov,reference,fastq)
    
    print()
    sum_good = 0
    for key,value in transcript_cov.items():
        sum_good += len(value)
        print(key,len(value))
    print("Total transcripts after filtering:", sum_good)
    subsamples = get_subsamples(transcript_cov, args.depth_dist, args.nr_isoforms,reference,new_reference)
    print("NEWFRE",new_reference)
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
