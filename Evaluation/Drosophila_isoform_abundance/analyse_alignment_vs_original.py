import argparse
import sys, os
import random
import pysam
def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].replace(" ", "_"), [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, (''.join(seqs), None)  # yield a fasta record
            if not last: break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs));  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, (seq, None)  # yield a fasta record instead
                break
def get_abundance_aligned_reads(sam_file):

    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    references = SAM_file.references
    valid_genes = set(["SIRV1","SIRV2","SIRV3","SIRV4", "SIRV5", "SIRV6","SIRV7"])


    transcript_dict={}

    for read in SAM_file.fetch(until_eof=True):
        if (read.flag == 0 or read.flag == 16):
            transcript_id = read.reference_name
            if not transcript_id in transcript_dict:
                transcript_dict[transcript_id]=[]
            transcript_dict[transcript_id].append(read.query_name)
            print(transcript_id,":",read.query_name)
    return transcript_dict
    #for id,reads in transcript_dict.items():
        #outfile.write("{0}:\n {1}\n".format(id, reads))
def main(args):
    #print(args)
    dirname = os.path.dirname(args.outfile)
    #print(args.sirv_ref)
    #print(args.fastq)
    #mkdir_p(dirname)
    #file_exists=os.path.exists(os.path.join(dirname,"minimap.txt"))
    #if not file_exists:
    result_dict={}
    transcript_dict = get_abundance_aligned_reads(args.alignment)
    predictions = {acc: (seq) for acc, (seq, qual) in readfq(open(args.fastq, 'r'))}
    for acc, value in predictions.items():
        if acc in transcript_dict:
            result_dict[acc]=transcript_dict[acc]
        else:
            result_dict[acc]=[]
    out=open(args.outfile,"w")
    for key,value in result_dict.items():
        out.write("{}: {}\n".format(key,value))
    out.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser("Parses alignments and subsamples a fastq file based on alignments.")
    parser.add_argument('--fastq', type=str, help='fastq file. ')
    parser.add_argument('--alignment', type=str, help='fastq file. ')
    parser.add_argument('--outfile', type=str, help='Fastq file. ')
    args = parser.parse_args()
    main(args)