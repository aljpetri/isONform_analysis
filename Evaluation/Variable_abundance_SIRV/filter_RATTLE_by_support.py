import argparse
import os,glob
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
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


def filter_and_write(full_outpath,ab,transcriptome):
    supp_dict={}
    consensus_file = open(full_outpath, "w")
    for acc,seq in transcriptome.items():
        #print(acc)
        reads_field=acc.split("_")[2]
        nr_reads=int(reads_field.split("=")[1])
        supp_dict[acc]=nr_reads

    for supp_id,supp in supp_dict.items():
        if supp>ab and supp_id in transcriptome:
            seq_to_write=transcriptome[supp_id]
            consensus_file.write("@{0}\n{1}\n+\n{2}\n".format(supp_id, seq_to_write,"+" * len(seq_to_write)))
def main(args):
    print("Hello World")
    #yourpath = args.indir
    clust_count = 0
    other_count = 0
    path=args.indir
    outfolder=args.outfolder
    #path = '/home/alexanderpetri/isONform_analysis/Analysis/SIRV/VarAbundance/variable_abundance/isONform/'
    #outfolder='/home/alexanderpetri/isONform_analysis/Analysis/SIRV/VarAbundance/variable_abundance/filtered'
    # Check whether the specified path exists or not
    isExist = os.path.exists(outfolder)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(outfolder)
    print(path)
    abundance_list=[2,5,10,15,20,25,50]
    unique_cl_ids = set()
    for entry in os.listdir(path):
        if entry.endswith('.fastq') and entry.startswith('rattle_res_50_'):
            name = entry.split(".")[0]
            id=name.split("_")[-1]
            unique_cl_ids.add(id)
    for u_id in unique_cl_ids:
        transcriptome_name='rattle_res_50_'+u_id+'.fastq'
        full_path_transcriptome=os.path.join(path,transcriptome_name)
        transcriptome = {acc: (seq) for acc, (seq, qual) in readfq(open(full_path_transcriptome, 'r'))}
        for ab in abundance_list:
            outfile_name="transcriptome_50_"+u_id+"_"+str(ab)+".fastq"
            full_outpath=os.path.join(outfolder,outfile_name)
            filter_and_write(full_outpath,ab,transcriptome)

        #if os.path.isfile(os.path.join(path, entry)):
        #    print("Hallo ")
    #for filename in glob.glob(os.path.join(path, '*.fastq')):
    #    print(filename)
    #    name=filename.split(".")[0]
    #    id=name.split("_")[-1]
    #    print(id)





if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        "Parses a fasta file and output all sequences longer than min_size to fasta format to /path/to/fasta_file_filtered.fa.")
    #parser.add_argument('original', type=str, help='CSV file. ')
    #parser.add_argument('rattle_name', type=str, help='csv file. ')
    #parser.add_argument('rattle_output', type=str, help='.fastq file. ')
    parser.add_argument('--indir',type=str, help='input folder ')
    parser.add_argument('--outfolder', type=str, help='output folder ')
    #parser.add_argument('--isoabundances',type=int, help="The minimum numbers of ")
    #parser.add_argument('dummy_outfile', type=str, help='csv file. ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')

    args = parser.parse_args()

    main(args)