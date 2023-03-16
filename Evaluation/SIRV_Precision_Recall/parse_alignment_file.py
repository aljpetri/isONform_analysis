import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd



def main(args):
    print("Hello World")
    yourpath = args.indir
    clust_count=0
    other_count=0
    for root, dirs, files in os.walk(yourpath, topdown=False):
        for name in files:
            file1 = open(os.path.join(root, name), 'r')
            print("NAME",name)
            Lines = file1.readlines()
            id_dict={}
            for line in Lines[1:]:
                line_entries=line.split(",")
                isON_id=line_entries[0]
                SIRV_id=line_entries[1]
                if not SIRV_id in id_dict:
                    id_dict[SIRV_id]=[]
                    id_dict[SIRV_id].append(isON_id)
                else:
                    id_dict[SIRV_id].append(isON_id)
            for key,value in id_dict.items():
                if len(value)>1:
                    print("key:",key)
                    print(*value, sep = ", ")
                    cluster_set=set()
                    for val in value:
                        cl_id=val.split("_")[0]
                        print(cl_id)
                        cluster_set.add(cl_id)
                    if len(cluster_set)>1:
                        print("isONclust")
                        this_clusters=len(cluster_set)-1
                        clust_count+=this_clusters
                        print("clust_count",this_clusters)
                        other_count+=len(value)-this_clusters
                    else:
                        print("other")
                        other_count+=len(value)
    print("isONclust:",clust_count)
    print("OTHER_count",other_count)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        "Parses a fasta file and output all sequences longer than min_size to fasta format to /path/to/fasta_file_filtered.fa.")
    #parser.add_argument('original', type=str, help='CSV file. ')
    #parser.add_argument('rattle_name', type=str, help='csv file. ')
    #parser.add_argument('rattle_output', type=str, help='.fastq file. ')
    parser.add_argument('--indir',type=str, help='csv file. ')
    #parser.add_argument('--outfolder', type=str, help='csv file. ')
    #parser.add_argument('--max_isoforms',type=int, help="The maximum number of isoforms we generate in our files")
    #parser.add_argument('dummy_outfile', type=str, help='csv file. ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')

    args = parser.parse_args()

    main(args)


