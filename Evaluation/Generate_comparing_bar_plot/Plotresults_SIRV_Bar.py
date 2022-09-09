import numpy as np
import matplotlib.pyplot as plt
from csv import reader
import argparse
import sys
import os
import csv

def read_isonform_file(isonform_file):
        # open file in read mode
        isonform_dict = {}
        # open file in read mode
        with open(isonform_file) as f:
            lines = f.readlines()
        for line in lines:
            if not line.startswith("read"):
                row = line.split(",")
                print(row)
                value = int(row[7].replace("\n", ""))
                print(value)
                if not row[1] in isonform_dict:
                    isonform_dict[row[1]] = value

                else:
                    old_count=isonform_dict[row[1]]
                    print("OLD", old_count)
                    new_count=old_count+value
                    isonform_dict.update({row[1]:new_count})
                    print("NEW",isonform_dict[row[1]])
        return isonform_dict
def read_original_analysis(filename):
    # open file in read mode
    original_dict={}
    # open file in read mode
    with open(filename) as f:
        lines = f.readlines()
    for line in lines:
        row=line.split(",")
        value=int(row[1].replace("\n",""))
        original_dict[row[0]]=value
    return original_dict
def read_rattle_analysis(rattle_name,rattle_output):
    support_dict={}
    print("RATTLE_OUTPUT",rattle_output)
    with open(rattle_output) as r_out:
        lines = r_out.readlines()
        for line in lines:
            if line.startswith("@"):
                row=line.split(" ")
                #TODO: This could possibly detect an error(seemingly they output two " ")
                #print(row)
                name=row[0].replace("@","")
                val=row[3].split("=")
                value=int(val[1])
                support_dict[name]=value
    #print(support_dict)
    rattle_dict = {}
    # open file in read mode
    with open(rattle_name) as f:
        lines = f.readlines()
    for line in lines[1:]:
        row = line.split(",")
        value=int(support_dict[row[0].replace(" ","")])
        if not row[1] in rattle_dict:
            rattle_dict[row[1].replace(" ","")] = value
        else:
            old_count=rattle_dict[row[1]]
            new_count=old_count+value
            up_dict={row[1]:new_count}
            rattle_dict.update(up_dict)
    print("RATTLE DICTUS",rattle_dict)
    return rattle_dict
def plot_results(id_list,original_list,rattle_list,isonform_list,gene_id,outfolder):
    # set width of bar
    barWidth = 0.25
    fig = plt.subplots(figsize=(12, 8))
    print("ISONF",isonform_list)
    #print("RATTLE LISTA",rattle_list)
    # set height of bar
    #IT = [12, 30, 1, 8, 22]
    #ECE = [28, 6, 16, 5, 10]
    #CSE = [29, 3, 24, 25, 17]

    # Set position of bar on X axis
    br1 = np.arange(len(original_list))
    br2 = [x + barWidth for x in br1]
    br3 = [x + barWidth for x in br2]

    # Make the plot
    plt.bar(br1, original_list, color='r', width=barWidth,
            edgecolor='grey', label='original')
    plt.bar(br2, rattle_list, color='g', width=barWidth,
            edgecolor='grey', label='rattle')
    plt.bar(br3, isonform_list, color='b', width=barWidth,
            edgecolor='grey', label='isONform')

    # Adding Xticks
    title="Found isoforms with abundances for gene id"+str(gene_id)
    plt.title(title)
    plt.xlabel('Sirv ID', fontweight='bold', fontsize=15)
    plt.ylabel('Read abundance', fontweight='bold', fontsize=15)
    plt.xticks([r + barWidth for r in range(len(original_list))],id_list)

    plt.legend()
    #plt.show()
    filename=str(gene_id)+'.png'
    plt.savefig(os.path.join(outfolder,filename))
def find_max_gene_id(original_dict):
    max=0
    for key in original_dict:
        tmp_key=key.replace("SIRV","")
        id=int(tmp_key)
        if id>max:
            max=id
    print(int(max/100))
    return int(max/100)
def get_gene_id(key):
    tmp_key = key.replace("SIRV", "")
    id = int(tmp_key)
    return int(id/100)

def unite_infos(original_dict,rattle_dict,isonform_dict,outfolder):
    print("ISONDICT",isonform_dict)
    max_gene_id=find_max_gene_id(original_dict)
    original_list = [[] for i in range(1,max_gene_id+1)]
    id_list=[[] for i in range(1,max_gene_id+1)]
    rattle_list=[[] for i in range(1,max_gene_id+1)]
    isonform_list=[[] for i in range(1,max_gene_id+1)]
    isonform_count=0
    #print(x)
    #id_list=[]
    #original_list=[]
    #rattle_list=[]
    for key,value in original_dict.items():
        gkey=get_gene_id(key)
        id_list[gkey-1].append(key)
        original_list[gkey-1].append(value)
        if key in rattle_dict:
            rattle_list[gkey-1].append(rattle_dict[key])
        else:
            rattle_list[gkey-1].append(0)
        if key in isonform_dict:
            isonform_list[gkey-1].append(isonform_dict[key])
            isonform_count+=isonform_dict[key]
        else:
            isonform_list[gkey-1].append(0)
    for i in range(max_gene_id):
        print("IDLIST",id_list[i])
        print("ORIGLIST",original_list[i])
        print("RATTLE LIST",rattle_list[i])
        plot_results(id_list[i],original_list[i],rattle_list[i],isonform_list[i],i+1,outfolder)
    print("Nr of reads support for Isonform:",isonform_count)
def main(args):
    isonform_dict=read_isonform_file(args.isonform_analysis)
    print("ISO DICT",isonform_dict)
    original_key_list=[]
    original_dict=read_original_analysis(args.original)
    for key in original_dict:
        original_key_list.append(key)
    rattle_dict=read_rattle_analysis(args.rattle_name,args.rattle_output)
    print(rattle_dict)
    outfolder=args.outfolder.replace("dummy.txt","")
    unite_infos(original_dict,rattle_dict,isonform_dict,outfolder)
    with open(args.outfolder, 'w') as f:
        f.write("dummy")
    #print(original_key_list)
    #print(csv_obj)
    #SIRV_dict=record_lines(csv_obj)
    #print(SIRV_dict)
    #write_Infos(SIRV_dict,args.outfile)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        "Parses a fasta file and output all sequences longer than min_size to fasta format to /path/to/fasta_file_filtered.fa.")
    parser.add_argument('original', type=str, help='CSV file. ')
    parser.add_argument('rattle_name', type=str, help='csv file. ')
    parser.add_argument('rattle_output', type=str, help='.fastq file. ')
    parser.add_argument('isonform_analysis',type=str, help='csv file. ')
    parser.add_argument('outfolder', type=str, help='csv file. ')
    #parser.add_argument('dummy_outfile', type=str, help='csv file. ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')

    args = parser.parse_args()

    main(args)
