from csv import reader
import argparse
import sys
import os
import csv
import fileinput


"""def read_all_batchfiles(folder):
    for batchfile in os.listdir(folder):
        print(batchfile)
        if batchfile.endswith("_batch.csv"):"""

def read_csv(filename):
    # open file in read mode
    rows=[]
    # open file in read mode
    with open(filename) as f:
        lines = f.readlines()
    for line in lines:
        row=line.split(",")
        #print(row)
        #print(type(row[1]),",",row[1])
        if not row[0]=='read':
            rows.append(row)
    SIRV_dict=parse_rows(rows)
    return SIRV_dict
def parse_rows(rows):
    SIRV_dict={}
    for row in rows:
        id=row[1].replace(" ","")
        if id in SIRV_dict:
            val=SIRV_dict[id]
            val+=1
            SIRV_dict[id]=val
        else:
            SIRV_dict[id]=1
    return SIRV_dict
def record_lines(csv_obj):
    SIRV_dict={}
    for line in csv_obj:
        if str(line[1])!="ref":
            #print(line[1])
            if line[1] in SIRV_dict:
                prev_list=SIRV_dict[line[1]]
                prev_list.append(line[0])
                SIRV_dict[line[1]]=prev_list
            else:
                lista=[]
                lista.append(line[0])
                SIRV_dict[line[1]]=lista
    return SIRV_dict
def write_Infos(SIRV_dict,outfile):
    with open(outfile, 'w') as f:
        print("Writing")
        for key,value in SIRV_dict.items():
            #print(type(key))
            f.write(str(key+","+str(value)+"\n"))
        print("Writing done")
def main(args):
    #outfile="/home/alexanderpetri/isONform_analysis/Analysis/SIRVtest/Final.txt"
    #for line in fileinput.input():
        #print("File name is: ", fileinput.filename())
        #with open(fileinput.filename()) as f:
         #   lines = f.readlines()
        #fileinput.nextfile()
    print("Hello World")
    #outfile=args.outfile
    #infile=args.infile
    SIRV_dict=read_csv(args.infile)
    print(SIRV_dict)
    #with open(args.outfile, 'w') as f:
    #    f.write("Hello")
    #print("Hello World")
    write_Infos(SIRV_dict,args.outfile)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="De novo error correction of long-read transcriptome reads",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile',type=str,help='input_csv file.')
    parser.add_argument('outfile', type=str, help='output_csv file. ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')
    args = parser.parse_args()

    main(args)