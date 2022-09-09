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
        rows.append(row)
    return rows
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
        f.write("Hello World")
        for key,value in SIRV_dict.items():
            #print(type(key))
            f.write(str(key+","+str(len(value))+"\n"))
def main():
    outfile="/home/alexanderpetri/isONform_analysis/Analysis/SIRVtest/Final.txt"

    for line in fileinput.input():
        inputname=fileinput.filename()
        print("File name is: ", inputname)
        fileinput.nextfile()

    #read_all_batchfiles(args.folder)
    #csv_obj=read_csv(args.analysis)
    #print(csv_obj)
    #SIRV_dict=record_lines(csv_obj)
    #print(SIRV_dict)
    SIRV_dict={}
    write_Infos(SIRV_dict,outfile)
if __name__ == '__main__':
    #parser.add_argument('outfile', type=str, help='csv file. ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')


    main()
