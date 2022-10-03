import argparse
import os
import itertools
import sys
from pathlib import Path
def read_counting_file(counting_file,alignment_res):
    path=str(counting_file)
    infoslist = path.split("/")[-1]
    id_list=infoslist.split("_")
    isoform_number=id_list[1]
    run_id=id_list[2]
    print(infoslist)
    file1 = open(counting_file, 'r')
    lines = file1.readlines()
    print(counting_file)
    for line in lines:
        if not line.startswith("read"):
            line_list=line.split(",")
            this_id=line_list[0]
            mapped_id=line_list[1]
            print("Line",this_id,",",mapped_id)
def read_mapping_file(mapping_name,batch_mappings_id):
    #mappingname =  "mapping"+str(batch_id)+".txt"
    incoming_ctr=0
    #print("MAPPING",batch_id,mappingname)
    with open(mapping_name) as g:
        for id, reads in itertools.zip_longest(*[g] * 2):
            #print(id, reads)
            inter_id = id.replace('\n', '')
            #id = int(inter_id.replace('consensus', ''))
            # print("ID",id)
            #print(reads)
            reads = reads.replace('\n', '')
            reads = reads.replace('[', '')
            reads = reads.replace(']', '')
            # reads=reads.replace('\'','')
            reads = reads.replace("'", "")
            readslist = reads.split(", ")
            #print(readslist)
            # print(len(readslist))
            batch_mappings_id[id] = readslist

            incoming_ctr+=len(readslist)
            #print(batch_mappings_id)
    #if batch_id == 5:
        #print(batch_mappings_id)
    print("INCOMING",incoming_ctr, ", ",mapping_name)
    return batch_mappings_id
def calculate_tp_fn_fp(id_dict,infos_dict):
    tp=0
    fn=0
    fp=0
    tp_entries={}
    fn_entries={}
    print(id_dict)
    for id, nr_reads in id_dict.items():
        if len(id_dict)<2:
            if not id in infos_dict:
                infos_dict[id]=nr_reads
            elif nr_reads>infos_dict[id]:
                up_dict={id:nr_reads}
                fn+=infos_dict[id]
                infos_dict.update(up_dict)
        else:
            max_elem=0
            max_id=None
            for id,nr_reads in id_dict.items():
                print(id,nr_reads)
                if nr_reads>max_elem:
                    max_elem=nr_reads
                    max_id=id
            for id,nr_reads in id_dict.items():
                if id != max_id:
                    fp+=nr_reads
    for id,nr_reads in infos_dict.items():
        tp+=nr_reads
    return tp,fn,fp
def analyse_mappings(batch_mappings_id,mapping_infos_dict):
    infos_dict={}
    for cons_id,supp_reads in batch_mappings_id.items():
        id_dict={}
        for supp_read in supp_reads:
            id=supp_read.split("|")[2]
            print(id)
            ison_ident=id.split(".")[0]
            if not ison_ident in id_dict:
                id_dict[ison_ident]=1
            else:
                id_count=id_dict[ison_ident]
                id_count+=1
                repl_dict={ison_ident:id_count}
                id_dict.update(repl_dict)
        tp,fn,fp = calculate_tp_fn_fp(batch_mappings_id,infos_dict)
        mapping_infos_dict[id]=id_dict
        return tp,fn,fp
def main(args):
    print("Evaluating performance on Simulated Reads")
    alignment_res={}
    batch_mappings_id={}
    mapping_infos_dict={}
    directory = args.fastq_folder #os.fsencode(args.fastq_folder)
    pat = Path(directory)
    file_list = list(pat.rglob('counts_*.csv'))
    incoming_cter_glob=0
    for file in file_list:
        tp_all = 0
        fp_all = 0
        fn_all = 0
        if not file=="/home/alexanderpetri/isONform_analysis/Para_out_500_September/transcriptome_mapping.txt":
            #read_mapping_file(file,batch_mappings_id)
            read_counting_file(file,alignment_res)
            """tp,fn,fp=analyse_mappings(batch_mappings_id,mapping_infos_dict)
            tp_all+=tp
            fn_all+=fn
            fp_all+=fp
        print("File:",file)
        print("TP:",tp_all,", FN:",fn_all,", FP: ",fp_all)
        precision=tp_all/(tp_all+fp_all)
        print("Precision: ",precision)
        recall=tp_all/(tp_all+fn_all)
        print("Recall:",recall)"""
    #other_file_list = list(pat.rglob('*_mapping_low_abundance.txt'))
    #incoming_cter_glob_low = 0
    #for o_file in other_file_list:
    #    incoming_cter_glob_low+=read_mapping_file(o_file)
    #print("LOW ABUNDANCE",incoming_cter_glob_low)
    #sum=incoming_cter_glob+incoming_cter_glob_low
    #print("SUM",sum)

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description="De novo clustering of long-read transcriptome reads",
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('--fastq_folder', type=str, default=False,
                            help='Path to input fastq folder with reads in clusters')
        args = parser.parse_args()
        main(args)