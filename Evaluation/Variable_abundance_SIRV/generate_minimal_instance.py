import argparse
import os
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

def main(args):

    clust_count=0
    other_count=0
    all_reads = {acc: ( seq, qual) for i, (acc, (seq, qual)) in
                 enumerate(readfq(open(args.fastq, 'r')))}
    #print(all_reads)
    id_list1=['104:801|ERR3588903.874408_63195041-90e9-429c-88d3-ca42bba308d0/1_strand=-', '104:764|ERR3588903.1153030_e6a68ae6-d80c-4f09-af98-eeeee2965576/1_strand=-', '106:760|ERR3588903.1506361_57e65bfa-75a5-415c-abc9-949c9e5bcfa0/1_strand=-', '106:740|ERR3588903.1570652_637b4631-60e1-411f-bee9-180efd487ce0/1_strand=-', '103:699|ERR3588903.611008_e6c612b8-86d7-408b-bfb4-9660d34eac7f/1_strand=+', '111:687|ERR3588903.411902_56b6af38-7ebb-40e1-925f-c183faaf4298/1_strand=-', '102:696|ERR3588903.1274621_e9953f4e-cd7a-4771-8938-d845ded72e93/1_strand=+']
    #id_list2=['112:1090|ERR3588903.1649107', '118:1062|ERR3588903.1106035', '105:1047|ERR3588903.871778', '102:1060|ERR3588903.252682', '124:1075|ERR3588903.207875', '106:1045|ERR3588903.534482', '109:1049|ERR3588903.1338655', '110:1059|ERR3588903.741123', '110:1071|ERR3588903.1654493', '112:1055|ERR3588903.1062833', '104:1030|ERR3588903.629950', '110:1062|ERR3588903.642184', '107:1050|ERR3588903.842598', '108:1053|ERR3588903.1625799', '114:1058|ERR3588903.268645', '117:1050|ERR3588903.1088582', '109:1050|ERR3588903.609091', '103:1048|ERR3588903.1119890', '103:1057|ERR3588903.1444517', '105:1059|ERR3588903.1029159', '115:1042|ERR3588903.496246', '101:1042|ERR3588903.977106', '117:1058|ERR3588903.1320610', '109:1028|ERR3588903.731657', '108:1044|ERR3588903.1158992', '103:1063|ERR3588903.121639', '113:1073|ERR3588903.474275', '107:1059|ERR3588903.979891', '108:1048|ERR3588903.1311752', '114:1034|ERR3588903.522653', '115:1052|ERR3588903.531211', '109:1039|ERR3588903.1557135', '112:1061|ERR3588903.1108528', '113:1075|ERR3588903.846587', '105:1082|ERR3588903.1432681', '118:1077|ERR3588903.621522', '106:1033|ERR3588903.1300304', '92:1050|ERR3588903.894889', '100:1063|ERR3588903.1228511', '100:1041|ERR3588903.299087', '105:1052|ERR3588903.585817', '111:1062|ERR3588903.180268', '109:1051|ERR3588903.1134494', '111:1051|ERR3588903.791401', '98:1016|ERR3588903.1494709', '120:1073|ERR3588903.699562', '106:1062|ERR3588903.994913', '98:1050|ERR3588903.257193', '780:1736|ERR3588903.34911', '102:1063|ERR3588903.1374257', '106:1067|ERR3588903.1658063', '113:1042|ERR3588903.922085', '111:1071|ERR3588903.738063', '88:1025|ERR3588903.1055655', '106:1027|ERR3588903.330514', '111:1007|ERR3588903.81656', '103:1029|ERR3588903.1488144', '104:1020|ERR3588903.1360487', '107:1036|ERR3588903.1289463', '115:1059|ERR3588903.305370', '116:1054|ERR3588903.1604964', '97:1015|ERR3588903.686170', '111:1041|ERR3588903.527175', '116:1042|ERR3588903.192817', '105:1077|ERR3588903.1275790', '105:1038|ERR3588903.55691', '110:1063|ERR3588903.1645377', '105:1041|ERR3588903.1351983', '112:1077|ERR3588903.1126409', '110:1061|ERR3588903.1624215', '102:1055|ERR3588903.738296', '101:1044|ERR3588903.1008006', '108:1016|ERR3588903.852773', '107:1056|ERR3588903.1014927', '117:1038|ERR3588903.1637298', '120:1077|ERR3588903.389537', '106:1063|ERR3588903.916120', '110:1049|ERR3588903.330683', '103:1055|ERR3588903.686679', '103:1037|ERR3588903.518049', '104:1026|ERR3588903.1010735', '102:1030|ERR3588903.1350253', '110:1048|ERR3588903.491855', '118:1070|ERR3588903.395866', '112:1063|ERR3588903.1190832', '109:1054|ERR3588903.141821', '114:1064|ERR3588903.363859', '115:1058|ERR3588903.120508', '129:1074|ERR3588903.951025', '118:1054|ERR3588903.1061564', '114:1054|ERR3588903.339153', '117:1057|ERR3588903.1040486', '106:1040|ERR3588903.990631', '100:1042|ERR3588903.529026', '106:1047|ERR3588903.74953', '100:1054|ERR3588903.135055', '119:1051|ERR3588903.143462', '117:1064|ERR3588903.762660', '112:1067|ERR3588903.536882', '100:1028|ERR3588903.498022', '111:1040|ERR3588903.651684', '108:1039|ERR3588903.1422773', '112:1047|ERR3588903.909811', '113:1060|ERR3588903.1257699', '104:1056|ERR3588903.399830', '111:1055|ERR3588903.573940', '110:1044|ERR3588903.454871', '102:1063|ERR3588903.1024228', '105:1050|ERR3588903.77246', '103:1040|ERR3588903.1435541', '113:1056|ERR3588903.1671679', '102:1037|ERR3588903.791563', '122:1056|ERR3588903.1418833', '104:1050|ERR3588903.51957', '128:1069|ERR3588903.751145', '109:1034|ERR3588903.858275', '104:1045|ERR3588903.571525', '108:1080|ERR3588903.1530725', '104:1027|ERR3588903.491906', '105:1039|ERR3588903.47229', '110:1037|ERR3588903.1150980', '106:967|ERR3588903.479313', '100:1015|ERR3588903.522444', '111:1022|ERR3588903.1526361', '122:1045|ERR3588903.15072']
    id_list2=[]
    #id_list1=[]
    new_reads={}
    for id in id_list1:
        new_reads[id]=all_reads[id]
    for id2 in id_list2:
        new_reads[id2]=all_reads[id2]
    outfile = open(os.path.join(args.outfolder, "1.fastq"), "w")
    for acc, ( seq, qual) in new_reads.items():
        outfile.write("@{0}\n{1}\n+\n{2}\n".format(acc, seq, qual))
    outfile.close()
    #print("isONclust:",clust_count)
    #print("OTHER_count",other_count)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        "Parses a fasta file and output all sequences longer than min_size to fasta format to /path/to/fasta_file_filtered.fa.")
    #parser.add_argument('original', type=str, help='CSV file. ')
    #parser.add_argument('rattle_name', type=str, help='csv file. ')
    #parser.add_argument('rattle_output', type=str, help='.fastq file. ')
    parser.add_argument('--fastq',type=str, help='fastq file. ')
    parser.add_argument('--outfolder', type=str, help='folder to write the file into ')
    #parser.add_argument('--max_isoforms',type=int, help="The maximum number of isoforms we generate in our files")
    #parser.add_argument('dummy_outfile', type=str, help='csv file. ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')

    args = parser.parse_args()

    main(args)
