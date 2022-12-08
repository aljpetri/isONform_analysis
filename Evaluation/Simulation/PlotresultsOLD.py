import argparse
import os
import matplotlib.pyplot as plt
import numpy as np

def read_analysis(isonform_name):
    support_dict={}
    print("RATTLE_OUTPUT",isonform_name)
    fp=0
    with open(isonform_name) as r_out:
        lines = r_out.readlines()
        for line in lines:
            if not line.startswith("q_acc"):
                row=line.split(",")
                print(row)
                #TODO: This could possibly detect an error(seemingly they output two " ")
                #print(row)
                name_list=[]
                name_list.append(row[0])
                print(name_list)
                val=row[1].split("|")[-1]
                if val=="full":
                    val=1
                value=int(val)
                if not value in support_dict:
                    support_dict[value]=name_list
                else:
                    fp+=1
                    new_names=support_dict[value]
                    new_names.extend(name_list)
                    new_entry={value: new_names}
                    support_dict.update(new_entry)
    print(support_dict)
    isonform_dict = {}
    # open file in read mode
    with open(isonform_name) as f:
        lines = f.readlines()
    for line in lines[1:]:
        row = line.split(",")

        if not row[1] in isonform_dict:
            isonform_dict[row[1].replace(" ","")] = value
        else:
            old_count=isonform_dict[row[1]]
            new_count=old_count+value
            up_dict={row[1]:new_count}
            isonform_dict.update(up_dict)
    print("ISONFORM DICTUS",isonform_dict)
    return isonform_dict,fp

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


#def generate_isonform_statistics(curr_file):
def set_box_color(bp, color):
        plt.setp(bp['boxes'], color=color)
        plt.setp(bp['whiskers'], color=color)
        plt.setp(bp['caps'], color=color)
        plt.setp(bp['medians'], color=color)
def write_to_csv(isonform_dict,rattle_dict,outfolder):
    file_path=os.path.join(outfolder,"data.csv")
    f = open(file_path, "w")
    f.write("IsONform")
    for key,value in isonform_dict.items():
        f.write(">{0}\n{1}\n".format(key,value))
    f.write("Rattle")
    for key,value in rattle_dict.items():
        f.write(">{0}\n{1}\n".format(key,value))
    f.close()
def calculate_statistics(tp,fp,fn):
    if (tp+fp)>0:
        precision=tp/(tp+fp)
    else:
        precision=0
    if(tp+fn)>0:
        recall=tp/(tp+fn)
    else:
        recall=0
    return precision,recall
def plot_data(header,data1,data2,ticks,outfolder):
    plt.figure()

    bpl = plt.boxplot(data1, positions=np.array(range(len(data1))) * 2.0 - 0.4, widths=0.6)
    bpr = plt.boxplot(data2, positions=np.array(range(len(data2))) * 2.0 + 0.4, widths=0.6)
    set_box_color(bpl, '#D7191C')  # colors are from http://colorbrewer2.org/
    set_box_color(bpr, '#2C7BB6')

    # draw temporary red and blue lines and use them to create a legend
    plt.plot([], c='#D7191C', label='isONform')
    plt.plot([], c='#2C7BB6', label='RATTLE')
    plt.legend()
    plt.title(header)
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-2, len(ticks) * 2)
    plt.ylim(0, 1.1)
    plt.tight_layout()
    file=header+".png"
    filename=os.path.join(outfolder, file)
    plt.savefig(filename)
    plt.show()
def main(args):
    #generate_isonform_statistics(curr_file)
    max_int=15
    #cl_dir = os.path.join(outdir, cl_id)
    ison_precision_list=[None] * (max_int-1)
    ison_recall_list=[None] * (max_int-1)
    rattle_precision_list=[None] * (max_int-1)
    rattle_recall_list=[None] * (max_int-1)
    print("IN",args.indir)
    print("OUT",args.outfolder)
    isExist = os.path.exists(args.outfolder)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(args.outfolder)
    for analysisfile in os.listdir(args.indir):
        if analysisfile.startswith("results_ison_"):
            print("File",analysisfile)
            file_dir = os.path.join(args.indir, analysisfile)
            isonform_dict,fp=read_analysis(file_dir)
            tp=len(isonform_dict)
            print("TRUE POSITIVE",tp)

            filename=analysisfile.split(".")[0]
            nr_isos=int(filename.split("_")[2])
            run_id=int(filename.split("_")[3])
            print(nr_isos)
            print(run_id)
            fn=nr_isos-tp
            precision,recall=calculate_statistics(tp,fp,fn)
            print (ison_precision_list)
            if ison_precision_list[nr_isos-2]:
                #old_prec_list=ison_precision_dict[nr_isos]
                #old_prec_list.append(precision)
                #prec_repl_dict={nr_isos: old_prec_list}
                #ison_precision_dict.update(prec_repl_dict)
                #old_rec_list = ison_recall_dict[nr_isos]
                #old_rec_list.append(recall)
                #rec_repl_dict = {nr_isos: old_rec_list}
                #ison_recall_dict.update(rec_repl_dict)
                ison_precision_list[nr_isos-2].append(precision)
                ison_recall_list[nr_isos-2].append(recall)
            else:
                #list_prec=[]
                #list_prec.append(precision)
                #ison_precision_dict[nr_isos]=list_prec
                #list_rec = []
                #list_rec.append(recall)
                #ison_recall_dict[nr_isos] = list_rec
                ison_precision_list[nr_isos-2] = []
                ison_precision_list[nr_isos-2].append(precision)
                ison_recall_list[nr_isos - 2] = []
                ison_recall_list[nr_isos - 2].append(recall)
        elif analysisfile.startswith("results_rattle_"):
            print("File", analysisfile)
            file_dir = os.path.join(args.indir, analysisfile)
            rattle_dict, fp = read_analysis(file_dir)
            tp = len(rattle_dict)

            filename = analysisfile.split(".")[0]
            nr_isos = int(filename.split("_")[2])
            run_id = int(filename.split("_")[3])
            print(nr_isos)
            print(run_id)
            fn = nr_isos - tp
            precision, recall = calculate_statistics(tp, fp, fn)
            if rattle_precision_list[nr_isos-2]:
                #old_prec_list = rattle_precision_dict[nr_isos]
                #old_prec_list.append(precision)
                #prec_repl_dict = {nr_isos: old_prec_list}
                #rattle_precision_dict.update(prec_repl_dict)
                #old_rec_list = rattle_recall_dict[nr_isos]
                #old_rec_list.append(recall)
                #rec_repl_dict = {nr_isos: old_rec_list}
                #rattle_recall_dict.update(rec_repl_dict)
                rattle_precision_list[nr_isos-2].append(precision)
                rattle_recall_list[nr_isos - 2].append(recall)
            else:
                #list_prec = []
                #list_prec.append(precision)
                #rattle_precision_dict[nr_isos] = list_prec
                #list_rec = []
                #list_rec.append(recall)
                #rattle_recall_dict[nr_isos] = list_rec
                rattle_precision_list[nr_isos-2]=[]
                rattle_precision_list[nr_isos-2].append(precision)
                rattle_recall_list[nr_isos - 2] = []
                rattle_recall_list[nr_isos - 2].append(recall)
    print("ison_prec_dict", ison_precision_list, ", ison_rec_dict", ison_recall_list)
    print("rattle_prec_dict",rattle_precision_list,", rattle_rec_dict", rattle_recall_list)
    """plt.rcParams["figure.figsize"] = [7.50, 3.50]
    plt.rcParams["figure.autolayout"] = True
    
    fig, ax = plt.subplots()

    ax.boxplot(ison_precision_dict.values())
    ax.set_xticklabels(ison_precision_dict.keys())
"""
    write_to_csv(isonform_dict,rattle_dict,args.outfolder)
    ticks = ['2', '3', '4','5','6','7','8','9','10','11','12','13','14','15']
    plot_data("Precision", ison_precision_list, rattle_precision_list, ticks,args.outfolder)
    plot_data("Recall", ison_recall_list, rattle_recall_list, ticks,args.outfolder)
    #print(original_key_list)
    #print(csv_obj)
    #SIRV_dict=record_lines(csv_obj)
    #print(SIRV_dict)
    #write_Infos(SIRV_dict,args.outfile)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        "Parses a fasta file and output all sequences longer than min_size to fasta format to /path/to/fasta_file_filtered.fa.")
    #parser.add_argument('original', type=str, help='CSV file. ')
    #parser.add_argument('rattle_name', type=str, help='csv file. ')
    #parser.add_argument('rattle_output', type=str, help='.fastq file. ')
    parser.add_argument('--indir',type=str, help='csv file. ')
    parser.add_argument('--outfolder', type=str, help='csv file. ')
    #parser.add_argument('dummy_outfile', type=str, help='csv file. ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')

    args = parser.parse_args()

    main(args)