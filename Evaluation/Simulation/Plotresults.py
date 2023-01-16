import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib import pyplot
def read_analysis(isonform_name):
    support_dict={}
    #print("RATTLE_OUTPUT",isonform_name)
    fp=0
    tp=0
    with open(isonform_name) as r_out:
        lines = r_out.readlines()
        for line in lines:
            if not line.startswith("q_acc"):
                row=line.split(",")
                #print(row)
                #TODO: This could possibly detect an error(seemingly they output two " ")
                #print(row)
                name_list=[]
                name_list.append(row[0])
                #sprint("NAMELIST",name_list)
                val=row[1].split("|")[-1]
                if val=="full":
                    val=1
                value=val
                if not value in support_dict:
                    support_dict[value]=name_list
                else:
                    fp+=1
                    new_names=support_dict[value]
                    new_names.extend(name_list)
                    new_entry={value: new_names}
                    support_dict.update(new_entry)
    tp=len(support_dict.keys())
    print(support_dict.keys())
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
    #print("ISONFORM DICTUS",isonform_dict)
    return isonform_dict,fp,tp
def plot_with_seaborn_precision(ticks,outfolder):
    #sns.set_theme()
    sns.set_style("whitegrid")
    palette = {
        'isONform': 'tab:red',
        'RATTLE': 'tab:green',
        # "strobealign_mixed" : 'magenta'
    }
    plt.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.6)
    sns.set_style("whitegrid")
    precision_path = os.path.join(args.outfolder,"precision.csv")
    recall_path = os.path.join(args.outfolder, "recall.csv")
    print(precision_path)
    #print(recall_path)
    prec_data = pd.read_csv(precision_path,sep=",")
    #recall_data=pd.read_csv(recall_path)
    print(prec_data)
    #print(recall_data)

    g = sns.violinplot(data=prec_data, x="nr_isos", y="precision", hue="tool", style="type",
                    #kind="line",  # dashes = dashes,
                    col="dataset"#, hue_order=tool,  # hue="datastructure", style="datastructure",
                    # col_wrap=3, col_order=["SIM1", "SIM2", "SIM4"], palette=palette)
                    #col_order=["SIM3"], palette=palette)
                       )
    #g.fig.suptitle('Precision')
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    g.set(ylim=(0,1.1), xlim=(0-0.5,len(ticks)),xlabel='Number of isoforms', ylabel='Precision')
    g.set_title("Precision SIM")
    #g.set_xticklabels(rotation=60, labels=[50, 75, 100, 150, 200, 250, 300, 500])
    #g.tight_layout()
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    plt.legend(title='Method')
    plt.gcf().subplots_adjust(bottom=0.15,left=0.20)
    plt.savefig(os.path.join(outfolder, "precision_SIM.pdf"))
    plt.close()

def plot_with_seaborn_recall(ticks, outfolder):
        sns.set_style("whitegrid")
        plt.rcParams.update({'font.size': 18})
        sns.set(font_scale=1.6)
        sns.set_style("whitegrid")
        precision_path = os.path.join(args.outfolder, "precision.csv")
        recall_path = os.path.join(args.outfolder, "recall.csv")
        #print(precision_path)
        print(recall_path)
        #prec_data = pd.read_csv(precision_path, sep=",")
        recall_data = pd.read_csv(recall_path)
        #print(prec_data)
        print(recall_data)

        g = sns.violinplot(data=recall_data, x="nr_isos", y="recall", hue="tool", style="type",
                           # kind="line",  # dashes = dashes,
                           col="dataset"  # , hue_order=tool,  # hue="datastructure", style="datastructure",
                           # col_wrap=3, col_order=["SIM1", "SIM2", "SIM4"], palette=palette)
                           # col_order=["SIM3"], palette=palette)
                           )
        # g.fig.suptitle('Precision')
        # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
        # axes = g.axes
        g.set_title("Recall SIM")
        g.set(ylim=(0, 1.1), xlim=(0-0.6, len(ticks)),xlabel='Number of isoforms', ylabel='Recall')

        # g.set_xticklabels(rotation=60, labels=[50, 75, 100, 150, 200, 250, 300, 500])
        # g.tight_layout()
        # g.set(ylim=(95, 100))
        # ax.set_xticks([18,24,30,36])
        plt.legend(title='Method')
        plt.gcf().subplots_adjust(bottom=0.15,left=0.20)
        plt.savefig(os.path.join(outfolder, "recall_SIM.pdf"))
        plt.close()
def plot_results(id_list,original_list,rattle_list,isonform_list,gene_id,outfolder):
    # set width of bar
    barWidth = 0.25
    fig = plt.subplots(figsize=(12, 8))
    #print("ISONF",isonform_list)
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
def write_to_csv(recall_list,precision_list,outfolder,filename):
    file_path=os.path.join(outfolder,filename)
    f = open(file_path, "w")
    f.write("recall\n")
    for i,rec in enumerate(recall_list):
        f.write(">{0}\n{1}\n".format(filename,rec))
    f.write("precision\n")
    for j,prec in enumerate(precision_list):
        f.write(">{0}\n{1}\n".format(filename,prec))
    f.close()
def calculate_statistics(tp,fp,fn):
    print("PARAMS",tp,fp,fn)
    if (tp+fp)>0:
        precision=tp/(tp+fp)
    else:
        precision=0
    if(tp+fn)>0:
        recall=tp/(tp+fn)
    else:
        recall=0
    print("P",precision," R",recall)
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
    isExist = os.path.exists(args.outfolder)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(args.outfolder)
    recall_path = os.path.join(args.outfolder, "recall.csv")
    recall_file = open(recall_path, "w")
    # print("RATTLEPATHS",rattle_path)
    precision_path = os.path.join(args.outfolder, "precision.csv")
    precision_file = open(precision_path, "w")
    #print("SIONPATHS", ison_path)
    max_iso_nr=args.max_isoforms-1
    #we save the precision and recall values for isONform an Rattle in lists
    ison_precision_list=[None] * max_iso_nr
    ison_recall_list=[None] * max_iso_nr
    rattle_precision_list=[None] * max_iso_nr
    rattle_recall_list=[None] * max_iso_nr
    print("IN",args.indir)
    print("OUT",args.outfolder)
    precision_file.write("id,nr_isos,tool,precision\n")
    recall_file.write("id,nr_isos,tool,recall\n")

    for analysisfile in os.listdir(args.indir):
        if analysisfile.startswith("results_ison_"):
            print("File",analysisfile)
            file_dir = os.path.join(args.indir, analysisfile)
            isonform_dict,fp,tp=read_analysis(file_dir)
            #tp=len(isonform_dict)
            filename=analysisfile.split(".")[0]
            nr_isos=int(filename.split("_")[2])
            run_id=int(filename.split("_")[3])
            id = str(nr_isos) + "_" + str(run_id)
            #print(nr_isos)
            #print(run_id)
            fn=max(nr_isos-tp,0)
            precision,recall=calculate_statistics(tp,fp,fn)
            print("Precision", precision)
            print("Position", nr_isos - 2)
            print("Recall",recall)
            precision_file.write("{0},{1},{2},{3}\n".format(id, nr_isos, "isONform", precision))
            # rattle_file.write("precision\n")
            recall_file.write("{0},{1},{2},{3}\n".format(id, nr_isos, "isONform", recall))

            #print(ison_precision_list)
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
            ison_name = "ison_res_" + str(nr_isos) + "_" + str(run_id) + ".csv"

        elif analysisfile.startswith("results_rattle_"):
            print("File", analysisfile)

            file_dir = os.path.join(args.indir, analysisfile)
            rattle_dict, fp,tp = read_analysis(file_dir)
            #tp = len(rattle_dict)
            filename = analysisfile.split(".")[0]
            nr_isos = int(filename.split("_")[2])
            run_id = int(filename.split("_")[3])
            id = str(nr_isos) + "_" + str(run_id)
            #print(nr_isos)
            #print(run_id)
            fn = nr_isos - tp
            precision, recall = calculate_statistics(tp, fp, fn)
            print("Precision", precision)
            print("Position", nr_isos - 2)
            print("Recall", recall)
            precision_file.write("{0}, {1}, {2}, {3}\n".format(id, nr_isos, "RATTLE", precision))
            recall_file.write("{0}, {1}, {2}, {3}\n".format(id, nr_isos, "RATTLE", recall))

            #print("rattle_precision_list",rattle_precision_list)
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
            rattle_name="rattle_res_"+str(nr_isos)+"_"+str(run_id)+".csv"

    #write_to_csv(ison_precision_list, ison_recall_list, args.outfolder, ison_name,id)
    #write_to_csv(rattle_precision_list, rattle_recall_list, args.outfolder, rattle_name,id)
    #print("ison_prec_dict", ison_precision_list, ", ison_rec_dict", ison_recall_list)
    #print("rattle_prec_dict",rattle_precision_list,", rattle_rec_dict", rattle_recall_list)
    """plt.rcParams["figure.figsize"] = [7.50, 3.50]
    plt.rcParams["figure.autolayout"] = True
    
    fig, ax = plt.subplots()

    ax.boxplot(ison_precision_dict.values())
    ax.set_xticklabels(ison_precision_dict.keys())
"""
    precision_file.close()
    recall_file.close()
    #write_to_csv(ison_precision_list,ison_recall_list,args.outfolder,id)
    ticks = ['2', '3', '4','5','6','7','8','9','10','11','12','13','14','15']
    plot_with_seaborn_precision(ticks, args.outfolder)
    plot_with_seaborn_recall(ticks, args.outfolder)
    #plot_data("Precision", ison_precision_list, rattle_precision_list, ticks,args.outfolder)
    #plot_data("Recall", ison_recall_list, rattle_recall_list, ticks,args.outfolder)
    #print(original_key_list)
    #print(csv_obj)
    #SIRV_dict=record_lines(csv_obj)
    #print(SIRV_dict)
    #write_Infos(SIRV_dict,args.outfile)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        "Generates Graphs from the alignments.")
    #parser.add_argument('original', type=str, help='CSV file. ')
    #parser.add_argument('rattle_name', type=str, help='csv file. ')
    #parser.add_argument('rattle_output', type=str, help='.fastq file. ')
    parser.add_argument('--indir',type=str, help='csv file. ')
    parser.add_argument('--outfolder', type=str, help='csv file. ')
    parser.add_argument('--max_isoforms',type=int, help="The maximum number of isoforms we generate in our files")
    #parser.add_argument('dummy_outfile', type=str, help='csv file. ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')

    args = parser.parse_args()

    main(args)
