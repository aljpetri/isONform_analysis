import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib import pyplot
def Average(lst):
    sum=0
    for ele in lst:
        sum += ele
    return sum / len(lst)
def read_analysis(isonform_name):
    support_dict = {}
    # print("RATTLE_OUTPUT",isonform_name)
    fp = 0
    tp = 0
    perc_identity_list=[]
    reconstr=[]
    with open(isonform_name) as r_out:
        lines = r_out.readlines()
        for line in lines:
            if not line.startswith("q_acc"):
                row = line.split(",")
                # print(row)
                # print(row)
                name_list = []
                name_list.append(row[0])
                # sprint("NAMELIST",name_list)
                perc_identity_list.append(float(row[9]))
                reconstr.append(float(row[10]))
    #nr_errors=sum(error_nr_list)
    avg_perc_identity=Average(perc_identity_list)
    # print("ISONFORM DICTUS",isonform_dict)
    perfect_reconstrs=len([1 for n1 in reconstr if n1 >95 and n1 < 105])
    partial=len(reconstr)-perfect_reconstrs
    return avg_perc_identity,perfect_reconstrs,partial

def plot_with_seaborn_total_errors(ticks,outfolder):
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
    total_error_path = os.path.join(args.outfolder, "Total_Errors.csv")
    print(total_error_path)
    total_errors_data = pd.read_csv(total_error_path,sep=",")
    #print(prec_data)

    g = sns.swarmplot(data=total_errors_data, x="nr_isos", y="total_errors", hue="tool", style="type",
                    #kind="line",  # dashes = dashes,
                    col="dataset"#, hue_order=tool,  # hue="datastructure", style="datastructure",
                    # col_wrap=3, col_order=["SIM1", "SIM2", "SIM4"], palette=palette)
                    #col_order=["SIM3"], palette=palette)
                       )
    #g.fig.suptitle('Precision')
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    maximum=max(total_errors_data.loc[:,"total_errors"])
    print("MAXI",maximum)
    g.set(ylim=(0,maximum+1), xlim=(0-0.5,len(ticks)))
    g.set_title("Total Errors SIRV")
    #g.set_xticklabels(rotation=60, labels=[50, 75, 100, 150, 200, 250, 300, 500])
    #g.tight_layout()
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])

    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig(os.path.join(outfolder, "Total_Errors.pdf"))
    plt.close()
def plot_with_seaborn_reconstruction(ticks,outfolder):
    #sns.set_theme()
    #sns.set_style("whitegrid")
    palette = {
        'isONform': 'tab:red',
        'RATTLE': 'tab:green',
        # "strobealign_mixed" : 'magenta'
    }
    plt.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.6)
    #sns.set_style("whitegrid")
    precision_path = os.path.join(args.outfolder,"Reconstruction.csv")
    #recall_path = os.path.join(args.outfolder, "recall.csv")
    print(precision_path)
    #print(recall_path)
    prec_data = pd.read_csv(precision_path,sep=",")
    #recall_data=pd.read_csv(recall_path)
    print(prec_data)
    #print(recall_data)
    #sns.histplot(data=tips, x="day", hue="sex", multiple="dodge", shrink=.8)
    g = sns.histplot(data=prec_data, x="nr_isos", y="full", hue="tool"#,multiple="stack"#, style="type",
                    #kind="line",  # dashes = dashes,
                    #col="dataset"#, hue_order=tool,  # hue="datastructure", style="datastructure",
                    # col_wrap=3, col_order=["SIM1", "SIM2", "SIM4"], palette=palette)
                    #col_order=["SIM3"], palette=palette)
                       )
    #g.fig.suptitle('Precision')
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    #maximum = max(prec_data.loc[:, "full"])
    #print("MAXI", maximum)
    #g.set(ylim=(0, maximum + 1), xlim=(0-0.5, len(ticks)),xlabel='Number of isoforms', ylabel='Average error rate (%)')
    g.set_title("Average Error Rate SIRV")
    #g.set_xticklabels(rotation=60, labels=[50, 75, 100, 150, 200, 250, 300, 500])
    #g.tight_layout()
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    plt.legend(title='Method')
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig(os.path.join(outfolder, "Perc_reconstruction.pdf"))
    plt.close()
def plot_with_seaborn_perc_id(ticks, outfolder):
        #sns.set_style("whitegrid")
        plt.rcParams.update({'font.size': 18})
        sns.set(font_scale=1.6)
        #sns.set_style("whitegrid")
        #precision_path = os.path.join(args.outfolder, "precision.csv")
        error_rates_path = os.path.join(args.outfolder, "Perc_identity.csv")
        print(error_rates_path)
        error_rates_data = pd.read_csv(error_rates_path, sep=",")
        print(error_rates_data)

        g = sns.boxplot(data=error_rates_data, x="nr_isos", y="perc_identity", hue="tool"#, style="type",
                           # kind="line",  # dashes = dashes,
                           #col="dataset"  # , hue_order=tool,  # hue="datastructure", style="datastructure",
                           # col_wrap=3, col_order=["SIM1", "SIM2", "SIM4"], palette=palette)
                           # col_order=["SIM3"], palette=palette)
                           )
        # g.fig.suptitle('Precision')
        # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
        # axes = g.axes
        g.set_title("Percent identity")
        #maximum = max(error_rates_data.loc[:, "error_rates"])
        #print("MAXI", maximum)
        g.set(ylim=(85, 101), xlim=(0, len(ticks)),xlabel='Number of isoforms', ylabel='Identity (%)')

        # g.set_xticklabels(rotation=60, labels=[50, 75, 100, 150, 200, 250, 300, 500])
        # g.tight_layout()
        # g.set(ylim=(95, 100))
        # ax.set_xticks([18,24,30,36])
        plt.legend(title='Method')
        plt.gcf().subplots_adjust(bottom=0.15,left=0.20)
        plt.savefig(os.path.join(outfolder, "Perc_identity.pdf"))
        plt.close()
def plot_results(id_list, original_list, rattle_list, isonform_list, gene_id, outfolder):
    # set width of bar
    barWidth = 0.25
    fig = plt.subplots(figsize=(12, 8))
    # print("ISONF",isonform_list)
    # print("RATTLE LISTA",rattle_list)
    # set height of bar
    # IT = [12, 30, 1, 8, 22]
    # ECE = [28, 6, 16, 5, 10]
    # CSE = [29, 3, 24, 25, 17]

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
    title = "Found isoforms with abundances for gene id" + str(gene_id)
    plt.title(title)
    plt.xlabel('Sirv ID', fontweight='bold', fontsize=15)
    plt.ylabel('Read abundance', fontweight='bold', fontsize=15)
    plt.xticks([r + barWidth for r in range(len(original_list))], id_list)

    plt.legend()
    # plt.show()
    filename = str(gene_id) + '.png'
    plt.savefig(os.path.join(outfolder, filename))


# def generate_isonform_statistics(curr_file):
def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)


def write_to_csv(recall_list, precision_list, outfolder, filename):
    file_path = os.path.join(outfolder, filename)
    f = open(file_path, "w")
    f.write("recall\n")
    for i, rec in enumerate(recall_list):
        f.write(">{0}\n{1}\n".format(filename, rec))
    f.write("precision\n")
    for j, prec in enumerate(precision_list):
        f.write(">{0}\n{1}\n".format(filename, prec))
    f.close()


def calculate_statistics(tp, fp, fn):
    print("PARAMS", tp, fp, fn)
    if (tp + fp) > 0:
        precision = tp / (tp + fp)
    else:
        precision = 0
    if (tp + fn) > 0:
        recall = tp / (tp + fn)
    else:
        recall = 0
    print("P", precision, " R", recall)
    return precision, recall

def max_value(inputlist):
    return max([max(sublist) for sublist in inputlist])
def plot_data(header, data1, data2, ticks, outfolder):
    plt.figure()

    bpl = plt.boxplot(data1, positions=np.array(range(len(data1))) * 2.0 - 0.4, widths=0.6)
    bpr = plt.boxplot(data2, positions=np.array(range(len(data2))) * 2.0 + 0.4, widths=0.6)
    set_box_color(bpl, '#D7191C')  # colors are from http://colorbrewer2.org/
    set_box_color(bpr, '#2C7BB6')
    mx1=max_value(data1)
    print("MX1",mx1)
    mx2=max_value(data2)
    maximum=max(mx1,mx2)
    print(maximum)
    # draw temporary red and blue lines and use them to create a legend
    plt.plot([], c='#D7191C', label='isONform')
    plt.plot([], c='#2C7BB6', label='RATTLE')
    plt.legend()
    plt.title(header)
    plt.xticks(range(0, len(ticks) * 2, 2), ticks)
    plt.xlim(-2, len(ticks) * 2)
    plt.ylim(0, maximum)
    plt.tight_layout()
    file = header + ".png"
    filename = os.path.join(outfolder, file)
    plt.savefig(filename)
    plt.show()


def main(args):
    isExist = os.path.exists(args.outfolder)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(args.outfolder)
    reconstruction_path = os.path.join(args.outfolder, "Reconstruction.csv")
    reconstruction_file = open(reconstruction_path, "w")
    # print("RATTLEPATHS",rattle_path)
    perc_identity_path = os.path.join(args.outfolder, "Perc_identity.csv")
    perc_identity_file = open(perc_identity_path, "w")

    perc_identity_file.write("id,nr_isos,tool,perc_identity\n")
    reconstruction_file.write("id,nr_isos,tool,full,partial\n")
    isoformNumbers = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    # max_iso_nr=args.max_isoforms-1
    max_iso_nr = len(isoformNumbers)
    # we save the precision and recall values for isONform an Rattle in lists
    ison_avg_errors_list = [None] * max_iso_nr
    ison_avg_error_rate = [None] * max_iso_nr
    rattle_avg_errors_list = [None] * max_iso_nr
    rattle_avg_error_rate = [None] * max_iso_nr
    print("IN", args.indir)
    print("OUT", args.outfolder)
    nr_isoforms_dict = {}
    for i,number in enumerate(isoformNumbers):
        nr_isoforms_dict[number]=i
    for analysisfile in os.listdir(args.indir):
        if analysisfile.startswith("results_ison_"):
            print("File", analysisfile)
            file_dir = os.path.join(args.indir, analysisfile)
            avg_perc_identity,perfect_reconstrs,partial = read_analysis(file_dir)
            # tp=len(isonform_dict)
            filename = analysisfile.split(".")[0]
            nr_isos = int(filename.split("_")[2])
            run_id = int(filename.split("_")[3])
            id = str(nr_isos) + "_" + str(run_id)
            perc_identity_file.write("{0},{1},{2},{3}\n".format(id, nr_isos, "isONform", avg_perc_identity))
            # rattle_file.write("precision\n")
            reconstruction_file.write("{0},{1},{2},{3},{4}\n".format(id, nr_isos, "isONform", perfect_reconstrs,partial))
            index = nr_isoforms_dict[nr_isos]
            #print(avg_error_rate)
            #if ison_avg_errors_list[index]:
            #    ison_avg_error_rate[index].append(avg_error_rate)
            #    ison_avg_errors_list[index].append(avg_nr_errors)
            #else:
            #    ison_avg_error_rate[index] = []
            #    ison_avg_error_rate[index].append(avg_error_rate)
            #    ison_avg_errors_list[index] = []
            #    ison_avg_errors_list[index].append(avg_nr_errors)"""
            ison_name = "ison_res_" + str(nr_isos) + "_" + str(run_id) + ".csv"

        elif analysisfile.startswith("results_rattle_"):
            print("File", analysisfile)

            file_dir = os.path.join(args.indir, analysisfile)
            avg_perc_identity,perfect_reconstrs,partial = read_analysis(file_dir)
            # tp = len(rattle_dict)
            filename = analysisfile.split(".")[0]
            nr_isos = int(filename.split("_")[2])
            run_id = int(filename.split("_")[3])
            id = str(nr_isos) + "_" + str(run_id)
            # print(nr_isos)
            # print(run_id)
            perc_identity_file.write("{0},{1},{2},{3}\n".format(id, nr_isos, "RATTLE", avg_perc_identity))
            reconstruction_file.write("{0},{1},{2},{3},{4}\n".format(id, nr_isos, "RATTLE", perfect_reconstrs,partial))
            """
            index = nr_isoforms_dict[nr_isos]
            if rattle_avg_errors_list[index]:
                rattle_avg_errors_list[index].append(avg_nr_errors)
                rattle_avg_error_rate[index].append(avg_error_rate)
            else:
                rattle_avg_errors_list[index] = []
                rattle_avg_errors_list[index].append(avg_nr_errors)
                rattle_avg_error_rate[index] = []
                rattle_avg_error_rate[index].append(avg_error_rate)"""
            rattle_name = "rattle_res_" + str(nr_isos) + "_" + str(run_id) + ".csv"


    # print("ison_prec_dict", ison_precision_list, ", ison_rec_dict", ison_recall_list)
    # print("rattle_prec_dict",rattle_precision_list,", rattle_rec_dict", rattle_recall_list)
    """plt.rcParams["figure.figsize"] = [7.50, 3.50]
    plt.rcParams["figure.autolayout"] = True

    fig, ax = plt.subplots()

    ax.boxplot(ison_precision_dict.values())
    ax.set_xticklabels(ison_precision_dict.keys())
"""
    print(ison_avg_errors_list)
    print(rattle_avg_errors_list)
    print(ison_avg_error_rate)
    print(rattle_avg_error_rate)
    reconstruction_file.close()
    """tot_errors_file.close()"""

    ticks = ['5', '10', '15', '20', '25', '30', '35', '40', '45', '50']
    #plot_data("AverageErrorRate", ison_avg_error_rate, rattle_avg_error_rate, ticks, args.outfolder)
    #plot_data("TotalNrErrors", ison_avg_errors_list, rattle_avg_errors_list, ticks, args.outfolder)
    #plot_data("AverageErrorRate", ison_avg_error_rate, rattle_avg_error_rate, ticks, args.outfolder)
    #plot_with_seaborn_total_errors(ticks, args.outfolder)
    #plot_with_seaborn_error_rate(ticks, args.outfolder)
    perc_identity_file.close()
    plot_with_seaborn_perc_id(ticks,args.outfolder)
    plot_with_seaborn_reconstruction(ticks,args.outfolder)
    # print(original_key_list)
    # print(csv_obj)
    # SIRV_dict=record_lines(csv_obj)
    # print(SIRV_dict)
    # write_Infos(SIRV_dict,args.outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        "Generates Graphs from the alignments.")
    # parser.add_argument('original', type=str, help='CSV file. ')
    # parser.add_argument('rattle_name', type=str, help='csv file. ')
    # parser.add_argument('rattle_output', type=str, help='.fastq file. ')
    parser.add_argument('--indir', type=str, help='csv file. ')
    parser.add_argument('--outfolder', type=str, help='csv file. ')
    parser.add_argument('--max_isoforms', type=int, help="The maximum number of isoforms we generate in our files")
    # parser.add_argument('dummy_outfile', type=str, help='csv file. ')
    # parser.add_argument('max_size', type=int, help='Min size of reads.')

    args = parser.parse_args()

    main(args)
