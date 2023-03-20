from __future__ import print_function
import os,sys
import argparse
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
def plot_nics(nic_isON,nic_stringtie,nic_RATTLE):
    venn3([nic_isON,nic_stringtie,nic_RATTLE],('IsONform', 'StringTie2', 'RATTLE'))
    plt.rcParams.update({'font.size': 18})
    plt.title("Predicted NICs")
    plt.savefig("NIC.pdf")
    plt.show()
def plot_fsms(isONform_fsms, stringtie_is_fsms, rattle_fsms):
    venn3([isONform_fsms, stringtie_is_fsms, rattle_fsms], ('IsONform', 'StringTie2', 'RATTLE'))

    plt.rcParams.update({'font.size': 18})
    plt.title("Predicted FSMs")
    plt.savefig("FSM.pdf")
    plt.show()
def main(args):
    isONfile = open(args.ison_file, 'r')
    isONLines = isONfile.readlines()
    rattlefile = open(args.rattle_file, 'r')
    rattleLines = rattlefile.readlines()[1:]
    headerline=isONLines.pop(0)
    isONfilename=args.ison_file.split("/")
    outfolder=args.outfolder
    print(isONfilename)
    print(headerline)
    header_line_list=headerline.split(",")
    print(header_line_list[14])
    isONform_fsms=set()
    stringtie_is_fsms=set()
    rattle_fsms=set()
    ison_fsm_ctr=0
    stringtie_fsm_cter=0
    rattle_fsm_cter=0
    stringtie_ra_fsms = set()
    isONform_NA_ct=0
    stringtie_NA_ct=0
    rattle_NA_ct=0
    isONform_nic_count=0
    stringtie_nic_count=0
    rattle_nic_count=0
    isONform_ism_count = 0
    stringtie_ism_count = 0
    rattle_ism_count = 0
    isONform_nnc_count = 0
    stringtie_nnc_count = 0
    rattle_nnc_count = 0
    nic_isON=set()
    nic_RATTLE=set()
    nic_stringtie=set()

    for isONline in isONLines:
        this_line=isONline.split(",")
        #print(this_line)
        fsm_indicator=this_line[14]
        nic_indicator=this_line[15]
        ism_indicator=this_line[16]
        nnc_indicator=this_line[17]
        nic_id = this_line[20]
        if fsm_indicator=="1":
            new_fsm_id = this_line[21]
            if this_line[1]=="corrected":
                ison_fsm_ctr += 1
                if new_fsm_id !="NA":

                    if not new_fsm_id in isONform_fsms:
                        isONform_fsms.add(new_fsm_id)

                    #else:
                        #print("Double fsm ,corr",new_fsm_id)
                else:
                    isONform_NA_ct+=1
            elif this_line[1]=="original":
                stringtie_fsm_cter += 1
                if new_fsm_id != "NA":

                    if not new_fsm_id in stringtie_is_fsms:

                        stringtie_is_fsms.add(new_fsm_id)

                    #else:
                        #print("Double fsm, orig ",new_fsm_id)
                else:
                    stringtie_NA_ct+=1
        if nic_indicator=="1":
            if this_line[1]=="corrected":
                if not nic_id in nic_isON:
                    nic_isON.add(nic_id)
                isONform_nic_count+=1
            else:
                if not nic_id in nic_stringtie:
                    nic_stringtie.add(nic_id)
                stringtie_nic_count+=1
        if ism_indicator=="1":
            if this_line[1]=="corrected":
                isONform_ism_count+=1
            else:
                stringtie_ism_count+=1
        if nnc_indicator=="1":
            if this_line[1]=="corrected":
                isONform_nnc_count+=1
            else:
                stringtie_nnc_count+=1
    for rattleline in rattleLines:
        this_line = rattleline.split(",")
        #print(this_line)
        fsm_indicator = this_line[14]
        nic_indicator = this_line[15]
        ism_indicator = this_line[16]
        nnc_indicator = this_line[17]
        nic_id=         this_line[20]
        if fsm_indicator=="1":
            new_fsm_id = this_line[21]
            if this_line[1] == "corrected":
                if new_fsm_id != "NA":
                    rattle_fsm_cter += 1
                    if not new_fsm_id in rattle_fsms:
                        rattle_fsms.add(this_line[21])
                else:
                    rattle_NA_ct+=1
            elif this_line[1] == "original":
                if new_fsm_id != "NA":
                    #stringtie_fsm_cter += 1
                    if not new_fsm_id in stringtie_ra_fsms:
                        stringtie_ra_fsms.add(this_line[21])
        if ism_indicator=="1":
            if this_line[1] == "corrected":
                rattle_ism_count+=1
        if nic_indicator=="1":
            if this_line[1] == "corrected":
                if not nic_id in nic_RATTLE:
                    nic_RATTLE.add(nic_id)
                rattle_nic_count+=1
        if nnc_indicator=="1":
            if this_line[1] == "corrected":
                rattle_nnc_count+=1
    print("Number of fsms detected with isonform: ",len(isONform_fsms))
    print("Number of fsms detected with stringtie (is): ", len(stringtie_is_fsms))
    print("Number of fsms detected with rattle: ", len(rattle_fsms))
    print("Number of fsms detected with stringtie (ra): ", len(stringtie_ra_fsms))
    print()
    #print("stringtie equality",stringtie_ra_fsms==stringtie_is_fsms)
    in_both = isONform_fsms & stringtie_is_fsms
    in_botch =rattle_fsms & stringtie_ra_fsms
    ra_iso=rattle_fsms &isONform_fsms
    print("intersection isONform Stringtie",len(in_both))
    print("intersection Rattle Stringtie", len(in_botch))
    print("intersection Rattle isONform",len(ra_iso))

    is_no_st=isONform_fsms.difference(stringtie_is_fsms)
    unique_isON=is_no_st.difference(rattle_fsms)
    ra_no_st=rattle_fsms.difference(stringtie_is_fsms)
    unique_rattle=ra_no_st.difference(isONform_fsms)
    st_no_ison=stringtie_is_fsms.difference(isONform_fsms)
    unique_stringtie=st_no_ison.difference(rattle_fsms)
    in_all = in_both & in_botch
    print("Intersection all three", len(in_all))
    print()
    print("Unique FSMs isONform",len(unique_isON))
    print("Unique FSMs RATTLE", len(unique_rattle))
    print("Unique FSMs stringtie", len(unique_stringtie))

    """print()
    print("Overall lines FSM isON",ison_fsm_ctr)
    print("Overall lines FSM rattle", rattle_fsm_cter)
    print("Overall lines FSM Stringtie", stringtie_fsm_cter)
    print()
    print("IsONform NIC count", isONform_nic_count)
    print("stringtie NIC count", stringtie_nic_count)
    print("RATTLE NIC count", rattle_nic_count)
    print()

    
    print("IsONform ISM count",isONform_ism_count)
    print("stringtie ISM count",stringtie_ism_count)
    print("RATTLE ISM count", rattle_ism_count)
    print()

    print("IsONform NNC count", isONform_nnc_count)
    print("stringtie NNC count", stringtie_nnc_count)
    print("RATTLE NNC count", rattle_nnc_count)
    plt.rcParams.update({'font.size': 18})
    plot_fsms(isONform_fsms, stringtie_is_fsms, rattle_fsms,args.outfolder,id)
    plot_nics(nic_isON,nic_stringtie,nic_RATTLE)
    #venn3([isONform_fsms, stringtie_is_fsms, rattle_fsms], ('IsONform', 'StringTie2', 'RATTLE'))
    #plt.rcParams.update({'font.size': 18})
    #plt.title("Predicted FSMs")
    #plt.savefig("FSM.pdf")
    #plt.show()"""
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate analysis output file.")
    parser.add_argument('ison_file', type=str, help='Path to the evaluation of the isON Pipeline')
    parser.add_argument('rattle_file', type=str, help='Path to the evaluation of the RATTLE Pipeline')
    parser.add_argument('outfolder', type=str, help='Path to the evaluation of the isON Pipeline')
    args = parser.parse_args()
    main(args)
