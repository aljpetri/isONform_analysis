import argparse
import os

from recordclass import recordclass
def main(args):
    in_directory=""
    Spoa_attributes = recordclass('spoa_attributes', "id_count_sum id_counts consensus_names in_output")
    out_file_name=args.outfile
    #in_directory = args.in_dir
    spoa_count_name=args.spoacount
    batch_count_name=args.batchcount
    #out_file_object = open(os.path.join(in_directory,out_file_name, 'a')
    spoa_dict={}
    batch_dict={}
    spoa_file=open(os.path.join(in_directory, spoa_count_name),'r')
    #next(spoa_file)
    spoa_lines = spoa_file.readlines()[1:]

    for spoa_line in spoa_lines:
        spoa_entries=spoa_line.split(',')
        spoa_consensus_name=spoa_entries[0]
        spoa_id=int(spoa_entries[1])
        spoa_id_count=int(spoa_entries[7])
        if not spoa_id in spoa_dict:
            spoa_attributes=Spoa_attributes(spoa_id_count,[spoa_id_count],[spoa_consensus_name],False)
            spoa_dict[spoa_id]=spoa_attributes
        else:
            spoa_dict[spoa_id].id_count_sum+=spoa_id_count
            spoa_dict[spoa_id].id_counts.append(spoa_id_count)
            spoa_dict[spoa_id].consensus_names.append(spoa_consensus_name)
    batch_file = open(os.path.join(in_directory, batch_count_name), 'r')
    batch_lines = batch_file.readlines()
    for batch_line in batch_lines:
        batch_entries = batch_line.split(',')
        batch_id = batch_entries[0]
        batch_id_count = batch_entries[1]
        batch_dict[batch_id]=batch_id_count
    file_object = open(os.path.join(in_directory,out_file_name), 'a')
    file_object.write(args.spoacount,'\n')
    for batch_id, batch_id_count in batch_dict.items():
        if batch_id in spoa_dict:
            file_object.write(batch_id,": batch ",batch_id_count," spoa ",spoa_dict[batch_id].id_count_sum," : ", spoa_dict[batch_id].id_counts )
            spoa_dict[batch_id].in_output=True
        else:
            file_object.write(batch_id,": batch ",batch_id_count)
    for spoa_id, spoa_attributes in spoa_dict.items():
        if not spoa_attributes.in_output:
            file_object.write("spoa ",spoa_id," ",spoa_dict[batch_id].id_count_sum," : ", spoa_dict[batch_id].id_counts )
    file_object.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="IsONform analysis - Isoform abundance",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--spoacount', type=str, default=False,help='Name of the spoa count file')
    parser.add_argument('--batchcount', type=str, default=False, help='Name of the batch count file')
    parser.add_argument('--outfile', type=str, default=False, help='Name of the batch count file')
    #parser.add_argument('--output', dest="nr_cores", type=int, default=8, help='Name of the output folder')
    #parser.add_argument('--in_dir', type=str, default=False, help='folder in which the files are found also used as output path')
    args = parser.parse_args()
    main(args)