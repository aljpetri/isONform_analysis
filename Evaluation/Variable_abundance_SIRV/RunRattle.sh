#!/bin/bash
#Author: Alexander Petri
# Modifications by Kristoffer Sahlin
# This shell program is used to run longterm tests on the algorithm. 
# For each amount of isoforms, this script performs 5 individual runs of 
# generating test data and running the algorithm on it. To find the remaining 
# bugs the script also outputs the test data files, which the algorithm did 
# not work on correctly
# RUN script as: ./generateTestResults.sh Max isoform number nr runs /path/to/input/reference.fa output_root /path/to/isONform.py
if [ $# -lt 2 ]; then
    # TODO: print usage
    echo "./RunIsonform.sh  <Max isoform number> <nr runs> </path/to/input/reference.fa> <output_root> </path/to/Rattle>"
    exit 1
fi

filedirectory=$1
rattle_loc=$2
mkdir -p $filedirectory/Rattle

outputfile=$filedirectory/rattle/resultserror1.tsv
#if results.tsv already exists 
if [ -s $outputfile ]
then 
#delete all data from the file 
   > $outputfile
else
#add a file results.tsv to write into
   touch $outputfile
fi
#counter used to name error files correctly
errorcounter=0
#write nice header into the output file
#echo -e "Number of Isoforms \t Found isoforms by IsONform \n">>results.tsv
#ls
#define the file we want to use as indicator for our algos performance
rattle_dir=$filedirectory/Rattle

FILES=$filedirectory"/reads/*.fastq"
echo $FILES
start=$(date +%s)
for file in $FILES
do
	echo $file
	# cp $filedirectory/reads/reads.fq $filedirectory
		#we want to figure out how many reads were actually generated
		#read_amount=$(< $filedirectory/reads/reads.fq wc -l)
		#
		
	part2=$(basename "$file")
	echo "Part2"
	echo "$part2"
	arrIN=(${part2//./ })
	echo ${arrIN[0]}  
	fname=(${part2//./ })
	echo ${fname[0]}
	fnamesplit=(${fname[0]//_/ })
	i=${fnamesplit[1]}
	j=${fnamesplit[2]}
	number="${i}_${j}"
	echo "$i \t $j" 
	rattle_out_file=$rattle_dir/transcriptome$i_$j.fq
	read_amount=$(< $file wc -l)
	#As fastq entries have 4 lines divide by 4 to get the actual number of reads
	var=4
	true_read_amount=$((read_amount / var))
		#mv $filedirectory/reads.fq $filedirectory/reads_$number.fq
		#run IsONform
		#if e=True
	echo $isonform_dir/main.py
	echo "RATTLE cluster"
	$rattle_loc/rattle cluster -i $file -t 24 -o $rattle_dir  --iso --min-reads-cluster 0
	echo "RATTLE correct"
	$rattle_loc/rattle correct -i $file -c $rattle_dir"/clusters.out"  -o $rattle_dir -t 24 --min-reads 0
	echo "RATTLE polish"
	$rattle_loc/rattle polish -i $rattle_dir"/consensi.fq" -o $rattle_out_file -t24 
	echo "Polished"
	FILE=$rattle_dir/transcriptome.fq
	if test -f "$FILE"; then
		mv $rattle_dir"/transcriptome.fq" $rattle_dir/rattle_res_$number.fastq
	else
		mv $rattle_dir"/consensi.fq" $rattle_dir/rattle_res_$number.fastq
	fi
	#python $isonform_dir/main.py --fastq $file --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --delta_len 5 --outfolder $filedirectory/isONform
	#python -m pyinstrument main.py --fastq $filedirectory/reads_$number.fq --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --delta_len 5 --outfolder $filedirectory/isonform/
		#if e=False
		#python main.py --fastq $filedirectory/reads/isoforms.fa --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --delta_len 3 --outfolder out
	
done
end=$(date +%s)
echo "Elapsed Time: $(($end-$start)) seconds"	
