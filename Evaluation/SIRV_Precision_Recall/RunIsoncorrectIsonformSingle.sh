#!/bin/bash
#Author: Alexander Petri
# Modifications by Kristoffer Sahlin
# This shell program is used to run longterm tests on the algorithm. 
# For each amount of isoforms, this script performs 5 individual runs of 
# generating test data and running the algorithm on it. To find the remaining 
# bugs the script also outputs the test data files, which the algorithm did 
# not work on correctly
# RUN script as: ./generateTestResults.sh Max isoform number nr runs /path/to/input/reference.fa output_root /path/to/isONform.py
#TODO invoke isONcorrect (Problem: we will need to generate folders for each read to fulfill the isONcorrect input scheme)
if [ $# -lt 2 ]; then
    # TODO: print usage
    echo "./RunIsonform.sh  <Max isoform number> <nr runs> </path/to/input/reference.fa> <output_root> </path/to/isONform.py>"
    exit 1
fi

filedirectory=$1
isonform_dir=$2
actual_file=$3
iso_abundance=$4
echo $iso_abundance
mkdir -p $filedirectory
#mkdir -p $filedirectory/errors
mkdir -p $filedirectory/reads
#mkdir -p $filedirectory/isONform
outputfile=$filedirectory/CorrectForm.tsv
echo $filedirectory
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
file=out/mapping.txt
file=$filedirectory/isONform/mapping.txt
num_cores=8
FILES=$filedirectory"/reads/*.fastq"

	echo $file
	# cp $filedirectory/reads/reads.fq $filedirectory
		#we want to figure out how many reads were actually generated
		#read_amount=$(< $filedirectory/reads/reads.fq wc -l)
		#
		
	part2=$(basename "$actual_file")
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
	read_amount=$(< $actual_file wc -l)
	#As fastq entries have 4 lines divide by 4 to get the actual number of reads
	var=4
	true_read_amount=$((read_amount / var))
		#mv $filedirectory/reads.fq $filedirectory/reads_$number.fq
		#run IsONform
		#if e=True
	echo $isonform_dir/main.py
	isONcorrect  --fastq $actual_file  --outfolder $filedirectory/correction/
	mv $filedirectory/correction/corrected_reads.fastq $filedirectory/correction/corr_$number.fastq
	python $isonform_dir/main.py --fastq $filedirectory/correction/corr_$number.fastq --k 9 --w 20 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --delta_len 10 --outfolder $filedirectory/isONcorrform --iso_abundance $iso_abundance
	#FILE=$filedirectory/isONcorrform/cluster_merged.fastq
	#if test -f "$FILE"; then
	mv $filedirectory/isONcorrform/cluster_merged.fastq  $filedirectory/isONcorrform/cl_spoa_$number.fastq
	#else
	#	mv $filedirectory/isONcorrform/spoa0merged.fasta $filedirectory/isONcorrform/cl_spoa_$number.fasta
	#fi
	#mv $filedirectory/isONcorrform/spoa0merged.fastq $filedirectory/isONcorrform/cl_spoa_$number.fastq
	#python -m pyinstrument main.py --fastq $filedirectory/reads_$number.fq --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --delta_len 5 --outfolder $filedirectory/isonform/
		#if e=False
		#python main.py --fastq $filedirectory/reads/isoforms.fa --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --delta_len 3 --outfolder out
		had_issue=$?
		zero=0
		echo $had_issue
		#python main.py --fastq $filedirectory/reads/isoforms.fa --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --max_bubblesize 2 --delta_len 3 --outfolder testout 
		#we count the lines in mapping.txt, which is an output of IsONform. As we have a fasta file we have to divide this number by 2 to retreive the actual number of isoforms we generated

end=$(date +%s)
#echo "Elapsed Time: $(($end-$start)) seconds"	
