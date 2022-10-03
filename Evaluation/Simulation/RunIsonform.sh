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
    echo "./RunIsonform.sh  <Max isoform number> <nr runs> </path/to/input/reference.fa> <output_root> </path/to/isONform.py>"
    exit 1
fi
max_iso_nr=$1
nr_runs=$2
input_ref=$3
filedirectory=$4
isonform_dir=$5
mkdir -p $filedirectory
mkdir -p $filedirectory/errors
mkdir -p $filedirectory/reads
mkdir -p $filedirectory/isONform
outputfile=$filedirectory/IsONform.tsv
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


FILES=$filedirectory"/reads/*.fastq"
echo $FILES
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
	read_amount=$(< $file wc -l)
	#As fastq entries have 4 lines divide by 4 to get the actual number of reads
	var=4
	true_read_amount=$((read_amount / var))
		#mv $filedirectory/reads.fq $filedirectory/reads_$number.fq
		#run IsONform
		#if e=True
	echo $isonform_dir/main.py
	python $isonform_dir/main.py --fastq $file --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --delta_len 5 --outfolder $filedirectory/isONform --max_seqs 500
	#python -m pyinstrument main.py --fastq $filedirectory/reads_$number.fq --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --delta_len 5 --outfolder $filedirectory/isonform/
		#if e=False
		#python main.py --fastq $filedirectory/reads/isoforms.fa --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --delta_len 3 --outfolder out
		had_issue=$?
		zero=0
		echo $had_issue
		#python main.py --fastq $filedirectory/reads/isoforms.fa --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --max_bubblesize 2 --delta_len 3 --outfolder testout 
		#we count the lines in mapping.txt, which is an output of IsONform. As we have a fasta file we have to divide this number by 2 to retreive the actual number of isoforms we generated
		#TODO: we have to change $file to another variable holding the name "transcriptome_mapping.txt"
		res=$(< "$filedirectory/isONform/cluster_mapping.txt" wc -l)
		cp $filedirectory/isONform/cluster_mapping.txt $filedirectory/isONform/cl_map_$number.txt
		result=$(($res/2))
		if [[ "$had_issue" != "$zero" ]]
		then
		errorcounter=$((errorcounter+1))
		result=-1
		cp $file $filedirectory/errors/error_${errorcounter}.fq
		#cp $filedirectory/reads/reads.fq $filedirectory/errors/error_${errorcounter}_$number.fq

		fi
		#echo "$result"
		echo -e "$i \t $j \t $result \t $true_read_amount \t" >> $outputfile
		num2=2
		otherres=$((result * num2))
done
