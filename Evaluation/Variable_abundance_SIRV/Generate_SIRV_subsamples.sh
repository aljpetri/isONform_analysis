#!/bin/bash
#Author: Alexander Petri
# Modifications by Kristoffer Sahlin
# This shell program is used to run longterm tests on the algorithm. 
# For each amount of isoforms, this script performs 5 individual runs of 
# generating test data and running the algorithm on it. To find the remaining 
# bugs the script also outputs the test data files, which the algorithm did 
# not work on correctly
# RUN script as: ./generateTestResults.sh /path/to/input/reference.fa /path/to/output/folder
if [ $# -lt 2 ]; then
    # TODO: print usage
    echo "./generateTestResults.sh  </path/to/input/reference.fa> <output_root>"
    exit 1
fi
max_iso_nr=$1
nr_runs=$2
input_ref=$3
filedirectory=$4
alignment=$5
reference=$6
mkdir -p $filedirectory
#mkdir -p $filedirectory/errors
mkdir -p $filedirectory/reads
mkdir -p $filedirectory/references
outputfile=$filedirectory/resultserror1.tsv
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
file=$filedirectory/isonform/mapping.txt
isoformNumbers=(5 10 15 20 25 30 35 40 45 50)
#iterate over different numbers of isoforms
for i in ${isoformNumbers[@]}; do
	echo $i
	#echo "Generating $i TestIsoforms" >>results.tsv
	#for each amount of isoforms we would like to run 15 tests to make sure our algo works stable
	for((j=1;j<=$nr_runs;j++))
	do
		echo $i_$j
		outputs=0
		#python generateTestCases.py --ref $input_ref --sim_genome_len 1344 --nr_reads 20 --outfolder testout --coords 50 100 150 200 250 300 350 400 450 500 --probs 0.4 0.4 0.4 0.4 0.4 --n_isoforms 2 --e True
		#run the test case generation script with the parameters needed

		number="${i}_${j}"
		############ COMMENT THE FOLLOWING TWO LINES FOR BUGFIXING ON IDENTICAL READ FILES ############
		###############################################################################################
		python sirv_subsample_isoforms.py --fastq $input_ref  --outfile $filedirectory/reads/reads_$number.fastq  --depth_dist exp --nr_isoforms $i --alignments $alignment --sirv_ref $reference
		#cp $filedirectory/isoforms/reads.fq $filedirectory/reads/reads.fq
		mv $filedirectory/reads/reference.fasta $filedirectory/references/reference_$number.fasta
	done
done
touch dummyfile

