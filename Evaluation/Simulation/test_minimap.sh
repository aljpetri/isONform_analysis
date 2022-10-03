#!/bin/bash
rm /home/alexanderpetri/isONform_analysis/Analysis/SimulationsSnake/alignment15_1.sam
rm /home/alexanderpetri/isONform_analysis/Analysis/SimulationsSnake/counts_15_1.csv
minimap2 -ax map-ont -k 4 -w 1 -p 1 --end-bonus 100 /home/alexanderpetri/isONform_analysis/Analysis/SimulationsSnake/isoforms/isoforms_15_1.fa /home/alexanderpetri/isONform_analysis/Analysis/SimulationsSnake/isONcorrform/cl_spoa_15_1.fastq > /home/alexanderpetri/isONform_analysis/Analysis/SimulationsSnake/alignment15_1.sam
python get_error_rates_original.py /home/alexanderpetri/isONform_analysis/Analysis/SimulationsSnake/isoforms/isoforms_15_1.fa /home/alexanderpetri/isONform_analysis/Analysis/SimulationsSnake/alignment15_1.sam /home/alexanderpetri/isONform_analysis/Analysis/SimulationsSnake/counts_15_1.csv
