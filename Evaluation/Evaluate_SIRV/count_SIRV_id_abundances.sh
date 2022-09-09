#!/bin/bash

#Author: Alexander Petri
# Modifications by Kristoffer Sahlin

# This shell program is used to verify algorithm outputs. 


# RUN script as: ./count_SIRV_id_abundances.sh <Inputfile> <Outputfile>

if [ $# -lt 1 ]; then
    # TODO: print usage
    echo " <inputfile> <outputfile>"
    exit 1
fi
filename=$1
outputfile=$2
array=("SIRV618" "SIRV617" "SIRV608" "SIRV611" "SIRV614" "SIRV609" "SIRV616" "SIRV606" "SIRV602" "SIRV607" "SIRV615" "SIRV605" "SIRV610" "SIRV613" "SIRV601" "SIRV612" "SIRV604" "SIRV603")
echo "Starting the counting"
for i in "${array[@]}"; do   # The quotes are necessary here
    output="$( grep -wc "$i" $filename)"
    echo -e "$i \t $output \t">>$outputfile
done
