#!/bin/bash
# acquire data in a parallel manner
# need a 'human3.txt' with all ftp inside

NAME=human3
OUTPUT=./raw
mkdir -p $OUTPUT
>"${NAME}_ftp.txt"
while read p; do
echo "ftp://$p" >> "${NAME}_ftp.txt"
done < "${NAME}.txt"
cd $OUTPUT
cat "../${NAME}_ftp.txt" | parallel wget -c -o ./wget_log_3.txt {}   #-c support resuming downloading after connection failure, -o is to export output into a file, -q means not showing the information
cd ..
