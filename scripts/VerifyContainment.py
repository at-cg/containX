#Usage : python3 <zerocoveragebed> <reads-headers> <containedReadIds>
#checks unmapped bed intervals, evaluates reads sampled from them and ensures that they were contained

import sys
import csv
import subprocess

contained=set()
with open(sys.argv[3], 'r') as f:
    Lines = f.read().splitlines() 
    for line in Lines:
        contained.add(line)

#convert headers into annonated bed file 
with open('temp.bed', "w") as outfile:
    cmd = 'cat ' + sys.argv[2] + ' | awk -F \'[=,-]\' \'{print $9\"\\t\"$5\"\\t\"$6\"\\t\"$0}\' | sort -k 1,1 -k2,2n -k3,3nr | tr -d \'>\' ' 
    subprocess.run(cmd, stdout=outfile, shell=True)

with open(sys.argv[1], 'r') as f:
    Lines = f.read().splitlines() 
    for line in Lines:
        #print unmapped bed interval to a file
        with open('temp2.bed', "w") as outfile:
            cmd = 'echo \"' + line + '\"';
            subprocess.run(cmd, stdout=outfile, shell=True)
        with open('temp3.bed', "w") as outfile:
            cmd = 'bedtools intersect -wa -a temp.bed -b temp2.bed'
            subprocess.run(cmd, stdout=outfile, shell=True)
        countOverlappingReads=0
        countContained=0
        listNonContainedReads=[]
        ContainedReads=[]
        with open("temp3.bed", 'r') as f2:
            reader = csv.reader(f2, delimiter='\t')
            for row in reader:
                countOverlappingReads+=1
                if row[3] in contained:
                    countContained+=1
                    ContainedReads.append(row[3])
                else:
                    listNonContainedReads.append(row[3])
        print(line, "\t", "countreads=",countOverlappingReads,",contained=",countContained,"ContainList=",",".join(ContainedReads),",notContainedList=",",".join(listNonContainedReads), sep="", flush=True) 

