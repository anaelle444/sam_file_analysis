#script 1 : is the file an actual sam file before running the rest 

#!/bin/bash

file=$1

if [ -f "$file" ]; then
    if [[ "$file" == *.sam ]]; then
        echo "'$file' is a sam file";
    else
        echo "'$file' is a regular file but not a sam file"
    fi
else
    echo "'$file' is not a regular file or does not exist"
fi

#add a boolean variable if sam file -> TRUE and rest of the workflow ok

#script 2 :

#!/usr/bin/python3
#-*- coding : utf-8 -*-

############### IMPORT MODULES ###############

import os
import sys
import re 

path="/Users/anaelle/Desktop/HAI724I/mapping.sam"


#FUNCTIONS TO :

    
##PARSE THE SAM FILE 
#to have columns
#this function return a dictionary containing parsed fields

def parse_sam(file_path):
    # Fields
    tags_main = [
        "QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"
    ]
    data = {field: [] for field in tags_main}

    with open(file_path, "r") as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")

            data["QNAME"].append(fields[0])
            data["FLAG"].append(int(fields[1]))
            data["RNAME"].append(fields[2])
            data["POS"].append(int(fields[3]))
            data["MAPQ"].append(int(fields[4]))
            data["CIGAR"].append(fields[5])
            data["RNEXT"].append(fields[6])
            data["PNEXT"].append(int(fields[7]))
            data["TLEN"].append(int(fields[8]))
            data["SEQ"].append(fields[9])
            data["QUAL"].append(fields[10])

            optional_fields = fields[11:]
            for field in optional_fields:
                try:
                    tag, type_, value = field.split(":", 2)  # Diviser en 3 parties : TAG, TYPE, et VALUE
                except ValueError:
                    # Si le champ n'a pas trois parties, on l'ignore ou on le signale
                    print(f"Tag malformé ou manquant : {field}")
                    continue

                # Handling xtra tags if not found
                if tag not in data:
                    data[tag] = []
                data[tag].append(value)
    return data

reads=parse_sam(path)


##CNONVERT THE FLAG INTO BINARY
def flagBinary(flag) :
    flagB = bin(int(flag)) # Transform the integer into a binary.
    flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagB = list(flagB) 
    if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
        add = 12 - len(flagB) # Difference between the maximal flag size (12) and the length of the binary flag.
        for t in range(add):
            flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.
    return flagB

###add the binary flag to the dico : 
reads["BINARY_FLAG"]=[]
for flag in reads["FLAG"]:
    binary_flag = flagBinary(flag)  # Compute binary flag
    reads["BINARY_FLAG"].append(binary_flag)

#total number of reads :
def total_reads(reads):
    return len(reads["FLAG"])

#how many reads for each flag : 
def list_reads_flags(reads):
    flag_count={}
    for flag in reads["FLAG"]:
        flag_count[flag] = flag_count.get(flag, 0) + 1
    return flag_count

#how many reads for this flag : THIS WORKS!!!    
def reads_having(reads, FLAG):
    reads_having = 0
    for flag in reads["FLAG"]:
        if flag & FLAG :
            reads_having += 1
    return reads_having


#partially mapped
def is_partially_mapped(cigar):         # Vérifie si le CIGAR contient des softs clips S ou des hard clips H
    return "S" in cigar or "H" in cigar

def count_partially_mapped(reads):
    count_partially=0
    for cigar in reads["CIGAR"]:
        if is_partially_mapped(cigar):
            count_partially += 1
    return count_partially

 


