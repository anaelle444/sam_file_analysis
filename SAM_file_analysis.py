#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 22:08:06 2024

@author: anaelle
"""

import pandas as pd
import operator
import matplotlib.pyplot as plt
import sys


#file = "/Users/anaelle/Desktop/HAI724I/mapping.sam" #path to the sam file
file = sys.argv[1]

##parsing the sam file

def read_parse(file):
    '''
    This function extracts infotmation from the sam file for analysing. It extracts the qname, flag, ref, pos, mapq, cigar and seq 
    but one can easily extract other columns by adding them into the tags_main list and the loop under the following format : data["col"].append(field[x])
    with col being the name of the column added in the tags_main list and x being the index of the coloumn in the original sam file
    '''
    tags_main = [
        "QNAME", "FLAG", "REF", "POS", "MAPQ", "CIGAR", "SEQ"
    ]
    data = {field: [] for field in tags_main}

    with open(file, "r") as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            data["QNAME"].append(fields[0])
            data["FLAG"].append(int(fields[1]))
            data["REF"].append(fields[2])
            data["POS"].append(int(fields[3]))
            data["MAPQ"].append(int(fields[4]))
            data["CIGAR"].append(fields[5])
            data["SEQ"].append(fields[9])
    return data

def flagBinary(flag):
    '''
    This function converts a SAM FLAG value to its binary representation and returns the bits in reverse order : the least significant bit is on the right.
    It then reverses the order of the bits with [::-1] and returns the result as a list of characters (with the bit 1 on the left for easier understanding), where each character represents a bit ('0' or '1'). 
    '''
    flagB = bin(int(flag))[2:]
    flagB = flagB.zfill(12)
    return list(flagB[::-1])

def add_binary_to_dict(data):
    '''
    This function adds the binary representation given by flagBinary() to the dictionary we got with read_parse(). It adds a new column in the dictionary and then iterate throught the "FLAG" column to convert the FLAG value into its binary representaion, and adds it to the new column.
    '''
    data["BINARY_FLAG"]=[]
    for flag in data["FLAG"]:
        binary_flag = flagBinary(flag)
        data["BINARY_FLAG"].append(binary_flag)


##functions to count : 

def count_any1(data, col, threshold): # operator : &
    '''
    this function counts the number of reads in the dictionary "data" where the value from the "col" colunm is contained (operator &) in the "threshold" 

    '''
    count=0
    for i, value in enumerate(data[col]):  # loop through the column indicated in parameter
        if value & threshold:                     
            count += 1
    # print(count, "reads where", col, "&", threshold )
    return count

def count_any2(data, col, threshold, op): # any operator
    '''
    this function counts the number of reads in the dictionary "data" wwhere the value from the "col" column satisfies the condition defined by the operator "op" and the "threshold".
    '''
    count = 0
    for value in data[col]:  # loop through the column values
        if op(value, threshold):  # apply the operator dynamically
            count += 1
    return count


##functions to filter : 

def filter_sup_anything(data, col, threshold): # filter out < threshold
    '''
        this function creates a new dictionary (filtered_data) containing only the reads satisfying the filtering conditions
        filter out reads < threshold
    '''    
    tags_main = ["QNAME", "FLAG", "REF", "POS", "MAPQ", "CIGAR", "SEQ"]
    filtered_data = {field: [] for field in tags_main if field in data}
    
    for i, value in enumerate(data[col]):  # loop through the column indicated in parameter
        if value >= threshold:             # apply threshold
            for key in filtered_data.keys():      # add corresponding data
                filtered_data[key].append(data[key][i])
                
    return filtered_data

def filter_inf_anything(data, col, threshold): # filter out > threshold

    '''
        this function creates a new dictionary (filtered_data) containing only the reads satisfying the filtering conditions
        filter out reads > threshold
    ''' 
    tags_main = ["QNAME", "FLAG", "REF", "POS", "MAPQ", "CIGAR", "SEQ"]
    filtered_data = {field: [] for field in tags_main if field in data}

    for i, value in enumerate(data[col]):
        if value <= threshold:  
            for key in filtered_data.keys():  
                filtered_data[key].append(data[key][i])

    return filtered_data

def filter_operator(data, col, threshold, op):  # filter dynamic operator CHECK WITH FORMAT OF VALUES!!!
    '''
        this function creates a new dictionary (filtered_data) containing only the reads where the value in the specified column satisfies the condition defined by the operator "op" and the "threshold"
        operators : operator. /le for <= /eq for == /ge for >= /and_ for &
    '''   
    tags_main = ["QNAME", "FLAG", "REF", "POS", "MAPQ", "CIGAR", "SEQ"]
    filtered_data = {field: [] for field in tags_main if field in data}
    
    for i, value in enumerate(data[col]):  # loop through the column indicated in parameter
        if op(value, threshold):             # filtering condition on value depending on the threshold and the operator specified by user
            for key in filtered_data.keys():      # add corresponding data
                filtered_data[key].append(data[key][i])
                
    return filtered_data
    

##functions to count reads 

def total_reads(reads):
    ''' this function counts the total number of reads in the dictionary as the number of values in "FLAG" is also the number of reads (each read has a FLAG value to it)
    '''
    return len(reads["FLAG"])

def count_bitAND(data, col, threshold):
    ''' this function counts the number of entries in the specified column that satisfy the condition where the value in the column
        bitwise ANDed with the threshold is non-zero. The function returns the total count of such entries.       
    '''
    count=0
    for i, value in enumerate(data[col]):                                       # loop through the column indicated in parameter
        if value & threshold:                     
            count += 1
    # print(count, "reads where", col, "&", threshold )
    return count

def count_any_op(data, col, threshold, op):
    '''
        this function counts the number of entries in the specified column where the value satisfies the condition defined by the dynamic operator "op" and the "threshold"

    '''
    count = 0
    for value in data[col]:  # Loop through the column values
        if op(value, threshold):  # Apply the operator dynamically
            count += 1
    return count


def count_reads_having(reads, FLAG):
    '''     this function counts the number of reads where their flag value & the FLAG (bitwise AND operation)
    '''
    c_reads_having = 0
    for flag in reads["FLAG"]:
        if flag & FLAG :
            c_reads_having += 1
    return c_reads_having


##functions to si if partially mapped 

def is_partially_mapped(cigar):
    '''    check if the CIGAR value contains 'S' or 'H' for soft and hard clips, meaning the reads is cliped, so partially mapped
    '''
    return "S" in cigar or "H" in cigar

def count_partially_mapped(reads):
    '''     this function counts the number of partially mapped reads in a dictionary based on the CIGAR value, using the is_partially_mapped() function
    '''
    count_partially=0
    for cigar in reads["CIGAR"]:
        if is_partially_mapped(cigar):
            count_partially += 1
    return count_partially

##function for graphs

def plot_read_proportions(var1, var2, var3): #i wanted to extract the name of the variable to pass it as strings in the labels but couldn't figure it out sorry
    name1=str(var1)
    name2=str(var2)
    name3=str(var3)
    labels = [name1, name2, name3]
    sizes = [var1, var2, var3]
    colors = ['#4CAF50', '#FF6347', '#FFD700']  
    explode = (0.1, 0, 0) 

    plt.figure(figsize=(7, 7))
    plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=True, startangle=140)
    plt.axis('equal')  # Equal aspect ratio ensures that pie chart is drawn as a circle.
    plt.title(f'Proportions of Reads: {name1}, {name2}, {name3}')  # Correctly format the title
    plt.show()


def summarise_results(file):
    sam_data = read_parse(file)  # load and parse the sam file
    
    summary = {
        "Description": [],
        "Count": [],
        "Percentage": []
    }

    total_reads_count = total_reads(sam_data)
    def calculate_percentage(count):
        return (count / total_reads_count * 100)

#mapping : total = unmapped + partially mapped + mapped
    total_reads_count = total_reads(sam_data)
    summary["Description"].append("Total Reads")
    summary["Count"].append(total_reads_count)
    summary["Percentage"].append(100)

    unmapped_reads = count_reads_having(sam_data, 4)
    summary["Description"].append("Unmapped Reads")
    summary["Count"].append(unmapped_reads)
    summary["Percentage"].append(calculate_percentage(unmapped_reads))

    partially_mapped = count_partially_mapped(sam_data)
    summary["Description"].append("Partially Mapped Reads")
    summary["Count"].append(partially_mapped)
    summary["Percentage"].append(calculate_percentage(partially_mapped))

    mapped_reads = total_reads_count - unmapped_reads - partially_mapped
    summary["Description"].append("Mapped Reads")
    summary["Count"].append(mapped_reads)
    summary["Percentage"].append(calculate_percentage(mapped_reads))

#read pairs
    properly_paired_reads = count_any2(sam_data, "FLAG", 99, operator.eq) + count_any2(sam_data, "FLAG", 83, operator.eq)
    summary["Description"].append("Properly Paired Reads")
    summary["Count"].append(properly_paired_reads)
    summary["Percentage"].append(calculate_percentage(properly_paired_reads))

#mapping quality score
    under_30_mapq = filter_operator(sam_data, "MAPQ", 30, operator.le)
    high_quality_reads = total_reads_count - len(under_30_mapq["MAPQ"])
    summary["Description"].append("High Quality Reads (MAPQ >= 30)")
    summary["Count"].append(high_quality_reads)
    summary["Percentage"].append(calculate_percentage(high_quality_reads))

    ten_to_30_mapq = filter_operator(under_30_mapq, "MAPQ", 10, operator.ge)
    medium_quality_reads = len(ten_to_30_mapq["MAPQ"])
    summary["Description"].append("Medium Quality Reads (10 <= MAPQ < 30)")
    summary["Count"].append(medium_quality_reads)
    summary["Percentage"].append(calculate_percentage(medium_quality_reads))

    low_quality_reads = count_any2(under_30_mapq, "MAPQ", 10, operator.lt)
    summary["Description"].append("Low Quality Reads (MAPQ < 10)")
    summary["Count"].append(low_quality_reads)
    summary["Percentage"].append(calculate_percentage(low_quality_reads))

#alignment 
    aligned_to_ref = count_any2(sam_data, "REF", "Reference", operator.eq)
    summary["Description"].append("Reads aligned to the Reference")
    summary["Count"].append(aligned_to_ref)
    summary["Percentage"].append(calculate_percentage(aligned_to_ref))

    aligned_to_star = count_any2(sam_data, "REF", "*", operator.eq)
    summary["Description"].append("Reads aligned to *")
    summary["Count"].append(aligned_to_star)
    summary["Percentage"].append(calculate_percentage(aligned_to_star))

#dataframe for viewing
    summary_df = pd.DataFrame(summary)

    print(summary_df)

summarise_results(file)
