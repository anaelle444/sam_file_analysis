# README for sam_file_analysis project

## overview

Project assignment for the SYSTEME class, Master in Bioinformatics of Montpellier, France. 

This bioinformatic project's objetcive focus is the processing and analysis of a SAM file. 
- Input = sam file (file.sam) 
- Output = text file with the summary of the sam file analysis (results.txt)

## execution

comman in terminal, executed under macOS Sequoia Version 15.1.1 (24B91) with Python 3.12.7 : 
- ./sam_project.sh path_to_sam_file/file.sam

## sam_project.sh

This bash script contains two parts : 
- part 1 : check if the file passed in argument is actually a sam file or not
- part 2 : if the file is a sam file then the script runs the python script (in the same repository) to analyse the reads

## SAM_file_analysis.py

This python script contains the analysis functions for the sam file : 

1. read_parse(file)  
   - Description: parses the SAM file and extracts key information such as query names, flags, reference sequences, positions, mapping quality, CIGAR strings, and sequences.  
   - Parameters:  
     - file: path to the SAM file.  
   - Returns:  
     - a dictionary containing the extracted data.

2. flagBinary(flag)  
   - Description: converts a SAM FLAG value to its binary representation and returns the bits in reverse order.  
   - Parameters:  
     - flag: the FLAG value as an integer.  
   - Returns:  
     - a list of characters representing the binary bits.

4. count_reads_having(reads, FLAG)  
   - Description: counts the number of reads that have a specific FLAG value using a bitwise AND operation.  
   - Parameters:  
     - reads: the dictionary containing parsed SAM data.  
     - FLAG: the FLAG value to check against.  
   - Returns:  
     - the count of reads matching the FLAG.

5. is_partially_mapped(cigar)  
   - Description: checks if the CIGAR string indicates that a read is partially mapped (contains 'S' or 'H' for soft or hard clips).  
   - Parameters:  
     - cigar: the CIGAR string.  
   - Returns:  
     - true if partially mapped, False otherwise.

6. count_partially_mapped(reads)  
   - Description: counts the number of partially mapped reads based on their CIGAR values.  
   - Parameters:  
     - reads: the dictionary containing parsed SAM data.  
   - Returns:  
     - the count of partially mapped reads.

7. calculate_percentage(count)  
   - Description: calculates the percentage of a given count relative to the total number of reads.  
   - Parameters:  
     - count: the count to calculate the percentage for.  
   - Returns:  
     - the calculated percentage.

8. summarise_results(file)  
   - Description: summarizes the results of the analysis, including counts and percentages of total, unmapped, partially mapped, and mapped reads, as well as quality metrics.  
   - Parameters:  
     - file: path to the SAM file.  
   - Returns:  
     - none (prints a summary DataFrame).
