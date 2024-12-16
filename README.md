# README for sam_file_analysis project

## overview

Project assignment for the SYSTEME class, Master in Bioinformatics of Montpellier, France. 

This bioinformatic project's objetcive focus is the processing and analysis of a SAM file. 
Input = sam file (file.sam)
Output = text file with the summary of the sam file analysis (results.txt)

## execution

./sam_project.sh path_to_sam_file/file.sam

## sam_project.sh

This bash script contains two parts : 
- part 1 : check if the file passed in argument is actually a sam file or not
- part 2 : if the file is a sam file then the script runs the python script (in the same repository) to analyse the reads

## SAM_file_analysis.py

This python script contains the analysis functions for the sam file : 

1. read_parse(file)
Description: Parses the SAM file and extracts key information such as query names, flags, reference sequences, positions, mapping quality, CIGAR strings, and sequences.
Parameters:
file: Path to the SAM file.
Returns: A dictionary containing the extracted data.
2. flagBinary(flag)
Description: Converts a SAM FLAG value to its binary representation and returns the bits in reverse order.
Parameters:
flag: The FLAG value as an integer.
Returns: A list of characters representing the binary bits.
3. add_binary_to_dict(data)
Description: Adds the binary representation of FLAG values to the parsed data dictionary.
Parameters:
data: The dictionary containing parsed SAM data.
4. count_reads_having(reads, FLAG)
Description: Counts the number of reads that have a specific FLAG value using a bitwise AND operation.
Parameters:
reads: The dictionary containing parsed SAM data.
FLAG: The FLAG value to check against.
Returns: The count of reads matching the FLAG.
5. is_partially_mapped(cigar)
Description: Checks if the CIGAR string indicates that a read is partially mapped (contains 'S' or 'H' for soft or hard clips).
Parameters:
cigar: The CIGAR string.
Returns: True if partially mapped, False otherwise.
6. count_partially_mapped(reads)
Description: Counts the number of partially mapped reads based on their CIGAR values.
Parameters:
reads: The dictionary containing parsed SAM data.
Returns: The count of partially mapped reads.
7. calculate_percentage(count)
Description: Calculates the percentage of a given count relative to the total number of reads.
Parameters:
count: The count to calculate the percentage for.
Returns: The calculated percentage.
8. summarise_results(file)
Description: Summarizes the results of the analysis, including counts and percentages of total, unmapped, partially mapped, and mapped reads, as well as quality metrics.
Parameters:
file: Path to the SAM file.
Returns: None (prints a summary DataFrame).

