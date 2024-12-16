#!/bin/bash

file="$1"

# part 1: check if the file is an actual SAM file before running the rest

controle=false

if [ -f "$file" ]; then
    if [[ "$file" == *.sam ]]; then
        echo "'$file' is a SAM file"
        controle=true
    else
        echo "'$file' is a regular file but not a SAM file"
    fi
else
    echo "'$file' is not a regular file or does not exist"
fi

# part 2: if the file passes the control, run the python analysis

if [ "$controle" = true ]; then
    python3 ./SAM_file_analysis.py "$file" > results.txt
else
    echo "'$file' is not a SAM file or does not exist; we cannot run the analysis"
fi
