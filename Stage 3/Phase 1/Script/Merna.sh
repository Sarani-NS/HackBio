#!/bin/bash

# Create directories if they don't exist
mkdir -p Merna
mkdir -p biocomputing

# Download the sequences (replace wget with curl)
curl -o sequence.fasta "https://www.ncbi.nlm.nih.gov/nuccore/X91251.1?report=fasta&log$=seqview&format=text"

# Count nucleotides (this assumes sequence.fasta is in the current directory)
num_A=$(sed 1d sequence.fasta | grep -o -i 'A' | wc -l)
num_G=$(sed 1d sequence.fasta | grep -o -i 'G' | wc -l)
num_T=$(sed 1d sequence.fasta | grep -o -i 'T' | wc -l)
num_C=$(sed 1d sequence.fasta | grep -o -i 'C' | wc -l)

# Output results to the console
echo "Number of A: $num_A"
echo "Number of G: $num_G"
echo "Number of T: $num_T"
echo "Number of C: $num_C"

# Optionally, save the results to a file
echo "Number of A: $num_A" >> results.txt
echo "Number of G: $num_G" >> results.txt
echo "Number of T: $num_T" >> results.txt
echo "Number of C: $num_C" >> results.txt
