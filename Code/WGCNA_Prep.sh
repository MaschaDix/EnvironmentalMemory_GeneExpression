# Preparing kallisto output (command line) for WGCNA in R

# 1. get a matrix containing tpm values of all treatments and all genes

# in /ibex/scratch/projects/c2116/mascha_thesis/EnvironmentalMemory/GeneExpression/quant

# get genes
cut -f1 A1/abundance.tsv > genes.txt

# get the tpm counts for each sample
for folder in A1 A1b A2 A2b B1 B2 B3 C1 C2 C3 D1 D2 D3 E1 E2 E3 F1 F2 F3 H1 H2 H3 I1 I1b I2 I2b; do
  cut -f 5 ${folder}/abundance.tsv > ${folder}_column5.txt;
done

# combine them all into one file
paste -d ',' genes.txt A1_column5.txt A1b_column5.txt A2_column5.txt A2b_column5.txt B1_column5.txt B2_column5.txt B3_column5.txt C1_column5.txt C2_column5.txt C3_column5.txt D1_column5.txt D2_column5.txt D3_column5.txt E1_column5.txt E2_column5.txt E3_column5.txt F1_column5.txt F2_column5.txt F3_column5.txt H1_column5.txt H2_column5.txt H3_column5.txt I1_column5.txt I1b_column5.txt I2_column5.txt I2b_column5.txt > combined_old_header.csv

# delete the header
tail -n +2 combined_old_header.csv > combined_no_header.csv

# make a new header with the treatment names
echo "Gene_ID, STP_C1, STP_C2, STP_C3, STP_C4, LTP_H1, LTP_H2, LTP_H3, N_H1, N_H2, N_H3, A_H1, A_H2, A_H3, LTP_C1, LTP_C2, LTP_C3, A_C1, A_C2, A_C3, N_C1, N_C2, N_C3, STP_H1, STP_H2, STP_H3, STP_H4" > headers.csv
cat headers.csv combined_no_header.csv > GE_WGCNA_input.csv

# clean up:
rm *_column5.txt
rm combined_*_header.csv
rm genes.txt
rm headers.csv
