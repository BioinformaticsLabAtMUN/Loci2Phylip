# Loci2Phylip
Python script for phylogenetic analysis

Loci2Phylip allows to extract loci present in samples based on several criteria. The script takes as input the filename with the loci, and the filename of the arguments file. It will create a file with the selected loci in [Phylip format](http://scikit-bio.org/docs/0.2.3/generated/skbio.io.phylip.html#r168) and a tab-delimited text file with a table indicating the number of loci in common between samples.

## 1. Requisites

Python 3

## 2. Usage

To run the script, one needs to type in the terminal:

```
python3.3 loci_to_phylip_2.py loci_file arguments_file filename_Phylip_output filename_txt_output
```

The loci file contains the loci sequence. The sample ID is indicated on the first column and the sequence is given on the second column. Columns are separed by a tab. Each locus is separated by //. For example, a loci file with three loci would look like this:

```
A11     CAGGCAAAANGGTGTCTTCGATGCT-GGGGAGCCAGGTGCCACGGTACGTGTGGGAGATTGTACCTGTATCAGGTAGGCTGATTCAGCCTGATTT
B12     CAGGCAAAAARGTGTCTTCGATGCTGGGGGAGCCAGTTG-CACGGTA-GTGTGGGAGATTGTACCTGTATCAGGTA-GCTGATTCA-CCTGATTT
//                -                         -                                                          |0|
A11     CAGGAGGTTTGTAGGACCT--GGAGGTCTCCCTGGACGGCACAGGACGACTTCGAGGAGTCTTTCCTTGTTGCAGCCGAGATCGGAAGAGCGGG
A4      CAGGAGGTTTGTAGGACCCCAAAAGGTCTCCCTGGACGTCGCAGGGAGAC-TCGGGGAGTCTTT-TTTGTTGCAGCCGAGATCGGAAGAGCGGG
F10     CAGGAGGTTTGTAGGACCC--GAAGGTCTCCCTGGACG-CACAGGA-GAC-TCGGGGAGTCTTT-CTTGTTGCAGCCGAGATCGGAAGAGCGG-
//                        -  --               - -    --       -          -                            |1|
A11     CAGATCTATCGAACACAAACATTGCCAGTTGGTAACAGATAAGCAAAGACAGCGCCACCGCTCAAACTGAATCTCAAACTGAGCACGA
G10     CAGATC-ATCGAACACAAACATTGCCAGTTGGTACCAGATAAGCAAAGACAGCGCCACCGCTCAAACTGAATCTCAAACTGAGCACGA
//
```
The argument file list first the IDs of the samples (one sample ID per line) and then it will list the conditions to select the loci (one condition per line). A sample condition is the following:
```
minimum of 5 (G8 and/or H8, E8, E4, H2, D2, E2, C7, E7)
```
This indicates that any loci in common between five samples including G8 and/or H8, and any combination of samples  E8, E4, H2, D2, E2, C7, E7 should be selected. A sample [argument file](arguments.txt).

The Phylip output file will have as header the number of samples and the total number of characters in the alignment. After that there will be a sample ID followed by the corresponding sequence for this sample in the next line.

The tab delimited file will be a square matrix indicating the number of loci in common betweenn any sample pair.

## 3. Citing

If you use this script, please cite:

María Esther Nieto-Blázquez, Lourdes Peña-Castillo and Julissa Roncal. (2020) Historical biogeography of Caribbean Podocarpus does not support the progression rule. Under review.
