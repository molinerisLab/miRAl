# miRAl
Optimized tool to perform the Needleman-Wunsch algorithm to align mature miRNA sequences.

### Installation using conda:
```conda install -c molinerislab miRAl```

### Description:
Needleman-Wunsch alignment given a list of couples of miRNAs in .tsv format structured as follows:
                    mature miRNA1     |    mature miRNA2
and a fasta file containing the sequences of the miRNAs perform the alignmetn between the sequences of the two miRNAs in the couple                                                                

The output file is a .tsv file structured as follows:
                    mature miRNA1     |    mature miRNA2     |    Needleman-Wunsch score
                    
```-f``` ```-a``` options add a second output file in a simil-FASTA format, where the header contains the names of the two aligned sequences while the "body" contains the alignment score and the alignment itself    

### Usage and options:
```
Usage: miRAl [OPTIONS] fasta.fa [fasta2.fa] < couples.tsv >output.tsv

...
```


[![DOI](https://zenodo.org/badge/852122116.svg)](https://zenodo.org/doi/10.5281/zenodo.13683069)
