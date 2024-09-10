# miRAl
Optimized tool to perform the Needleman-Wunsch algorithm to align mature miRNA sequences.

### Installation using conda:
```conda install -c molinerislab miRAl```

### Description:
This version of the Needleman-Wunsch algorithm has been optimized to align miRNA sequences, in particular it allows to weight the alignment of the seed region of the miRNA differently from the rest of the sequence.

### Usage and options:
```
everygff2tab [OPTIONS] ATTRIBUTE.. <GFF2/GFF3 >TSV

Options:
  -h, --help: show this help message and exit
  -l, --list-attributes: list all the available attributes
  -m STRING, --missing-value=STRING: write this STRING when an attribute is missing (default: NA)
  -g STRING, --gff-version=STRING: STRING indicating the gff version (GFF2, GFF3)(default: GFF2)

ATTRIBUTE: the tag identifying a given value in the "attribute" column in the GFF2/GFF3 file
GFF2/GFF3: GFF file to parse
```


[![DOI](https://zenodo.org/badge/852122116.svg)](https://zenodo.org/doi/10.5281/zenodo.13683069)
