# Igclone

# Automated generation of immune receptor chains for cloning 

This program was created to genenerate full length TCR or BCR chains for cloning purposes using information typically found in 10X or Adaptive datasets.
Receptor regions are compiled using FASTA sequences from
[IMGT](http://www.imgt.org).

# Getting Started

### Dependencies

This program was implemented for use with python version 3.8 

biopython -- `pip install biopython`
difflib -- `pip install difflib`
pandas -- `pip install pandas`

The provided IMGT FASTA files for region and species were download from [here](http://www.imgt.org/vquest/refseqh.html).  

### Installation

The following commands install and run the development version of Igclone:

```sh
git clone https://github.com/edavis71/Igclone
cd Igclone/
```

### Input

Input file should contain a sample identifier, CDR3 amino acid sequence, V and J gene names in the following format:

| id	| v |	cdr3	| j |
| --- | --- | --- | --- |
| tcr1_a	| TRAV8-3*01	| CAVGRGGSQGNLIF	| TRAJ42*01 |
| tcr1_b	| TRBV4-2*01	| CASSQYQSLVRGNNEQFF	| TRBJ2-1*01 |

An example input file can be found under chain_example.csv

### Output

For each input file a csv is returned in the following format:

| id	| v |	cdr3	| j | chain |
| --- | --- | --- | --- | --- |
| tcr1_a	| TRAV8-3*01	| CAVGRGGSQGNLIF	| TRAJ42*01 | MLLELIPLLGI....NLIFGKGTKLSVKP |
| tcr1_b	| TRBV4-2*01	| CASSQYQSLVRGNNEQFF	| TRBJ2-1*01 | MGCRLLCCAVLC....EQFFGPGTRLTVL |


## Feedback

Please send feedback to Emily Davis-Marcisak:
<edavis71@jhu.edu>

## Licensing

XXX is liscenced under the [MIT](link) license.

