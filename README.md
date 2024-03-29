# genefunction


Returns the summarized function of a gene or list of genes using UniProt.

# Usage
The input can either be a single gene, or an input file of genes whose functions will be output (along with the gene name) to the specified output file:

```
$ python3 genefunction.py [-h] [-gene_name GENE_NAME] [-input_file INPUT_FILE] [-output_file OUTPUT_FILE]                       

optional arguments:
  -h, --help            show this help message and exit
  -gene_name GENE_NAME  A single gene whose function will be printed in the
                        terminal (no input/output file required).
  -input_file INPUT_FILE
                        A plain text file of newline-delimited gene names.
  -output_file OUTPUT_FILE
                        The output plain text file of newline-delimited gene
                        names with their functions.
```

## Example
```
$ python3 genefunction.py -gene_name TP53
TP53: component of the EKC/KEOPS complex that is required for the formation of a threonylcarbamoyl 
group on adenosine at position 37 in tRNAs that read codons beginning with adenine
```
