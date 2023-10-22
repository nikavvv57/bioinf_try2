# Bio_files_processor
## Overview:
There are two functions:
1. convert_multiline_fasta_to_oneline.
2. select_genes_from_gbk_to_fasta.

## Instruction: 
### convert_multiline_fasta_to_oneline 
The function takes multiline fasta file as an input and convert the multiline sequences in one line.
- **input_fasta**: *str*. Path to multiline fasta file.
- **output_fasta**: *str*. Name of new one_line fasta file. This parameter isn't necessary. If there isn't output_fasta function will write a new file with the input file's name. Adds ".fasta" to the end of file name, if there isn't.
- **return**: *file*. Consists of one line sequences.

#### Example:
```
if __name__ == "__main__":
    input_fasta = "/Users/veronikavadehina/Downloads/example_multiline_fasta.fasta"
    convert_multiline_fasta_to_oneline(input_fasta)
```

### select_genes_from_gbk_to_fasta:
This function takes a gbk file, parse it and find name of our gene around which we want to see genes of interest. Returns file with gene of interest names and their translation sequences.
- **input_gbk**: *str*. Path to gbk file.
- **gene**: *str*. Gene name around which we want to see genes of interest.
- **n_before**: *int*. Numbers of genes of interest before the main gene. By default = 1. Strict > 0. 
- **n_after**: *int*. Numbers of genes of interest after the main gene. By default = 1. Strict > 0.
- **output_gbk_genes_fasta**: *str*. Name of output file. If there isn't name, function will write a new file with the input file's name. Adds ".fasta" to the end of file name, if there isn't.
- **return**: file.
#### Example:

```
if __name__ == "__main__":
    input_gbk = "/Users/veronikavadehina/Downloads/example.gbk"
    gene = 'dtpD'
    select_genes_from_gbk_to_fasta(input_gbk, gene, n_before=2, n_after=1, output_gbk_genes_fasta='gene_of_interest')
```   
