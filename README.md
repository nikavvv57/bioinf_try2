# ODBS(Obtaining Data from Biological Sequences)

## Overview:
The idea of this program is to make life easier for experimenters and bioinformaticians working with **fastq sequences**, **amino acid sequences** (both long and short) and **DNA and/or RNA sequences**. From the main script you can call any function from the three modules(FASTQ_tools, dna_rna_tools, amino_acids_tools). The ODBS tool will allow you to quickly obtain various characteristics of dna/rna/amino acid/fastq sequences. 

## Instruction: 
### Amino_acids 
The function **‘run_aminoacid_seq’** takes 1 amino acid sequence as input (both in three-letter and one-letter form, *str*), next you need to specify named arguments: function, record_type.
- The argument **sequence** entirely must be written in capital or lowercase letters without spaces or commas. Output saves register
- The argument **“function “** must be passed a string with the name of the action (see possible actions below) that needs to be performed. 
- The argument **“record_type”** indicates the form in which you present your sequence. If you are using a three-letter sequence, specify **“record_type= 3”**. If you are using a one -letter sequence, specify **“record_type= 1”** (by default). 

#### Example:
```
run_aminoacid_seq('ALAGLNGLU', function = 'count', record_type = 3)
```          
```
run_aminoacid_seq('glyvalala', function = 'count', record_type = 3)
```     

Also, if necessary specify a named argument **percent=True** (default False) for actions: determine_charge, determine_polarity (Look in the description of actions).

#### Example:
```
run_aminoacid_seq('LLYdD', function = 'determine_charge', record_type = 1, percent=True)
``` 

### FASTQ 
The function **‘filter_fastq’** takes fastq sequences(**seqs**, *dict*), where the key is - sequence name(*str*) and the value is a tuple of two strings: sequence(*str*) and quality(*str*). **‘filter_fastq’** filters reads according to different parameters and returns the dictionary, consisting only of those sequences that passed all the conditions. You need to specify named arguments: gc_bounds, length_bounds, quality_threshold.
- **gc_bounds**: *int, float or tuple*. If it's an int or float, then the function then the function discards reads with a GC composition lower than this value. If it's a tuple, then function will filter sequences in this range. 
- **length_bounds**: *int, float or tuple*. If it's an int or float, then the function discards reads with a length lower than this value. If it's a tuple, then function will filter sequences in this range. 
- **quality_threshold**: *int*. All reads with quality below the value will be discarded.
- **return**: *dict*. Consists of filtered sequences.

#### Example:
```
# dict of fastq files
seqs = {
    "Seq1": ("ATKLATGAT", "BBBBBBBB@"), 
    "Seq2": ("GCTAGCTAGCTA", "DDDDDDDDDDDD"),
    "Seq3": ("ATCGATCGATCG", "BBBBBBBBBBBB"),
}
```
#### Example:
```
filter_fastq(seqs, gc_bounds=(30, 70), length_bounds=(10, 15), quality_threshold=20)
```          
```
filter_fastq(seqs, gc_bounds=44.4, length_bounds=44.4, quality_threshold=20)
```    
### RNA or DNA tools
The **‘run_dna_rna_tools’** function takes as input an arbitrary number of arguments with DNA or/and RNA sequences (*str*), as well as the name of the **procedure** to be executed (which is the last argument, *str*). If one sequence is given as input, a string with the result is returned. If several sequences are supplied, a list of strings is returned. Output saves register. 
#### Example:
```
sequences = ['ATGcG', 'ATGUACCu', 'complement']
run_dna_rna_tools(*sequences)
```
## Troubleshooting:
### Amino_acids
The program works with 20 proteinogenic amino acids {'G', 'A', 'V', 'L', 'I', 'M', 'P', 'F', 'W', 'S', 'T', 'N', 'Q', 'Y', 'C', 'K', 'R', 'H', 'D', 'E'}. If you use a symbol that does not represent an amino acid, or wrong record type, the program will generate an **error** where you can see the first wrong symbol. 
- correct function launch:    
```run_aminoacid_seq('ALALEUILE', function = 'count', record_type = 3)```  
- incorrect function launch:      
```run_aminoacid_seq('ALAmIle', function = 'count', record_type = 3)```

### FASTQ
The program works only with DNA nucleotides {'A', 'T', 'G', 'C'}. If you use a symbol that does not represent in dna_aphabet, the program will generate an **error** where you can see the first wrong symbol.

### RNA or DNA tools
- The program works only with DNA or RNA nucleotides dna = {'A', 'T', 'G', 'C'}, rna = {'A', 'U', 'G', 'C'} . If you use a symbol that is not represent in dna_aphabet or rna_alphabet the program will generate an **error**.
- If the RNA sequence is sent to the ```transcribe``` function, the program will generate an **error**. Since this function only accepts DNA as input.

## Possible actions:
### Amino_acids
1. **translate** - Translation of a one-letter amino acid sequence into a three-letter one (for better visual perception), and the reverse operation. Output: str
2. **count** - obtaining the length of the amino acid sequence. Output: int
3. **count_possible_number_of_disulfide_bonds** - counting the number of possible combinations of two different cysteines to form a disulfide bond. Output: int
4. **count_molecular_weight** - calculating the molecular weight of a protein. Output: int
5. **determine_charge** - counting the number of positive, negative and neutral amino acids in a protein. To get the output in percent, specify percent=True. Output: dict
6. **determine_polarity** - counting hybrophobic and hydrophilic amino acids in a protein. To get the output in percent, specify percent=True. Output: dict
7. **convert_amino_acid_seq_to_dna** - convert an amino acid sequence to the most likely DNA sequence. Output: str
8. **summary** - a summary of all information about the sequence (the result of executing all functions). Output: dict
#### Example:

```
run_aminoacid_seq('LLYdD', function = 'translate', record_type = 1)
```   
```
run_aminoacid_seq('ALAGLYALA', function = 'translate', record_type = 3)
```          
```
run_aminoacid_seq('LLYdD', function = 'count', record_type = 1)
```          
```
run_aminoacid_seq('ALAGLYALA', function = 'count_protein_length', record_type = 3)
```          
```
run_aminoacid_seq('LLYdD', function = 'summary', record_type = 1)
```          
```
run_aminoacid_seq('alaglyala', function = 'determine_charge', record_type = 3, percent=True)
```   

```
run_aminoacid_seq('LLYdD', function = 'determine_charge', record_type = 1, percent=True)
```
```
run_aminoacid_seq('alaglyala', function = 'determine_charge', record_type = 3)
```          
```
run_aminoacid_seq('LLYd', function = 'determine_charge', record_type = 1)
```         
### FASTQ
The only action is filtering fastq sequences according to this parameters:
1. ```filter_gc```
2. ```filter_len```
3. ```filter_qual```
### RNA or DNA tools
- ```transcribe``` — print the transcribed sequence. If an RNA sequence is supplied as input, then the program is interrupted and generates an error
- ```reverse``` — print the reversed sequence
- ```complement``` — print the complementary sequence
- ```reverse_complement``` — print the reverse complementary sequence
- ```gc_content``` — calculation of GC content in sequence

#### Example:
```
sequences = ['ATGcG', 'ATGUACCu', 'reverse']
run_dna_rna_tools(*sequences)
```
