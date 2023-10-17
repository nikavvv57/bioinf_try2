import os
from Nika_modules import parser_gbk


def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None):
    """
    This function takes multiline fasta file as an input and convert the multiline sequences in one line.
    :param input_fasta: str. path to multiline fasta file
    :param output_fasta: str. name of new one line fasta file. This parameter isn't necessary. 
    :return: file. 
    """
    try:
        with open(input_fasta, 'r') as input_file:
            fasta_data = input_file.read()
    except FileNotFoundError:
        return "Input fasta file not found."

    oneline_fasta = []
    fasta_parts = fasta_data.split('\n>')

    for part in fasta_parts:
        if part:
            lines = part.strip().split('\n')
            header = lines[0]
            sequence = ''.join(lines[1:])
            oneline_fasta.append(f'>{header}\n{sequence}')

    oneline_fasta = '\n'.join(oneline_fasta)

    if not os.path.exists("/Users/veronikavadehina/Desktop/Bioinf "
                          "23:24/Python_Nika_GitHub/bioinf_for_py/fasta_oneline_rewrite_files"):
        os.makedirs("/Users/veronikavadehina/Desktop/Bioinf "
                    "23:24/Python_Nika_GitHub/bioinf_for_py/fasta_oneline_rewrite_files")

    if output_fasta is None:
        output_fasta = os.path.basename(input_fasta)

    if not output_fasta.endswith(".fasta"):
        output_fasta += ".fasta"

    output_file_fasta_path = os.path.join("/Users/veronikavadehina/Desktop/Bioinf "
                                          "23:24/Python_Nika_GitHub/bioinf_for_py/fasta_oneline_rewrite_files",
                                          output_fasta)
    if output_fasta:
        try:
            with open(output_file_fasta_path, 'w') as output_file:
                output_file.write(oneline_fasta)
            print(f"Fasta in one line saved to {output_file_fasta_path}")
        except Exception as e:
            print(f"Error: {e}")
            print("Failed to save one line fasta.")

    return oneline_fasta


def select_genes_from_gbk_to_fasta(input_gbk, gene, n_before=1, n_after=1, output_gbk_genes_fasta=None):
    """
    This function takes a gbk file, parse it and find name of our gene around which we want to see genes of interest.
    Returns file with gene of interest names and their translation sequences. 
    :param input_gbk: str. Path to gbk file. 
    :param gene: str. Gene name around which we want to see genes of interest.  
    :param n_before: int. Numbers of genes of interest before the main gene. 
    :param n_after: int. Numbers of genes of interest after the main gene.
    :param output_gbk_genes_fasta: str. Name of output file. 
    :return: file with gene of interest names and their translation sequences.
    """
    result = {}
    cds_list = parser_gbk.parse_gbk(input_gbk)
    for i, gen in enumerate(cds_list):
        if gene == gen[0]:
            target_index = i
            break
    for z in range(target_index - n_before, target_index):
        gene_interest_lower = cds_list[z][0]
        translation_interest_lower = cds_list[z][1]
        result[gene_interest_lower] = translation_interest_lower
    for h in range(target_index + 1, target_index + n_after + 1):
        gene_interest_upper = cds_list[h][0]
        translation_interest_upper = cds_list[h][1]
        result[gene_interest_upper] = translation_interest_upper

    if not os.path.exists("/Users/veronikavadehina/Desktop/Bioinf "
                          "23:24/Python_Nika_GitHub/bioinf_for_py/fasta_gbk_interest_gene"):
        os.makedirs("/Users/veronikavadehina/Desktop/Bioinf "
                    "23:24/Python_Nika_GitHub/bioinf_for_py/fasta_gbk_interest_gene")

    if output_gbk_genes_fasta is None:
        output_gbk_genes_fasta = os.path.basename(input_gbk)

    if not output_gbk_genes_fasta.endswith(".fasta"):
        output_gbk_genes_fasta += ".fasta"

    path = os.path.join("/Users/veronikavadehina/Desktop/Bioinf "
                        "23:24/Python_Nika_GitHub/bioinf_for_py/fasta_gbk_interest_gene",
                        output_gbk_genes_fasta)
    if output_gbk_genes_fasta:
        try:
            with open(path, 'w') as output_file:
                for gen_interest, gene_interest_translation in result.items():
                    output_file.write(f'>{gen_interest}\n{gene_interest_translation}\n')
            print(f"Fasta with interested genes saved to {path}")
        except Exception as e:
            print(f"Error: {e}")
            print("Failed to save fasta.")
    return result

