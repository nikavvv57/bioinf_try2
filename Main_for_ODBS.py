from Nika_modules import FASTQ_tools, dna_rna_tools, amino_acids_tools
from Nika_modules.amino_acids_tools import single_letter_alphabet, three_letter_alphabet


def filter_fastq(seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0):
    """
    Takes fastq sequences(dict), where is the key - sequence name(str), and the value is a tuple of two strings:
    sequence and quality. Returns the dictionary, consisting only of those sequences that passed all the conditions.
    :param seqs: dict
    :param gc_bounds: int, float or tuple. If it's an int or float, then the function discards reads
    with a GC composition lower than this value.
    :param length_bounds: int, float or tuple. If it's an int or float, then the function discards reads
    with a length lower than this value.
    :param quality_threshold: int. All reads with quality below the value will be discarded.
    :return: dict.
    """
    filtered_seqs = {}
    for name, (seq, qual) in seqs.items():
        if FASTQ_tools.filter_gc(seq, gc_bounds) and FASTQ_tools.filter_len(seq, length_bounds) \
                and FASTQ_tools.filter_qual(qual, quality_threshold):
            filtered_seqs[name] = seq
    return filtered_seqs


def run_dna_rna_tools(*sequences):
    """
    This run_dna_rna_tools function takes as input an arbitrary number of arguments with DNA or/and RNA sequences,
    as well as the name of the procedure to be executed (which is the last argument, str).
    Possible procedures performed by the program:
    - transcribe — print the transcribed sequence. If an RNA sequence is supplied as input,
    then the program is interrupted and generates an error.
    - reverse — print the reversed sequence
    - complement — print the complementary sequence
    - reverse_complement — print the reverse complementary sequence
    - gc_content—calculation of GC content in sequence,
    :param sequences: str, can be more than 1 seq, only 1 procedure.
    :return: str or list. If one sequence is given as input, a string with the result is returned.
    If several sequences are supplied, a list of strings is returned.
    """
    procedure = sequences[-1]  # Unpack procedure
    sequences = sequences[:-1]

    results = []
    dna_nucleotide_alphabet = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c'}
    rna_nucleotide_alphabet = {'A', 'U', 'G', 'C', 'a', 'u', 'g', 'c'}
    for sequence in sequences:
        for nucleotide in sequence:
            if nucleotide in dna_nucleotide_alphabet:
                continue
            elif nucleotide in rna_nucleotide_alphabet:
                continue
            else:
                raise ValueError('this is not a nucleotide sequence')
        if procedure == 'transcribe':
            result = dna_rna_tools.transcribe(sequence)
        elif procedure == 'reverse':
            result = dna_rna_tools.reverse(sequence)
        elif procedure == 'reverse_complement':
            result = dna_rna_tools.reverse_complement(sequence)
        elif procedure == 'complement':
            result = dna_rna_tools.complement(sequence)
        elif procedure == 'gc_content':
            result = dna_rna_tools.gc_content(sequence)
        else:
            result = "Invalid procedure"
        results.append(result)

    if len(results) == 1:
        return results[0]
    else:
        return results


def run_aminoacid_seq(sequence: str, function: str = 'summary', record_type: int = 1, percent: bool = False):
    """
    Performs the following list of operations:
    count - counts number of amino acid
    translate - converts one record type to another
    determine_charge - counts number or percent of amino acid with different charges
    determine_polarity - counts number or percent of amino acid with different polarity
    convert_amino_acid_seq_to_dna - takes an amino acid sequence as input
    and returns the optimal DNA sequence for E.coli
    count_possible_number_of_disulfide_bonds - counting the number of possible combinations of two different cysteines
    to form a disulfide bond
    count_molecular_weight - takes an amino acid sequence as input and returns the molecular weight of the protein
    summary - returns results of all functions (default)

    Arguments:
    - sequence:str - sequence for function
    - function:str - name of the function you need to perform. You can use: 'count', 'translate', 'summary'(default)
    - record_type:int - record type of your sequence. You can use 1(default) for single letter type or
    3 for three letter type
    - percent: bool - for determine_charge and determine_polarity shows result in percent (True) or in number (False).
    Default - False

    Return:
    - count - int
    - translate - str
    - determine_charge - dict
    - determine_polarity - dict
    - convert_amino_acid_seq_to_dna - str
    - count_possible_number_of_disulfide_bonds - int
    - count_molecular_weight - int
    - summary - dict

    """
    if record_type == 1:
        for amino_acid in sequence:
            if amino_acid not in single_letter_alphabet:
                raise ValueError(f"{amino_acid} in your sequence is not amino acid")
    elif record_type == 3:
        tuple_sequence = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
        for amino_acid in tuple_sequence:
            if amino_acid not in three_letter_alphabet:
                raise ValueError(f"{amino_acid} in your sequence is not amino acid")
        sequence = ''.join(tuple_sequence)
    if function == 'count':
        return amino_acids_tools.count(sequence, record_type)
    elif function == 'translate':
        return amino_acids_tools.translate(sequence, record_type)
    elif function == 'summary':
        return amino_acids_tools.summary(sequence, record_type, percent)
    elif function == 'determine_polarity':
        if record_type == 3:
            return amino_acids_tools.determine_polarity(amino_acids_tools.translate(sequence, record_type), percent)
        else:
            return amino_acids_tools.determine_polarity(sequence, percent)
    elif function == 'determine_charge':
        if record_type == 3:
            return amino_acids_tools.determine_charge(amino_acids_tools.translate(sequence, record_type), percent)
        else:
            return amino_acids_tools.determine_charge(sequence, percent)
    elif function == 'count_possible_number_of_disulfide_bonds':
        if record_type == 3:
            return amino_acids_tools.count_possible_number_of_disulfide_bonds(amino_acids_tools.translate(sequence,
                                                                                                          record_type))
        else:
            return amino_acids_tools.count_possible_number_of_disulfide_bonds(sequence)
    elif function == 'count_molecular_weight':
        if record_type == 3:
            return amino_acids_tools.count_molecular_weight(amino_acids_tools.translate(sequence, record_type))
        else:
            return amino_acids_tools.count_molecular_weight(sequence)
    elif function == 'convert_amino_acid_seq_to_dna':
        if record_type == 3:
            return amino_acids_tools.convert_amino_acid_seq_to_dna(amino_acids_tools.translate(sequence, record_type))
        else:
            return amino_acids_tools.convert_amino_acid_seq_to_dna(sequence)
