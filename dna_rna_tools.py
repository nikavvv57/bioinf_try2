def gc_content(sequence):
    # GC content count
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    total_bases = len(sequence)
    gc_content_percent = (gc_count / total_bases) * 100
    return gc_content_percent


def transcribe(sequence):
    # DNA to RNA transcription
    rna_seq = ""
    for base in sequence:
        if base == 'U':
            raise ValueError("Input sequence contains uracil. This function is for DNA sequences only.")
        elif base == "u":
            raise ValueError("Input sequence contains uracil. This function is for DNA sequences only.")
        elif base == "T":
            rna_seq += "U"
        elif base == "t":
            rna_seq += "u"
        else:
            rna_seq += base
    return rna_seq


def reverse(sequence):
    # Reverse DNA or RNA sequence
    return sequence[::-1]


def complement(sequence):
    # Build complement sequence
    complement_dict_dna = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                           'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    complement_dict_rna = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C',
                           'a': 'u', 'u': 'a', 'c': 'g', 'g': 'c'}
    if all(base in dna_nucleotide_alphabet for base in sequence):
        complement_seq = ''.join(complement_dict_dna.get(base, base) for base in sequence)
    else:
        complement_seq = ''.join(complement_dict_rna.get(base, base) for base in sequence)
    return complement_seq


def reverse_complement(sequence):
    # Build reverse complement sequence
    complement_dict_dna = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                           'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    complement_dict_rna = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C',
                           'a': 'u', 'u': 'a', 'c': 'g', 'g': 'c'}
    if all(base in dna_nucleotide_alphabet for base in sequence):
        reversed_complement_seq = ''.join(complement_dict_dna.get(base, base) for base in reversed(sequence))
    else:
        reversed_complement_seq = ''.join(complement_dict_rna.get(base, base) for base in reversed(sequence))
    return reversed_complement_seq


dna_nucleotide_alphabet = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c'}
rna_nucleotide_alphabet = {'A', 'U', 'G', 'C', 'a', 'u', 'g', 'c'}


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
    :param sequences: str, can be more than 1 seq, nd only 1 procedure.
    :return: str or list. If one sequence is given as input, a string with the result is returned.
    If several sequences are supplied, a list of strings is returned.
    """
    procedure = sequences[-1]  # Unpack procedure
    sequences = sequences[:-1]

    results = []

    for sequence in sequences:
        for nucleotide in sequence:
            if nucleotide in dna_nucleotide_alphabet:
                continue
            elif nucleotide in rna_nucleotide_alphabet:
                continue
            else:
                raise ValueError('this is not a nucleotide sequence')
        if procedure == 'transcribe':
            result = transcribe(sequence)
        elif procedure == 'reverse':
            result = reverse(sequence)
        elif procedure == 'reverse_complement':
            result = reverse_complement(sequence)
        elif procedure == 'complement':
            result = complement(sequence)
        elif procedure == 'gc_content':
            result = gc_content(sequence)
        else:
            result = "Invalid procedure"
        results.append(result)

    if len(results) == 1:
        return results[0]
    else:
        return results

