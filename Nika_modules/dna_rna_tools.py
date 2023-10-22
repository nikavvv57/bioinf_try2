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
