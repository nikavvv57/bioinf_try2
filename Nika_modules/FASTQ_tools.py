dna_nucleotide_alphabet = {'A', 'T', 'G', 'C'}


def filter_gc(seq, gc_bounds):
    for nuc in seq:
        if nuc not in dna_nucleotide_alphabet:
            raise ValueError(f"{nuc} this is not a nucleotide in your sequence")
    gc_content = (seq.count('G') + seq.count('C')) / len(seq) * 100
    if isinstance(gc_bounds, int):
        return gc_content <= gc_bounds
    elif isinstance(gc_bounds, float):
        return gc_content <= gc_bounds
    elif isinstance(gc_bounds, tuple):
        return gc_bounds[0] <= gc_content <= gc_bounds[1]


def filter_len(seq, length_bounds):
    if isinstance(length_bounds, int):
        return len(seq) <= length_bounds
    elif isinstance(length_bounds, float):
        return len(seq) <= length_bounds
    elif isinstance(length_bounds, tuple):
        return length_bounds[0] <= len(seq) <= length_bounds[1]


def filter_qual(qual, quality_threshold):
    avg_quality = sum(ord(q) - 33 for q in qual) / len(qual)
    return avg_quality >= quality_threshold
