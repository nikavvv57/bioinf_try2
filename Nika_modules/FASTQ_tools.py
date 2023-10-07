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


def filter_fastq(seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0):
    """
    Takes fastq sequences(dict), where is the key - sequence name(str) and the value is a tuple of two strings:
    sequence and quality. Returns the dictionary, consisting only of those sequences that passed all the conditions.
    :param seqs: dict
    :param gc_bounds: int,float or tuple. If it's an int or float, then the function then the function discards reads
    with a GC composition lower than this value.If it's a tuple, then function will filter sequences in this range. 
    :param length_bounds: int,float or tuple. If it's an int or float, then the function then the function discards reads
    with a length lower than this value.If it's a tuple, then function will filter sequences in this range. 
    :param quality_threshold: int. All reads with quality below the value will be discarded.
    :return: dict.
    """
    filtered_seqs = {}
    for name, (seq, qual) in seqs.items():
        if filter_gc(seq, gc_bounds) and filter_len(seq, length_bounds) and filter_qual(qual, quality_threshold):
            filtered_seqs[name] = seq
    return filtered_seqs


