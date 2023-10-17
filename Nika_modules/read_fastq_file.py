def read_fastq(input_path):
    """
    This function read fastq file and returns dictionary, where key is a name of sequence, and value is a tuple(
    sequence, comment and quality)
    :param input_path: str. path to fastq file.
    :return: dict.
    """
    fastq_dict = {}

    with open(input_path, 'r') as fastq_file:
        lines = fastq_file.readlines()

        for i in range(0, len(lines), 4):
            name_seq = lines[i].strip()
            sequence = lines[i + 1].strip()
            comment = lines[i + 2].strip()
            quality = lines[i + 3].strip()

            fastq_dict[name_seq] = (sequence, comment, quality)

    return fastq_dict
