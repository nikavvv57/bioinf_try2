import os


def save_filtered_fastq(filtered_seqs, input_path, output_filename=None):
    """
    This function should save filtered reads in new fastq file. It creates a folder with the name -
    fastq_filtrator_results, if there isn't such directory in your laptop.
    :param filtered_seqs: dict. From filter_fastq tool.
    :param input_path: str. Path to non filtered fastq file with reads, it's necessary for
    Main_for_ODBS.py.
    :param output_filename: str. If there is no output_filename parameter function saves new fastq file with the name
    of non filtered file. And adds .fastq to the end of file name, if there isn't.
    :return: fastq file. 
    """
    if not os.path.exists("fastq_filtrator_results"):
        os.makedirs("fastq_filtrator_results")

    if output_filename is None:
        output_filename = os.path.basename(input_path)

    if not output_filename.endswith(".fastq"):
        output_filename += ".fastq"

    output_file_path = os.path.join("fastq_filtrator_results", output_filename)

    try:
        with open(output_file_path, "w") as output_file:
            for name, (seq, comment, qual) in filtered_seqs.items():
                output_file.write(f"{name}\n{seq}\n{comment}\n{qual}\n")

        print(f"Filtered data saved to {output_file_path}")
    except Exception as e:
        print(f"Error: {e}")
        print("Failed to save the filtered data.")

