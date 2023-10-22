import os


def parse_gbk(input_gbk):
    """
    This function takes gbk file and extracts the gene names along with their corresponding translational sequences. 
    :param input_gbk: str. path to gbk file. 
    :return: list[list[str]].
    """
    with open(input_gbk, 'r') as gbk_file:
        lines = gbk_file.readlines()

    cds = []
    for i, line in enumerate(lines):
        if 'CDS' in line:
            if '/gene' in lines[i + 1]:
                gen_name = lines[i + 1].split('gene=')[-1].split()[0].strip('\"')
            else:
                gen_name = 'noname'
            translation = ''
            sub_lines = lines[i:len(lines) + 1]
            for j, lin in enumerate(sub_lines):
                if '/translation' in lin:
                    start_trans = j + 1
                    translation += lin.split('/translation=')[-1].strip().strip('"\n')
                    while '"' not in sub_lines[start_trans]:
                        translation += sub_lines[start_trans].strip().strip('"\n')
                        start_trans += 1
                    translation += sub_lines[start_trans].strip().strip('"\n')
                    break
            cds.append([gen_name, translation])
    return cds

