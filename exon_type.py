# GENCODE augmented DB header
# 0_chr1_11869_12227_+, 1_ENSG00000223972.4, 2_ENST00000456328.2, 3_pseudogene, 4_DDX11L1, 5_1, 6_MISSING_START,
# 7_MISSING_STOP, 8_12010, 9_13670
from __future__ import division
augmentedDBpath = "C:/Users/mjk140030/Desktop/GTEx/toyFile.txt"
outputPath = augmentedDBpath[:-4] + "_ver2.txt"
with open(augmentedDBpath) as f:
    length = sum(1 for _ in f)
print("Determining Exon Type...")


def get_type(line):
    line = map(int, line)
    exstart = line[0]
    exstop = line[1]
    TSS = line[2]
    start_codon = line[3]
    stop_codon = line[4]
    TES = line[5]
    vector = [0, 0, 0, 0, 0]
    if exstart < TSS:
        print(line)
    elif exstop > TES:
        print(line)
    else:
        if exstart <= start_codon - 3:  # ex start before start codon
            vector[0] = 1
            if exstop <= start_codon - 3:  # exstop before start codon
                vector[0] = 1
                return vector
            elif start_codon - 3 <= exstop <= start_codon:  # ex stop in start codon
                vector[1] = 1
                return vector
            elif start_codon < exstop <= stop_codon - 3:  # exstop in CDS
                vector[2] = 1
                print(vec_to_str(vector))
                quit()
                return vector
            elif stop_codon - 3 < exstop <= stop_codon:  # exstop in stop codon
                vector[3] = 1
                return vector
            elif exstop >= stop_codon:  # exstop stop after stop codon
                vector[4] = 1
                return vector
        elif stop_codon - 3 <= exstart <= start_codon:  # ex start in start codon
            vector[1] = 1
            if start_codon - 3 <= exstop <= start_codon:
                vector[1] = 1
                return vector
            elif start_codon < exstop <= stop_codon - 3:
                vector[2] = 1
                return vector
            elif stop_codon - 3 < exstop <= stop_codon:
                vector[3] = 1
                return vector
            elif exstop >= stop_codon:
                vector[4] = 1
                return vector
        elif start_codon + 1 <= exstart <= stop_codon - 3:  # ex start in cds
            vector[2] = 1
            if start_codon < exstop <= stop_codon - 3:
                vector[2] = 1
                return vector
            elif stop_codon - 3 < exstop <= stop_codon:
                vector[3] = 1
                return vector
            elif exstop >= stop_codon:
                vector[4] = 1
                return vector
        elif stop_codon - 3 < exstart <= stop_codon:  # ex start in stop codon
            vector[3] = 1
            if stop_codon - 3 < exstop <= stop_codon:
                vector[3] = 1
                return vector
            elif exstop >= stop_codon:
                vector[4] = 1
                return vector
        else:  # ex start after stop codon
            if exstop <= TES:
                vector[4] = 1
                return vector


def vec_to_str(vector):
    vector = map(int, list(vector))
    if len(vector) < 5:
        print("VECTOR ERROR: INCORRECT LENGTH")
        return "ERROR"
    if vector.count(1) > 2:
        print("VECTOR ERROR: MORE THAN 2 HITS!")
        return "ERROR"
    if vector.count(1) == 2 :  # more than one exon type, have to fill in between the ones
        first1 = vector.index(1)
        last1 = vector[first1+1:].index(1) + (5-len(vector[first1+1:]))
        for x in range(first1, last1+1):
            vector[x] = 1
    exonType = ""
    if vector[0] == 1:  # must be 5UTR
        exonType += "5UTR,"
    if vector[1] == 1:  # must be SAC
        exonType += "SAC,"
    if vector[2] == 1:
        exonType += "CDS,"
    if vector[3] == 1:
        exonType += "SPC,"
    if vector[4] == 1:
        exonType += "3UTR,"
    return exonType[:-1]  # remove the extra comma at the end


with open(augmentedDBpath) as f:
    line_count = 1
    for line in f:  # for each line determine exon type
        if line_count % 1200 == 0:
            val = (line_count / length) * 100
            print(str(round(val, 2)) + "%")
        elements = line.split()
        id_elem = elements[0].split('_')
        strand = id_elem[3]
        exon_start = int(id_elem[1])
        exon_stop = int(id_elem[2])
        gene_5_end = int(elements[8])
        gene_3_end = int(elements[9])
        start_codon_end_pos = elements[6]
        stop_codon_end_pos = elements[7]
        if ',' in start_codon_end_pos:  # number of start is 2
            first_start = start_codon_end_pos.split(",")[0]
            second_start = start_codon_end_pos.split(",")[1]
            if strand == "+":
                new_line1 = [exon_start, exon_stop, gene_5_end, first_start, stop_codon_end_pos, gene_3_end]
                new_line2 = [exon_start, exon_stop, gene_5_end, second_start, stop_codon_end_pos, gene_3_end]
                exon_type1 = get_type(new_line1)
                exon_type2 = get_type(new_line2)
                exon_type = vec_to_str(exon_type1)
                exon_type += "_"+vec_to_str(exon_type2)
                new_line = line[:-1] + "\t" + exon_type + "\n"
                with open(outputPath, 'a') as f:
                    f.write(new_line)
                continue
            else:
                new_line1 = [exon_start, exon_stop, gene_5_end, stop_codon_end_pos, first_start, gene_3_end]
                new_line2 = [exon_start, exon_stop, gene_5_end, stop_codon_end_pos, second_start, gene_3_end]
                exon_type1 = get_type(new_line1)
                exon_type2 = get_type(new_line2)
                exon_type = vec_to_str(exon_type1)
                exon_type += "_" + vec_to_str(exon_type2)
                new_line = line[:-1] + "\t" + exon_type + "\n"
                with open(outputPath, 'a') as f:
                    f.write(new_line)
                continue
        elif ',' in stop_codon_end_pos:  # number of stops is 2
            first_stop = stop_codon_end_pos.split(",")[0]
            second_stop = stop_codon_end_pos.split(","[1])
            if strand == "+":
                new_line1 = [exon_start, exon_stop, gene_5_end, start_codon_end_pos, first_stop, gene_3_end]
                new_line2 = [exon_start, exon_stop, gene_5_end, start_codon_end_pos, second_stop, gene_3_end]
                exon_type1 = get_type(new_line1)
                exon_type2 = get_type(new_line2)
                exon_type = vec_to_str(exon_type1)
                exon_type += "_" + vec_to_str(exon_type2)
                new_line = line[:-1] + "\t" + exon_type + "\n"
                with open(outputPath, 'a') as f:
                    f.write(new_line)
                continue
            else:
                new_line1 = [exon_start, exon_stop, gene_5_end, first_stop, start_codon_end_pos, gene_3_end]
                new_line2 = [exon_start, exon_stop, gene_5_end, second_stop, start_codon_end_pos, gene_3_end]
                exon_type1 = get_type(new_line1)
                exon_type2 = get_type(new_line2)
                # turn exon type vector into string here
                new_line = line[:-1] + "\t" + exon_type + "\n"
                with open(outputPath, 'a') as f:
                    f.write(new_line)
                continue
        elif start_codon_end_pos == "NONCODING":  # no start
            start_codon_end_pos = gene_5_end
            if strand == "+":
                new_line = [exon_start, exon_stop, gene_5_end, start_codon_end_pos, stop_codon_end_pos, gene_3_end]
                exon_type = get_type(new_line)
            else:
                new_line = [exon_start, exon_stop, gene_5_end, stop_codon_end_pos, start_codon_end_pos, gene_3_end]
                exon_type = get_type(new_line)
            # turn exon type vector into string here
            new_line = line[:-1] + "\t" + exon_type + "\n"
            with open(outputPath, 'a') as f:
                f.write(new_line)
            continue
        elif stop_codon_end_pos == "MISSING_STOP":  # no stop
            stop_codon_end_pos = gene_3_end
            if strand == "+":
                new_line = [exon_start, exon_stop, gene_5_end, start_codon_end_pos, stop_codon_end_pos, gene_3_end]
                exon_type = get_type(new_line)
            else:
                new_line = [exon_start, exon_stop, gene_5_end, stop_codon_end_pos, start_codon_end_pos, gene_3_end]
                exon_type = get_type(new_line)
            # turn exon type vector into string here
            new_line = line[:-1] + "\t" + exon_type + "\n"
            with open(outputPath, 'a') as f:
                f.write(new_line)
            continue
        elif stop_codon_end_pos == "MISSING_STOP" and start_codon_end_pos == "NONCODING":  # no stop AND no start
            stop_codon_end_pos = gene_3_end
            start_codon_end_pos = gene_5_end
            if strand == "+":
                new_line = [exon_start, exon_stop, gene_5_end, start_codon_end_pos, stop_codon_end_pos, gene_3_end]
                exon_type = get_type(new_line)
            else:
                new_line = [exon_start, exon_stop, gene_5_end, stop_codon_end_pos, start_codon_end_pos, gene_3_end]
                exon_type = get_type(new_line)
            exon_type = vec_to_str(exon_type1)
            exon_type += "_" + vec_to_str(exon_type2)
            new_line = line[:-1] + "\t" + exon_type + "\n"
            with open(outputPath, 'a') as f:
                f.write(new_line)
            continue
        else:  # only 1 start and only 1 stop present
            if strand == "+":
                new_line = [exon_start, exon_stop, gene_5_end, start_codon_end_pos, stop_codon_end_pos, gene_3_end]
                exon_type = vec_to_str(get_type(new_line))
            else:
                new_line = [exon_start, exon_stop, gene_5_end, stop_codon_end_pos, start_codon_end_pos, gene_3_end]
                exon_type = vec_to_str(get_type(new_line))
            new_line = line[:-1] + "\t" + exon_type + "\n"
            with open(outputPath, 'a') as f:
                f.write(new_line)

        line_count += 1

print("Done. Output written to: " + outputPath)
