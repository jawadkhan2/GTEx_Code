# GENCODE augmented DB header
# 0_chr1_11869_12227_+, 1_ENSG00000223972.4, 2_ENST00000456328.2, 3_pseudogene, 4_DDX11L1, 5_1, 6_MISSING_START,
# 7_MISSING_STOP, 8_12010, 9_13670
from __future__ import division
import os
import smtplib
from email.mime.text import MIMEText

sender = "script.info1@gmail.com"
sender_pass = "info4script"


def notify_user(recipient, notification):
    session = smtplib.SMTP('smtp.gmail.com', 587)
    session.starttls()
    msg = MIMEText(notification)
    msg['Subject'] = "Script Update"
    msg['From'] = sender
    msg['To'] = recipient
    session.login(sender, sender_pass)
    session.sendmail(sender, recipient, msg.as_string())
    session.quit()

body = """
Exon Type Script finished.
"""

augmentedDBpath = "C:/Users/mjk140030/Desktop/GTEx/v19_gencode_table.txt"
outputPath = augmentedDBpath[:-4] + "_ver2.txt"
if os.path.isfile(outputPath):
    with open(outputPath, 'w') as a:
        print("Previous file was cleared.")
with open(augmentedDBpath) as f:
    length = sum(1 for _ in f)
print("Determining Exon Type...")
coding_gene_type = [ "G_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_V_gene", "TR_C_gene", "TR_J_gene",
                     "TR_V_gene", "TR_D_gene", "protein_coding" ]

def get_type(ln_list):
    ln_list = map(int, ln_list)
    exstart = ln_list[0]
    exstop = ln_list[1]
    tss = ln_list[2]  # Transcription start site
    start_codon = ln_list[3]
    stop_codon = ln_list[4]
    tes = ln_list[5]  # Transcription end site
    vector = [0, 0, 0, 0, 0]
    spacer = "-" * 25
    if exstart < tss:
        print("ERROR: EXON OUTSIDE OF TSS")
        print(line)
    elif exstop > tes:
        print("ERROR: EXON OUTSIDE OF TES")
        print(line)
    else:
        if exstart <= start_codon - 3:  # ex start before start codon
            vector[0] = 1
            if exstop <= start_codon - 3:  # exstop before start codon
                vector[0] = 1
            elif start_codon - 3 <= exstop <= start_codon:  # ex stop in start codon
                vector[1] = 1
            elif start_codon < exstop <= stop_codon - 3:  # exstop in CDS
                vector[2] = 1
            elif stop_codon - 3 < exstop <= stop_codon:  # exstop in stop codon
                vector[3] = 1
            elif exstop >= stop_codon:  # exstop stop after stop codon
                vector[4] = 1
        elif stop_codon - 3 <= exstart <= start_codon:  # ex start in start codon
            vector[1] = 1
            if start_codon - 3 <= exstop <= start_codon:
                vector[1] = 1
            elif start_codon < exstop <= stop_codon - 3:
                vector[2] = 1
            elif stop_codon - 3 < exstop <= stop_codon:
                vector[3] = 1
            elif exstop >= stop_codon:
                vector[4] = 1
        elif start_codon + 1 <= exstart <= stop_codon - 3:  # ex start in cds
            vector[2] = 1
            if start_codon < exstop <= stop_codon - 3:
                vector[2] = 1
            elif stop_codon - 3 < exstop <= stop_codon:
                vector[3] = 1
            elif exstop >= stop_codon:
                vector[4] = 1
        elif stop_codon - 3 < exstart <= stop_codon:  # ex start in stop codon
            vector[3] = 1
            if stop_codon - 3 < exstop <= stop_codon:
                vector[3] = 1
            elif exstop >= stop_codon:
                vector[4] = 1
        else:  # ex start after stop codon
            if exstop <= tes:
                vector[4] = 1
        #print(str(tss) + spacer + str(start_codon) + spacer + str(stop_codon) + spacer + str(tes) + "\t"
        #    + vec_to_str(vector))
        #print(str(exstart) + "\t" + str(exstop))
        #quit()
        return vector


def vec_to_str(vector):
    vector = map(int, list(vector))
    if len(vector) < 5:
        print("VECTOR ERROR: INCORRECT LENGTH")
        return "ERROR"
    if vector.count(1) > 2:
        print("VECTOR ERROR: MORE THAN 2 HITS!")
        return "ERROR"
    if vector.count(1) == 2:  # more than one exon type, have to fill in between the ones
        first1 = vector.index(1)
        last1 = vector[first1+1:].index(1) + (len(vector)-len(vector[first1+1:]))
        for x in range(first1, last1+1):
            vector[x] = 1
    et = ""  # stands for exon type not extra terrestrial
    if vector[0] == 1:  # must be 5UTR
        et += "5UTR,"
    if vector[1] == 1:  # must be SAC
        et += "SAC,"
    if vector[2] == 1:
        et += "CDS,"
    if vector[3] == 1:
        et += "SPC,"
    if vector[4] == 1:
        et += "3UTR,"
    return et[:-1]  # remove the extra comma at the end


with open(augmentedDBpath) as f:
    line_count = 1
    for line in f:  # for each line determine exon type
        if line_count % 1200 == 0:
            val = (line_count / length) * 100
            print(str(round(val, 2)) + "%")
        elements = line.split()
        id_elem = elements[0].split('_')
        strand = id_elem[3]
        if strand != "+" and strand != "-":
            print("ERROR: STRAND IS NEITHER POSITIVE OR NEGATIVE!")
            print(line)
            line_count += 1
            continue
        exon_start = int(id_elem[1])
        exon_stop = int(id_elem[2])
        gene_5_end = int(elements[8])
        gene_3_end = int(elements[9])
        start_codon_end_pos = elements[6]
        stop_codon_end_pos = elements[7]
        gene_type = elements[3]
        if ',' in start_codon_end_pos:  # number of start is 2
            first_start = start_codon_end_pos.split(",")[0]
            second_start = start_codon_end_pos.split(",")[1]
            if stop_codon_end_pos == "MISSING_STOP":
                stop_codon_end_pos = gene_3_end
            if strand == "+":
                ln_lst1 = [exon_start, exon_stop, gene_5_end, first_start, stop_codon_end_pos, gene_3_end]
                ln_lst2 = [exon_start, exon_stop, gene_5_end, second_start, stop_codon_end_pos, gene_3_end]
            else:
                ln_lst1 = [exon_start, exon_stop, gene_5_end, stop_codon_end_pos, first_start, gene_3_end]
                ln_lst2 = [exon_start, exon_stop, gene_5_end, stop_codon_end_pos, second_start, gene_3_end]
            exon_type1 = get_type(ln_lst1)
            exon_type2 = get_type(ln_lst2)
            exon_type = vec_to_str(exon_type1)
            exon_type += "_" + vec_to_str(exon_type2)
            new_line = line[:-1] + "\t" + exon_type + "\n"
            with open(outputPath, 'a') as o:
                o.write(new_line)
        elif ',' in stop_codon_end_pos:  # number of stops is 2
            first_stop = stop_codon_end_pos.split(",")[0]
            second_stop = stop_codon_end_pos.split(",")[1]
            if start_codon_end_pos == "MISSING_START":
                start_codon_end_pos = gene_5_end
            if strand == "+":
                ln_lst1 = [exon_start, exon_stop, gene_5_end, start_codon_end_pos, first_stop, gene_3_end]
                ln_lst2 = [exon_start, exon_stop, gene_5_end, start_codon_end_pos, second_stop, gene_3_end]
            else:
                ln_lst1 = [exon_start, exon_stop, gene_5_end, first_stop, start_codon_end_pos, gene_3_end]
                ln_lst2 = [exon_start, exon_stop, gene_5_end, second_stop, start_codon_end_pos, gene_3_end]
            exon_type1 = get_type(ln_lst1)
            exon_type2 = get_type(ln_lst2)
            exon_type = vec_to_str(exon_type1)
            exon_type += "_" + vec_to_str(exon_type2)
            new_line = line[:-1] + "\t" + exon_type + "\n"
            with open(outputPath, 'a') as o:
                o.write(new_line)
        elif stop_codon_end_pos == "MISSING_STOP" and start_codon_end_pos == "MISSING_START":  # no stop AND no start
            if gene_type in coding_gene_type:
                exon_type = "CDS"
            else:
                exon_type = "NONCODING"
            new_line = line[:-1] + "\t" + exon_type + "\n"
            with open(outputPath, 'a') as o:
                o.write(new_line)
        elif start_codon_end_pos == "MISSING_START":  # no start
            # will always be a 3 UTR
            start_codon_end_pos = gene_5_end
            if strand == "+":
                ln_lst = [exon_start, exon_stop, gene_5_end, start_codon_end_pos, stop_codon_end_pos, gene_3_end]
            else:
                ln_lst = [exon_start, exon_stop, gene_5_end, stop_codon_end_pos, start_codon_end_pos, gene_3_end]
            exon_type = vec_to_str(get_type(ln_lst))
            new_line = line[:-1] + "\t" + exon_type + "\n"
            with open(outputPath, 'a') as o:
                o.write(new_line)
        elif stop_codon_end_pos == "MISSING_STOP":  # no stop
            # will always be a 5 UTR
            stop_codon_end_pos = gene_3_end
            if strand == "+":
                ln_lst = [exon_start, exon_stop, gene_5_end, start_codon_end_pos, stop_codon_end_pos, gene_3_end]
            else:
                ln_lst = [exon_start, exon_stop, gene_5_end, stop_codon_end_pos, start_codon_end_pos, gene_3_end]
            exon_type = vec_to_str(get_type(ln_lst))
            new_line = line[:-1] + "\t" + exon_type + "\n"
            with open(outputPath, 'a') as o:
                o.write(new_line)
        else:  # only 1 start and only 1 stop present
            if strand == "+":
                ln_lst = [exon_start, exon_stop, gene_5_end, start_codon_end_pos, stop_codon_end_pos, gene_3_end]
            else:
                ln_lst = [exon_start, exon_stop, gene_5_end, stop_codon_end_pos, start_codon_end_pos, gene_3_end]
            exon_type = vec_to_str(get_type(ln_lst))
            new_line = line[:-1] + "\t" + exon_type + "\n"
            with open(outputPath, 'a') as o:
                o.write(new_line)
        line_count += 1

print("Done. Output written to: " + outputPath)
notify_user("mohammed.khan2014@gmail.com", body)
