from __future__ import division
exons = []
start_codons = {}
stop_codons = {}
genes = {}
filePath = "C:/Users/mjk140030/Desktop/GTEx/toy_gtf.gtf"
outputFilePath = "C:/Users/mjk140030/Desktop/GTEx/toy_gtf2.txt"
with open(filePath) as f:  # setup lists
    print("Opened file. Setting up dictionary definitions..."),
    for line in f:
        if "\texon\t" in line:
            exons.append(line)
        elif "\tstart_codon\t" in line:
            elements = line.replace('\"', '').replace(';', '').split()
            key = elements[11]  # transcriptID
            if key in start_codons:  # if transcriptID already in dictionary,
                value = start_codons[key] + ',' + elements[4]
            else:
                value = elements[4]  # end
            start_codons[key] = value
        elif "\tstop_codon\t" in line:
            elements = line.replace('\"', '').replace(';', '').split()
            key = elements[11]  # transcriptID
            if key in stop_codons:  # if transcriptID already in dictionary,
                value = stop_codons[key] + ',' + elements[4]
            else:
                value = elements[4]  # end
            stop_codons[key] = value
        elif "\tgene\t" in line:
            elements = line.replace('\"', '').replace(';', '').split()
            key = elements[9]  # gene_id
            value = elements[3] + "," + elements[4]  # gene_start,gene_stop
            genes[key] = value
print("done.")
# example elements array
# ['0chr1', '1HAVANA', '2exon', '11869', '12227', '5.', '6+', '7.', '8gene_id', '9ENSG00000223972.4', '10transcript_id',
# '11ENST00000456328.2', '12gene_type', '13pseudogene', '14gene_status', '15KNOWN', '16gene_name', '17DDX11L1',
# '18transcript_type', '19processed_transcript', '20transcript_status', '21KNOWN', '22transcript_name', '23DDX11L1-002',
#  '24exon_number', '1', '26exon_id', '27ENSE00002234944.1', '28level', '2', '30tag', 'basic', 'havana_gene',
# '33OTTHUMG00000000961.2', 'havana_transcript', '35OTTHUMT00000362751.1']


def getStart(transcriptID):  # will return 2 start positions with comma as delimiter if two starts per transcript
    try:
        return start_codons[transcriptID]
    except Exception:
        return "MISSING_START"


def getStop(transcriptID):  # will return 2 stop positions with comma as delimiter if two stops per transcript
    try:
        return stop_codons[transcriptID]
    except Exception:
        return "MISSING_STOP"

print("Creating file..."),
length = len(exons)
for x in range(0, length):
    if x % 1200 == 0:
        val = (x/length) * 100
        print(str(round(val, 2))+"%")
    elem = exons[x].replace('\"', '').replace(';', '').split()
    start_codon = getStart(elem[11])
    stop_codon = getStop(elem[11])
    gene_start = genes[elem[9]].split(",")[0]
    gene_stop = genes[elem[9]].split(",")[1]
    line = elem[0] + "_" + elem[3] + "_" + elem[4] + "_" + elem[6] + "\t" + elem[9] + "\t" + elem[11] + "\t" + elem[
        13] + "\t" + elem[17] + "\t" + elem[
               25] + "\t" + start_codon + "\t" + stop_codon + "\t" + gene_start + "\t" + gene_stop + "\n"
    with open(outputFilePath, "a") as f:
      f.write(line)
print("Done. Output written to file: "+outputFilePath)
