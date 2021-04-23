import csv
from Bio import SeqIO
import re
import helper

# preps for ANOVA data format.
def pullTaxonomy(fasta_in):
    with open("taxon_info.csv", 'w', newline='') as outFile:
        taxon_writer = csv.writer(outFile)
        taxon_writer.writerow(["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species Formatted"])
        for record in SeqIO.parse(fasta_in, "fasta"):
            sequence = record.seq
            header = record.description
            domain = header.split(';')[0].split(" ")[1]
            phylum = header.split(';')[1]
            classx = header.split(';')[2]
            order = header.split(';')[3]
            family = header.split(';')[4]
            genus = header.split(';')[5]
            species = header.split(';')[6]
            species_formatted = re.sub('[/]', '-', species)
            species_formatted = re.sub('[ ]', '_', species_formatted)

            taxon_writer.writerow([domain, phylum, classx, order, family, genus, species, species_formatted])
    outFile.close()

def formatDataforANOVA(inFileScore, inFileTaxon):
    taxon_dict = {}
    with open(inFileTaxon, 'r') as taxonomy:
        readTaxonomy = csv.reader(taxonomy)
        for line in readTaxonomy:
            taxon_dict[line[7]] = line[2]
        print(taxon_dict)

    with open(inFileScore, 'r') as scores:
        readScores = csv.reader(scores)
        classxBCTScores = {}
        classxBKScores = {}
        line_num = 0
        classx = ''
        for line in readScores:
            try:
                classx = taxon_dict[line[0]]

            except:
                print(line[0] + " not found")

            # Create phylaScores entry for phylum if not present
            if not helper.dict_contains(classx, classxBCTScores) and len(classx) > 0:
                classxBCTScores[classx] = []
            if not helper.dict_contains(classx, classxBKScores) and len(classx) > 0:
                classxBKScores[classx] = []

            BCTscore = 0
            BKscore = 0
            if line_num != 0:
                BCTscore = int(line[2]) + int(line[3])
                BKscore = int(line[2]) + int(line[4])
            line_num += 1

            if len(classx) > 0:
                classxBCTScores[classx].append(BCTscore)
                classxBKScores[classx].append(BKscore)

        print(classxBKScores)
        print(classxBCTScores)


        with open("BK_ANOVA.csv", 'w', newline='') as Out_BK:
            writeBK = csv.writer(Out_BK)
            for c in classxBKScores.keys():
                row = []
                row.append(c)
                for x in classxBKScores[c]:
                    row.append(x)
                writeBK.writerow(row)
        Out_BK.close()
        with open("BCT_ANOVA.csv", 'w', newline='') as Out_BCT:
            writeBCT = csv.writer(Out_BCT)
            for c in classxBCTScores.keys():
                row = []
                row.append(c)
                for x in classxBCTScores[c]:
                    row.append(x)
                writeBCT.writerow(row)
        Out_BCT.close()


def formatDataforPhylogeny(inFileSilva):
    outfileSpecies = open("Silva_16s_Sequences_Species.fasta", 'w')
    outfileClass = open("Silva_16s_Sequences_Class.fasta", 'w')
    line_num = 0
    for record in SeqIO.parse(inFileSilva, "fasta"):
        header = record.description
        sequence = record.seq
        species = header.split(';')[6]
        species_formatted = re.sub('[/]', '-', species)
        species_formatted = re.sub('[ ]', '_', species_formatted)
        classx = header.split(';')[2]
        line_num += 1

        outfileSpecies.write(">" + species_formatted + "_#" + str(line_num) + "\n" + str(sequence) + "\n")
        outfileClass.write(">" + classx + "_" + str(line_num) + "\n" + str(sequence) + "\n")







#pullTaxonomy("Silva_16s_Sequences.fasta")
#formatDataforANOVA("scored_results_60.csv", "taxon_info.csv")
formatDataforPhylogeny("Silva_16s_Sequences.fasta")