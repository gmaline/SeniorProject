from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys
import os
import csv
import helper

# Name: buildSpeciesDB
# Summary: This builds a local nucleotide blastdb for the input species.
# Parameters: species_name - describes the name of the current species for analysis
# Returns: db_name - the name of the db that can be queried
def buildSpeciesDB(species_name):
    # Species the file location of the fasta file to be made into db
    file_location = "Species\\"
    fasta_file = species_name + ".fsa"
    # Writes over the current blast db each time to save space
    db_name = "temp_blastdb"

    # Creates the db with the right input file and the temp output db names
    blastn_db = "makeblastdb -in " + file_location + fasta_file + " -dbtype nucl -input_type fasta -out " + db_name
    os.system(blastn_db)
    return db_name

# Name: alignSequence
# Summary: This aligns a specied query protein against our newly made blast db.
# Parameters: filename - describes the name of the current file to query, will be one of the
#                         multiple pathway proteins.
#             db_name - the name of the db to query against
# Returns: out_file - the name of the blast result xml file
def alignSequence(species_name, protein_name, db_name):
    # Run the blastn command on the local database
    out_file_location = "BlastResults\\"
    out_file = out_file_location + species_name + protein_name + "_blastresults.xml"

    query_file = "Proteins\\" + protein_name + ".fasta"
    tblastn_cmd = NcbitblastnCommandline(query=query_file, db=db_name, outfmt=5, out=out_file)
    stdout, stderr = tblastn_cmd()

    return out_file

# Name: openOutFile
# Summary: Creates the output file and writes its header.
# Parameters: csv_out - the name of the output file
# Returns: NA
def openOutFile(csv_out):
    with open("CSVResults\\" + csv_out, 'w', newline='') as out_file:
        writer =csv.writer(out_file)
        header = ["alignment#", "hsp#", "species_name", "query_protein",
                  "hit_title", "%identity", "%query_cover", "evalue",
                  "align_length", "score", "query_start", "query_end",
                  "sbjct_start", "sbjct_end", "query", "sbjct", "identities",
                  "gaps", "positives", "bits", "out_file"]
        writer.writerow(header)

# Name: parseResults
# Summary: Parses XML results to create a csv output of significant hits.
# Parameters: xmlFile - the blast results file to consider
#             protein_name - the name of the pathway protein query.
#             species_name - the name of the species subject.
#             csv_out - the csv file to write results out to.
# Returns: NA
def parseResults(xmlFile, species_name, protein_name,  csv_out):
    result_handle = open(xmlFile)
    blast_records = NCBIXML.parse(result_handle)

    sequence = SeqIO.read("Proteins\\" + protein_name + ".fasta", "fasta")

    with open("CSVResults\\" + csv_out, 'a', newline='') as out_file:
        writer =csv.writer(out_file)

        for record in blast_records:
            alignment_num = 0
            for alignment in record.alignments:
                alignment_num += 1
                hsp_num = 0
                for hsp in alignment.hsps:
                    hsp_num += 1
                    if hsp.expect < .01:
                        list = []
                        list.append(alignment_num)
                        list.append(hsp_num)
                        list.append(species_name)
                        list.append(protein_name)
                        list.append(alignment.title)
                        list.append(str(hsp.identities / len(sequence.seq) * 100))
                        list.append(str(hsp.align_length / len(sequence.seq) * 100))
                        list.append(hsp.expect)
                        list.append(hsp.align_length)
                        list.append(hsp.score)
                        list.append(hsp.query_start)
                        list.append(hsp.query_end)
                        list.append(hsp.sbjct_start)
                        list.append(hsp.sbjct_end)
                        list.append(hsp.query)
                        list.append(hsp.sbjct)
                        list.append(hsp.identities)
                        list.append(hsp.gaps)
                        list.append(hsp.positives)
                        list.append(hsp.bits)
                        list.append(xmlFile)

                        writer.writerow(list)
    out_file.close()

def openOutFile2():
    with open("scored_results.csv", 'w', newline='') as outfile:
        write_out = csv.writer(outfile)

        header = ["species", "core", "BCT", "BK", "thiolase",
                  "beta hydroxybutyryl-CoA dehydrogenase",
                  "crotonase", "butyryl-CoA dehydrogenase",
                  "electron transfer flavoprotein alpha-subunit",
                  "electron transfer flavoprotein beta-subunit",
                  "butyryl-COA-transferase", "phoasephate butyryltransferase",
                  "butyrate kinase", "definite_BCT", "definite_BK"]
        write_out.writerow(header)

def openOutFile3(csv_out):
    with open(csv_out, 'w', newline='') as outfile:
        write_out = csv.writer(outfile)

        header = ["species", "core", "BCT", "BK", "thiolase",
                  "beta hydroxybutyryl-CoA dehydrogenase",
                  "crotonase", "butyryl-CoA dehydrogenase",
                  "electron transfer flavoprotein alpha-subunit",
                  "electron transfer flavoprotein beta-subunit",
                  "butyryl-COA-transferase", "phoasephate butyryltransferase",
                  "butyrate kinase", "definite_BCT", "definite_BK"]
        write_out.writerow(header)

# Name: parseCSV
# Summary: Parses CSV results to score each bacteria.
# Parameters: csv_in - the csv file to score species from
# Returns: NA
def parseCSV(csv_in, csv_out, threshold):
    protein_namesCore = ["thiolase", "beta hydroxybutyryl-CoA dehydrogenase",
                         "crotonase", "butyryl-CoA dehydrogenase",
                         "electron transfer flavoprotein alpha-subunit",
                         "electron transfer flavoprotein beta-subunit"]
    protein_namesBCT = ["butyryl-COA-transferase"]
    protein_namesBK = ["phoasephate butyryltransferase", "butyrate kinase"]
    # the cutoff for identity that will be considered for a protein score

    temp_dict = {}

    current_species = ""
    temp_dict_prot = {}
    template_temp_dict_prot = {}
    for protein in protein_namesBK:
        temp_dict_prot[protein] = ""
        template_temp_dict_prot[protein] = ""
    for protein in protein_namesBCT:
        temp_dict_prot[protein] = ""
        template_temp_dict_prot[protein] = ""
    for protein in protein_namesCore:
        temp_dict_prot[protein] = ""
        template_temp_dict_prot[protein] = ""


    with open(csv_in, 'r') as in_file:
        scorer = csv.reader(in_file, delimiter=',')
        row_num = 0
        for row in scorer:
            if row_num > 0:
                species_name = row[2]

                if not helper.dict_contains(species_name, temp_dict):
                    temp_dict[species_name] = []

                # Save results in a temporary list and create new entry
                temp_list = temp_dict[species_name]
                temp_list.append([row[3], float(row[5])])
                temp_dict[species_name] = temp_list
            row_num += 1
    in_file.close()

    with open(csv_out, 'a', newline='') as outfile:
        write_out = csv.writer(outfile)

        for key in temp_dict.keys():
            # If there is a protein that reaches the threshold
            if reaches_threshold(temp_dict[key], threshold):
                total_scoreCore = 0
                total_scoreBCT = 0
                total_scoreBK = 0
                protein_dict = {
                    "thiolase": 0,
                    "beta hydroxybutyryl-CoA dehydrogenase": 0,
                    "crotonase": 0,
                    "butyryl-CoA dehydrogenase": 0,
                    "electron transfer flavoprotein alpha-subunit": 0,
                    "electron transfer flavoprotein beta-subunit": 0,
                    "butyryl-COA-transferase": 0,
                    "phoasephate butyryltransferase": 0,
                    "butyrate kinase": 0
                }
                definite_BCT = "no"
                definite_BK = "no"

                # Score the available proteins by pathway
                for protein in protein_namesCore:
                    if has_protein(temp_dict[key], protein, threshold):
                        total_scoreCore += 1
                    protein_dict[protein] = highest_identity(temp_dict[key], protein)

                for protein in protein_namesBCT:
                    if has_protein(temp_dict[key], protein, threshold):
                        total_scoreBCT += 1
                    protein_dict[protein] = highest_identity(temp_dict[key], protein)

                for protein in protein_namesBK:
                    if has_protein(temp_dict[key], protein, threshold):
                        total_scoreBK += 1
                    protein_dict[protein] = highest_identity(temp_dict[key], protein)

                if (total_scoreBCT + total_scoreCore == 7):
                    definite_BCT = "yes"
                if (total_scoreBK + total_scoreCore == 8):
                    definite_BK = "yes"
                row = [key, total_scoreCore, total_scoreBCT, total_scoreBK,
                       protein_dict["thiolase"],
                       protein_dict["beta hydroxybutyryl-CoA dehydrogenase"],
                       protein_dict["crotonase"],
                       protein_dict["butyryl-CoA dehydrogenase"],
                       protein_dict["electron transfer flavoprotein alpha-subunit"],
                       protein_dict["electron transfer flavoprotein beta-subunit"],
                       protein_dict["butyryl-COA-transferase"],
                       protein_dict["phoasephate butyryltransferase"],
                       protein_dict["butyrate kinase"],
                       definite_BCT, definite_BK]
                write_out.writerow(row)
    outfile.close()











# Name: parseCSVOldVersion
# Summary: Parses CSV results to score each bacteria. This version does
#           not keep track of any results that do not reach threshold.
# Parameters: csv_in - the csv file to score species from
# Returns: NA
def parseCSVOldVersion(csv_in):
    protein_namesCore = ["thiolase", "beta hydroxybutyryl-CoA dehydrogenase",
                     "crotonase", "butyryl-CoA dehydrogenase",
                     "electron transfer flavoprotein alpha-subunit",
                     "electron transfer flavoprotein beta-subunit"]
    protein_namesBCT = ["butyryl-COA-transferase"]
    protein_namesBK = ["phoasephate butyryltransferase", "butyrate kinase"]
    # the cutoff for identity that will be considered for a protein score
    threshold = 60
    temp_dict = {}


    with open(csv_in, 'r') as in_file:
        scorer = csv.reader(in_file, delimiter=',')
        row_num = 0
        for row in scorer:
            atThreshold = False
            if row_num > 0:
                species_name = row[2]
                if float(row[5]) >= threshold:
                    atThreshold = True

                # Check if it matches the threshold and if it does add it
                if float(row[5]) >= threshold:
                    # Create an empty list if there is not already an entry
                    if not helper.dict_contains(species_name, temp_dict):
                        temp_dict[species_name] = []

                    # Save results in a temporary list and create new entry
                    temp_list = temp_dict[species_name]
                    temp_list.append([row[3], float(row[5])])
                    temp_dict[species_name] = temp_list

            row_num += 1
    in_file.close()

    with open("scored_results.csv", 'a', newline='') as outfile:
        write_out = csv.writer(outfile)
        for key in temp_dict:
            total_scoreCore = 0
            total_scoreBCT = 0
            total_scoreBK = 0
            protein_dict = {
                "thiolase": 0,
                "beta hydroxybutyryl-CoA dehydrogenase": 0,
                "crotonase": 0,
                "butyryl-CoA dehydrogenase": 0,
                "electron transfer flavoprotein alpha-subunit": 0,
                "electron transfer flavoprotein beta-subunit": 0,
                "butyryl-COA-transferase": 0,
                "phoasephate butyryltransferase": 0,
                "butyrate kinase": 0
            }
            definite_BCT = "no"
            definite_BK = "no"

            for protein in protein_namesCore:
                if has_protein(temp_dict[key], protein, threshold):
                    total_scoreCore += 1
                    protein_dict[protein] = highest_identity(temp_dict[key], protein)


            for protein in protein_namesBCT:
                if has_protein(temp_dict[key], protein, threshold):
                    total_scoreBCT += 1
                    protein_dict[protein] = highest_identity(temp_dict[key], protein)


            for protein in protein_namesBK:
                if has_protein(temp_dict[key], protein, threshold):
                    total_scoreBK += 1
                    protein_dict[protein] = highest_identity(temp_dict[key], protein)

            if (total_scoreBCT + total_scoreCore == 7):
                definite_BCT = "yes"
            if (total_scoreBK + total_scoreCore == 8):
                definite_BK = "yes"
            row = [key, total_scoreCore, total_scoreBCT, total_scoreBK,
                    protein_dict["thiolase"],
                    protein_dict["beta hydroxybutyryl-CoA dehydrogenase"],
                    protein_dict["crotonase"],
                    protein_dict["butyryl-CoA dehydrogenase"],
                    protein_dict["electron transfer flavoprotein alpha-subunit"],
                    protein_dict["electron transfer flavoprotein beta-subunit"],
                    protein_dict["butyryl-COA-transferase"],
                    protein_dict["phoasephate butyryltransferase"],
                    protein_dict["butyrate kinase"],
                   definite_BCT, definite_BK]
            write_out.writerow(row)
    outfile.close()





def has_protein(list, check_protein, threshold):
    for x in list:
        if x[0] == check_protein and x[1] > threshold:
            return True
    return False



def highest_identity(list, check_protein):
    highest = 0
    for x in list:
        if x[0] == check_protein:
            if x[1] > highest:
                highest = x[1]
    return highest

def reaches_threshold(list, threshold):
    for x in list:
        if x[1] > threshold:
            return True
    return False