from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
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


# Name: parseResults
# Summary: Parses XML results of each species - protein aligmnemnt
#           to create a csv output of significant hits.
# Parameters: xmlFile - the blast results file to consider
#             protein_name - the name of the pathway protein query.
#             species_name - the name of the species subject.
#             csv_out - the csv file to write results out to.
# Returns: NA
def parseResults(xmlFile, species_name, protein_name, csv_out):
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



# Name: score
# Summary: Parses CSV results to score each bacteria.
# Parameters: csv_in - the csv file to score species from
# Returns: NA
def score(csv_in, csv_out, threshold):
    protein_namesCore = helper.protein_namesCore
    protein_namesBCT = helper.protein_namesBCT
    protein_namesBK = helper.protein_namesBK

    # Create a dict entry for each protein in our list.
    temp_dict_prot = {}
    for protein in protein_namesBK:
        temp_dict_prot[protein] = ""
    for protein in protein_namesBCT:
        temp_dict_prot[protein] = ""
    for protein in protein_namesCore:
        temp_dict_prot[protein] = ""

    temp_dict = {}    # will hold a list of scores for species.
    # Read in the csv results file
    with open(csv_in, 'r') as in_file:
        scorer = csv.reader(in_file, delimiter=',')
        row_num = 0
        # For each row that is not the header, score results.
        for row in scorer:
            if row_num > 0:
                species_name = row[2]
                protein_name = row[3]
                protein_identity = float(row[5])

                # create an entry for that
                if not helper.dict_contains(species_name, temp_dict):
                    temp_dict[species_name] = []

                # Save results in a temporary list and create new entry
                temp_list = temp_dict[species_name]
                temp_list.append([protein_name, protein_identity])
                temp_dict[species_name] = temp_list
            row_num += 1
    in_file.close()

    # After results are saved as a species dictionary list of lists of proteins
    # write the results out to a file.
    with open(csv_out, 'a', newline='') as outfile:
        write_out = csv.writer(outfile)

        for key in temp_dict.keys():
            # If there is a protein that reaches the threshold
            if helper.reaches_threshold(temp_dict[key], threshold):
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
                    if helper.has_protein(temp_dict[key], protein, threshold):
                        total_scoreCore += 1
                    protein_dict[protein] = helper.highest_identity(temp_dict[key], protein)

                for protein in protein_namesBCT:
                    if helper.has_protein(temp_dict[key], protein, threshold):
                        total_scoreBCT += 1
                    protein_dict[protein] = helper.highest_identity(temp_dict[key], protein)

                for protein in protein_namesBK:
                    if helper.has_protein(temp_dict[key], protein, threshold):
                        total_scoreBK += 1
                    protein_dict[protein] = helper.highest_identity(temp_dict[key], protein)

                # Create labels for capability.
                if (total_scoreBCT + total_scoreCore == 7):
                    definite_BCT = "yes"
                if (total_scoreBK + total_scoreCore == 8):
                    definite_BK = "yes"

                # Format the row to write to the csv.
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


# Name: run
# Summary: Performs alignments and outputs to "Output.csv"
# Parameters: letter - the species name starting letter for this round of runs.
# Returns: NA
def run(letter):
    csv_out = letter + ".csv"
    with open("Species\\files_info.csv", 'r') as species_in:
        read_species = csv.reader(species_in, delimiter=',')

        # Writes the header for the appendable output file.
        openOutFile(csv_out)
        # For all of the separate file names, create a blast db
        for row in read_species:
            if row[0][0] == letter:
                current_species = row[0]
                db = buildSpeciesDB(current_species)
                # For all of the protein files, run an alignment against the blastdb
                # parse the results and create a csv output file.
                for protein in helper.protein_names:
                    results = alignSequence(current_species, protein, db)
                    parseResults(results, current_species, protein, csv_out)
                    # Track progress in the console
                    print(current_species + ": " + protein)


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

# Name: openOutFileScored
# Summary: Creates the output file and writes its header.
# Parameters: csv_out - the name of the output file
# Returns: NA
def openOutFileScored(csv_out):
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

