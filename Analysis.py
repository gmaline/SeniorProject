from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys
import os
import csv

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

    query_file = protein_name + ".fasta"
    tblastn_cmd = NcbitblastnCommandline(query=query_file, db=db_name, outfmt=5, out=out_file)
    stdout, stderr = tblastn_cmd()

    return out_file

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

    sequence = SeqIO.read(protein_name + ".fasta", "fasta")

    with open("CSVResults\\" + csv_out, 'a', newline='') as out_file:
        writer =csv.writer(out_file)
        header = ["alignment#", "hsp#", "species_name", "query_protein",
                  "hit_title", "%identity", "%query_cover", "evalue",
                  "align_length", "score", "query_start", "query_end",
                  "sbjct_start", "sbjct_end", "query", "sbjct", "identities",
                  "gaps", "positives", "bits", "out_file"]
        writer.writerow(header)
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