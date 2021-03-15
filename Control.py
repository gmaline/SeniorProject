from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys
import os

def run_control():
    #Create the local blast database
    file_location = "Species\\"
    # Change THIS LINE to change the species for analysis and file output
    fasta_name = "Eubacterium_rectale_DSM_17629"
    fasta_file = fasta_name + ".fsa"
    db_name = fasta_name + "_blastdb"
    blastn_db = "makeblastdb -in " + file_location + fasta_file + " -dbtype nucl -input_type fasta -out " + db_name
    os.system(blastn_db)

    #Run the blastn command on the local database
    out_file_location = "BlastResults\\"
    out_file = out_file_location + fasta_name + "_blastresults.xml"

    query_file = "thiolase.fasta"
    tblastn_cmd = NcbitblastnCommandline(query=query_file, db=db_name, outfmt=5, out=out_file)
    stdout, stderr = tblastn_cmd()

    result_handle = open(out_file)
    blast_records = NCBIXML.parse(result_handle)

    sequence = SeqIO.read(query_file, "fasta")
    print("query length: " + str(len(sequence.seq)))

    hit_len = 0
    for record in blast_records:
        for alignment in record.alignments:
                for hsp in alignment.hsps:
                    list = []

                    print("***Alignment***")
                    print("title: " + str(alignment.title))
                    print("hit id: " + str(alignment.hit_id))
                    print("length: " + str(alignment.length))
                    # print("evalue: " + str(hsp.expect))
                    # print("score: " + str(hsp.score))
                    # print("bits: " + str(hsp.bits))
                    print("num_alignments: " + str(hsp.num_alignments))
                    print("identities: " + str(hsp.identities))
                    # print("positives: " + str(hsp.positives))
                    # print("gaps: " + str(hsp.gaps))
                    print("align_length: " + str(hsp.align_length))
                    print("stand: " + str(hsp.strand))
                    print("query: " + str(hsp.query))
                    print("query_start: " + str(hsp.query_start))
                    print("query_end: " + str(hsp.query_end))
                    # print("sbjct: " + str(hsp.sbjct))
                    # print("sbjct_start: " + str(hsp.sbjct_start))
                    # print("sbjct_end: " + str(hsp.sbjct_end))
                    print("Query Cover: " + str(hsp.align_length/len(sequence.seq) * 100))
                    print("Identity: " + str(hsp.identities/len(sequence.seq) * 100))
