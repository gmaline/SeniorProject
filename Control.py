from Bio.Blast.Applications import NcbiblastnCommandline
import sys
import os

def run_control():
    #Create the local blast database
    file_location = "Species\\"
    fasta_name = "Coprococcus_catus_GD-7"
    fasta_file = fasta_name + ".fsa"
    db_name = fasta_name + "_blastdb"
    #blastn_db = "makeblastdb -in " + file_location + fasta_file + " -dbtype nucl -input_type fasta -out " + db_name
    #os.system(blastn_db)

    #Run the blastn command on the local database
    out_file_location = "BlastResults\\"
    out_file = out_file_location + fasta_name + "_blastresults.xml"

    query_file = "butyrate_genes_sep1.fasta"
    blastn_cmd = NcbiblastnCommandline(query=query_file, db=db_name, evalue=0.001, outfmt=5, out=out_file)
    stdout, stderr = blastn_cmd()


