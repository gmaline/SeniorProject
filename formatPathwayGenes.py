from Bio import SeqIO
from Bio import Seq


# Name: separatePathwaySequences
# Summary: separates long CDS file into individual genes.
# Parameters: NA
# Returns: names - a list of the file names it creates.
def separatePathwaySequences():
    gene_locations_BTCnBK = {"thiolase": [1462, 2647],
                             "crotonase": [2696, 3488],
                             "beta hydroxybutyryl-CoA dehydrogenase": [3639, 4512],
                             "butyryl-CoA dehydrogenase": [4548, 5715],
                             "electron transfer flavoprotein beta-subunit": [5731, 6517],
                             "electron transfer flavoprotein alpha-subunit": [6568, 7606]}

    for record in SeqIO.parse("butyrate_genes_BTC&BK.fasta", "fasta"):
        sequence = record.seq

    protein_sequences = {}  # Will contain the resulting protein sequences.

    # Iterate through the locations dictionary and create entries for each gene with their
    # full sequences.
    for key in gene_locations_BTCnBK.keys():
        protein_sequences[key + ", BCT&BK"] = sequence[gene_locations_BTCnBK[key][0]:gene_locations_BTCnBK[key][1]]

    # Translate the separated gene sequences into proteins.
    for key in protein_sequences:
        sequence = protein_sequences[key]
        sequence = sequence.translate()
        protein_sequences[key] = sequence

    # Add the butyrate kinase protein sequence.
    for record in SeqIO.parse("butyrate_kinase_BK.fasta", "fasta"):
        sequence = record.seq
        protein_sequences["butyrate kinase, BK"] = sequence

    # Add the phosphate butyryltransferase protein sequence.
    for record in SeqIO.parse("phosphate_butyryltransferase_BK.fasta", "fasta"):
        sequence = record.seq
        protein_sequences["phoasephate butyryltransferase, BK"] = sequence

    # Add the butyryl COA: acetate CoA transferase translated gene sequence.
    for record in SeqIO.parse("butyryl_COA_transferase_BCT.fasta", "fasta"):
        sequence = record.seq
        temp = sequence[392:1733]
        sequence = temp.translate()
        protein_sequences["butyryl-COA-transferase, BCT"] = sequence

    # Write out to a file.
    outfile = open("butyrate_genes_sep1.fasta", 'w')
    for key in protein_sequences:
        protein_name = key.split(',')[0]
        outfile.write(">" + key + "\n" + str(protein_sequences[key]) + "\n")
        outfile_indvidual =  open("Proteins\\" + protein_name + ".fasta", 'w')
        outfile_indvidual.write(">" + key + "\n" + str(protein_sequences[key]) + "\n")
        outfile_indvidual.close()
    outfile.close()
