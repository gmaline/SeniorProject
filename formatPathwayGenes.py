from Bio import SeqIO


def separatePathwaySequences():
    gene_locations_BTCnBK = {"thiolase": [1462, 2647],
                             "crotonase": [2696, 3488],
                             "beta hydroxybutyryl-CoA dehydrogenase": [3639, 4512],
                             "butyryl-CoA dehydrogenase": [4548, 5715],
                             "electron transfer flavoprotein beta-subunit": [5731, 6517],
                             "electron transfer flavoprotein alpha-subunit": [6568, 7606]}

    for record in SeqIO.parse("butyrate_genes_BTC&BK.fasta", "fasta"):
        sequence = record.seq

    gene_sequences = {}  # Will contain the gene sequences broken up by their start and end positions

    # Iterate through the locations dictionary and create entries for each gene with their
    # full sequences.
    for key in gene_locations_BTCnBK.keys():
        gene_sequences[key + ", BCT&BK"] = sequence[gene_locations_BTCnBK[key][0]:gene_locations_BTCnBK[key][1]]

    for record in SeqIO.parse("butyrate_kinase_BK.fasta", "fasta"):
        sequence = record.seq
        gene_sequences["butyrate kinase, BK"] = sequence

    for record in SeqIO.parse("phosphate_butyryltransferase_BK.fasta", "fasta"):
        sequence = record.seq
        gene_sequences["phoasephate butyryltransferase, BK"] = sequence

    for record in SeqIO.parse("butyryl_COA_transferase_BCT.fasta", "fasta"):
        sequence = record.seq
        gene_sequences["butyryl-COA-transferase, BCT"] = sequence[392:1733]

    outfile = open("butyrate_genes_sep1.fasta", 'w')
    for key in gene_sequences:
        outfile.write(">" + key + "\n" + str(gene_sequences[key]) + "\n")
