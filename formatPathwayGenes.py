from Bio import SeqIO

def run():
    gene_locations1 = {'Fe-S oxidoreductase': [0, 1245],
                      'thiolase': [1462, 2647],
                      'crotonase': [2696, 3488],
                      'beta hydroxybutyryl-CoA dehydrogenase': [3639, 4512],
                      'butyryl-CoA dehydrogenase': [4548, 5715],
                      'electron transfer flavoprotein beta-subunit': [5731, 6517],
                      'electron transfer flavoprotein alpha-subunit': [6595, 7606],
                      'putative multidrug efflux pump': [7978, 8390]}

    gene_locations2 = {'Fe-S oxidoreductase': [0, 1245],
                        ''
    for record in SeqIO.parse("butyrate_genes_1.fasta", "fasta"):
        sequence = record.seq

    gene_sequences = {} #Will contain the gene sequences broken up by their start and end positions

    #Iterate through the locations dictionary and create entries for each gene with their
    # full sequences.
    for key in gene_locations1.keys():
        gene_sequences[key] = sequence[gene_locations1[key][0]:gene_locations1[key][1]]

    outfile = open("butyrate_genes_sep1.fasta", 'a')
    for key in gene_sequences:
        outfile.write(">" + key + "\n" + str(gene_sequences[key]) + "\n")