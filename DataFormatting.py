import formatSpecies
import formatPathwayGenes

# Run the method to separate out each species and collect a list of file names
# Output - a fasta file per species in the "Species" Folder
#          a csv file containing number of fasta records per file name.
species_file = 'Gastrointestinal_tract.cds.fsa'
file_names = formatSpecies.separateSpecies(species_file)
formatSpecies.fileInfo(file_names)

# Prepare the protein seequence files
# will be deposited in a "Proteins" Folder.
formatPathwayGenes.separatePathwaySequences()


