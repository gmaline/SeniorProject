import re
from Bio import SeqIO


# Name: separateSpecies
# Summary: parses large fasta file and separates according to species.
#          inserts separated fastas in a folder called "Species"
# Parameters: file - the name of the fasta file in this directory
# Returns: names - a list of the file names it creates.
def separateSpecies(file):
    names = []
    # Parse through the records in this large file and break them into individual files
    # for each species
    # Pulls the species name from the header using re
    for record in SeqIO.parse(file, "fasta"):
        speciesHeader = re.findall(' \[([^\[]*)\]$', str(record.description))
        # Pulls the last substring found that matches the regular expression
        thisSpecies = speciesHeader[len(speciesHeader) - 1]

        # Remove characters that would interfere with making this string a file name
        # Change slashes to dashes and spaces to underscores
        thisSpecies = re.sub('[/]', '-', thisSpecies)
        thisSpecies = re.sub('[ ]', '_', thisSpecies)

        # Preliminary step in keeping track of the file names + records per file.
        names.append(thisSpecies)

        # Write the sequence out to a file specific to the species
        outfile = open("Species\\" + thisSpecies + ".fsa", 'a')
        outfile.write(">" + str(record.description) + "\n" + str(record.seq) + "\n")
        outfile.close()
    return names

# Name: fileInfo
# Summary: creates metadata for the resulting separated fasta files and
#           how many records are in each.
# Parameters: names - a list of the names of all the files made
# Returns: NA
def fileInfo(names):
    files = {}
    #Create a record in the dictionary for each file and set to 0
    for name in names:
        files[name] = 0
    #Count the records in each file
    for name in names:
        files[name] += 1
    #Write out the file names and their count into a file called files_info.csv
    outfile = open("Species\\files_info.csv", 'a')
    for key in files.keys():
        outfile.write(key + ',' + str(files[key]) + '\n')
    outfile.close()



