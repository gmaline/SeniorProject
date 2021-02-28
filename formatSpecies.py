import re
from Bio import SeqIO
from pathlib import Path

def separateSpecies():
    file = Path(r'C:\Users\User\Desktop\Senior_Project\Gastrointestinal_tract.cds.fsa')
    check = '>CW1_2704 conserved hypothetical protein [else] [Bacteroides [easy xylanisolvens/ SD CC 2a]'
    search = re.findall(' \[(.*?)\]', check)
    result = search[len(search) - 1]
    result = re.sub('[/]', '-', result)
    print(result)
    names = []
    # Parse through the records in this large file and break them into individual files
    # for each species
    for record in SeqIO.parse(file, "fasta"):
       speciesHeader = re.findall(' \[(.*?)\]', str(record.description))
       thisSpecies = speciesHeader[len(speciesHeader) - 1]

       # Remove characters that would interfere with making this string a file name
       thisSpecies = re.sub('[/]', '-', thisSpecies)

       # Preliminary step in keeping track of the file names + records per file.
       names.append(thisSpecies)

        # Write the sequence out to a file specific to the species
       #outfile = open("Species\\" + thisSpecies + ".fsa", 'a')
       #outfile.write(">" + str(record.description) + "\n" + str(record.seq) + "\n")
       #outfile.close()
    return names

def file_info(names):
    files = {}
    for name in names:
        files[name] = 0
    for name in names:
        files[name] += 1
    outfile = open("Species\\files_info.csv", 'a')
    for key in files.keys():
        outfile.write(key + ',' + str(files[key]) + '\n')
    outfile.close()



