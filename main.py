# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import formatPathwayGenes
import formatSpecies
import Control
import Analysis


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    #formatPathwayGenes.separatePathwaySequences()
    #names = formatSpecies.separateSpecies()
    #formatSpecies.file_info(names)
    #Control.run_control()
    species_name = "Eubacterium_rectale_DSM_17629"
    protein_name = "thiolase"

    db = Analysis.buildSpeciesDB(species_name)
    results = Analysis.alignSequence(species_name, protein_name, db)

    Analysis.parseResults(results, species_name, protein_name, "out1.csv")


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
