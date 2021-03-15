# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import formatPathwayGenes
import formatSpecies
import Control
import Analysis
import csv


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # formatPathwayGenes.separatePathwaySequences()
    # names = formatSpecies.separateSpecies()
    # formatSpecies.file_info(names)
    # Control.run_control()

    protein_names = ["thiolase", "beta hydroxybutyryl-CoA dehydrogenase",
                   "crotonase", "butyryl-CoA dehydrogenase",
                    "electron transfer flavoprotein alpha-subunit",
                    "electron transfer flavoprotein beta-subunit",
                    "butyryl-COA-transferase", "phoasephate butyryltransferase",
                    "butyrate kinase"]

    # Controls: Run on 3/15/2021 at 4:37 pm
    csv_out = "control.csv"
    with open("Species\\control_files_info.csv", 'r') as species_in:
        read_species = csv.reader(species_in, delimiter=',')

        Analysis.openOutFile(csv_out)

        for row in read_species:
            current_species = row[0]
            db = Analysis.buildSpeciesDB(row[0])
            for protein in protein_names:
                results = Analysis.alignSequence(current_species, protein, db)
                Analysis.parseResults(results, current_species, protein, csv_out)
                print(current_species + ": " + protein)





