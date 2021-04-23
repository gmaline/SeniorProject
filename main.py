# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import formatPathwayGenes
import formatSpecies
import Control
import Analysis
import csv

def run(letter):
    csv_out = letter + ".csv"
    with open("Species\\files_info.csv", 'r') as species_in:
        read_species = csv.reader(species_in, delimiter=',')

        Analysis.openOutFile(csv_out)

        for row in read_species:
            if row[0][0] == letter:
                current_species = row[0]
                db = Analysis.buildSpeciesDB(row[0])
                for protein in protein_names:
                    results = Analysis.alignSequence(current_species, protein, db)
                    Analysis.parseResults(results, current_species, protein, csv_out)
                    print(current_species + ": " + protein)

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # formatPathwayGenes.separatePathwaySequences()
    # names = formatSpecies.separateSpecies()
    # formatSpecies.file_info(names)
    # Control.run_control()

    out = "scored_results_50.csv"
    Analysis.openOutFile3(out)
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for letter in alphabet:
        Analysis.parseCSV("CSVResults\\" + letter + ".csv", out, 50)

    out = "scored_results_60.csv"
    Analysis.openOutFile3(out)
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for letter in alphabet:
        Analysis.parseCSV("CSVResults\\" + letter + ".csv", out, 60)



    protein_names = ["thiolase", "beta hydroxybutyryl-CoA dehydrogenase",
                   "crotonase", "butyryl-CoA dehydrogenase",
                    "electron transfer flavoprotein alpha-subunit",
                    "electron transfer flavoprotein beta-subunit",
                    "butyryl-COA-transferase", "phoasephate butyryltransferase",
                    "butyrate kinase"]



    # Controls: Run on 3/15/2021 at 4:37 pm
    # csv_out = "control.csv"
    # with open("Species\\control_files_info.csv", 'r') as species_in:
    #     read_species = csv.reader(species_in, delimiter=',')
    #
    #     Analysis.openOutFile(csv_out)
    #
    #     for row in read_species:
    #         current_species = row[0]
    #         db = Analysis.buildSpeciesDB(row[0])
    #         for protein in protein_names:
    #             results = Analysis.alignSequence(current_species, protein, db)
    #             Analysis.parseResults(results, current_species, protein, csv_out)
    #             print(current_species + ": " + protein)

    # Species that start with A: Run on 3/15/2021 at 4:50 pm
    # csv_out = "A.csv"
    # with open("Species\\files_info.csv", 'r') as species_in:
    #     read_species = csv.reader(species_in, delimiter=',')
    #
    #     Analysis.openOutFile(csv_out)
    #
    #     for row in read_species:
    #         if row[0][0] == 'A':
    #             current_species = row[0]
    #             db = Analysis.buildSpeciesDB(row[0])
    #             for protein in protein_names:
    #                 results = Analysis.alignSequence(current_species, protein, db)
    #                 Analysis.parseResults(results, current_species, protein, csv_out)
    #                 print(current_species + ": " + protein)

    # # Species that start with B: Run on 3/15/2021 at 5:15 pm
    # letter = 'B'
    # csv_out = letter + ".csv"
    # with open("Species\\files_info.csv", 'r') as species_in:
    #     read_species = csv.reader(species_in, delimiter=',')
    #
    #     Analysis.openOutFile(csv_out)
    #
    #     for row in read_species:
    #         if row[0][0] == letter:
    #             current_species = row[0]
    #             db = Analysis.buildSpeciesDB(row[0])
    #             for protein in protein_names:
    #                 results = Analysis.alignSequence(current_species, protein, db)
    #                 Analysis.parseResults(results, current_species, protein, csv_out)
    #                 print(current_species + ": " + protein)

    # # Species that start with C: Run on 3/15/2021 at 5:35 pm
    # letter = 'C'
    # csv_out = letter + ".csv"
    # with open("Species\\files_info.csv", 'r') as species_in:
    #     read_species = csv.reader(species_in, delimiter=',')
    #
    #     Analysis.openOutFile(csv_out)
    #
    #     for row in read_species:
    #         if row[0][0] == letter:
    #             current_species = row[0]
    #             db = Analysis.buildSpeciesDB(row[0])
    #             for protein in protein_names:
    #                 results = Analysis.alignSequence(current_species, protein, db)
    #                 Analysis.parseResults(results, current_species, protein, csv_out)
    #                 print(current_species + ": " + protein)

    # # Species that start with D: Run on 3/15/2021 at 6:39 pm
    # letter = 'D'
    # csv_out = letter + ".csv"
    # with open("Species\\files_info.csv", 'r') as species_in:
    #     read_species = csv.reader(species_in, delimiter=',')
    #
    #     Analysis.openOutFile(csv_out)
    #
    #     for row in read_species:
    #         if row[0][0] == letter:
    #             current_species = row[0]
    #             db = Analysis.buildSpeciesDB(row[0])
    #             for protein in protein_names:
    #                 results = Analysis.alignSequence(current_species, protein, db)
    #                 Analysis.parseResults(results, current_species, protein, csv_out)
    #                 print(current_species + ": " + protein)

    # Species that start with E: Run on 3/15/2021 at 6:44 pm, 7:02
    # run('E')

    # Species that start with F&G: Run on 3/15/2021 at 7:03, 7:09
    #run('F')
    #run('G')

    # Species that start with H&I: Run on 3/15/2021 at 7:09, 7:20
    # run('H')
    # run('I')

    # Species that start with J&K: Run on 3/15/2021 at 7:21, 7:22
    # run('J')
    # run('K')

    # Species that start with L&M: Run on 3/15/2021 at 7:23, 7:30
    # run('L')
    # run('M')

    # Species that start with N&O: Run on 3/15/2021 at 7:31, 7:31
    #run('N')
    #run('O')

    # Species that start with P&Q: Run on 3/15/2021 at 7:32, 7:38
    # run('P')
    # run('Q')

    # Species that start with RSTUVW: Run on 3/15/2021 at 7:39, 7:46
    # run('R')
    # run('S')
    # run('T')
    # run('U')
    # run('V')
    # run('W')

    # Species that start with XYZ: Run on 3/15/2021 at 7:46, 7:47
    # run('X')
    # run('Y')
    # run('Z')
