import Align


alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
# Breaking alignment runs by letter to make it easier to
# start and stop if need be.
# After runs, there should be a folder called "CSVResults"
# with preliminary alignment results consolidated for each species.
for letter in alphabet:
    Align.run(letter)

# Choose a theshold and score the alignments
out = "scored_results_60.csv"
Align.openOutFileScored(out)
alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
for letter in alphabet:
    Align.score("CSVResults\\" + letter + ".csv", out, 60)

out = "scored_results_50.csv"
Align.openOutFileScored(out)
alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
for letter in alphabet:
    Align.score("CSVResults\\" + letter + ".csv", out, 50)
