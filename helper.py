
# Commonly used data for protein names
protein_names = ["thiolase", "beta hydroxybutyryl-CoA dehydrogenase",
                     "crotonase", "butyryl-CoA dehydrogenase",
                     "electron transfer flavoprotein alpha-subunit",
                     "electron transfer flavoprotein beta-subunit",
                     "butyryl-COA-transferase", "phoasephate butyryltransferase",
                     "butyrate kinase"]

protein_namesCore = ["thiolase", "beta hydroxybutyryl-CoA dehydrogenase",
                     "crotonase", "butyryl-CoA dehydrogenase",
                     "electron transfer flavoprotein alpha-subunit",
                     "electron transfer flavoprotein beta-subunit"]
protein_namesBCT = ["butyryl-COA-transferase"]
protein_namesBK = ["phoasephate butyryltransferase", "butyrate kinase"]


# Commonly used small helper methods.
def dict_contains(check_key, dict):
    for key in dict.keys():
        if key == check_key:
            return True
    return False

def has_protein(list, check_protein, threshold):
    for x in list:
        if x[0] == check_protein and x[1] > threshold:
            return True
    return False

def highest_identity(list, check_protein):
    highest = 0
    for x in list:
        if x[0] == check_protein:
            if x[1] > highest:
                highest = x[1]
    return highest

def reaches_threshold(list, threshold):
    for x in list:
        if x[1] > threshold:
            return True
    return False