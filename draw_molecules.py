### Draw line representations of the compounds, store in 300 x 300 png files.

# load libraries
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import PIL
import csv

# load the list of compounds
compounds = pd.read_csv("compounds.csv")

# generate rdkit molecules from the SMILES strings
compounds["mol"] = compounds["SMILES"].apply(Chem.MolFromSmiles)

# the molecule drawer occasionally outputs an error while still producing an image
# the try/except clause stops this from terminating the script
# we will save the compounds that throw errors and write them out for later use
errors = []

# draw each molecule, save to file
for i in range(len(compounds["mol"])):
    try:
        # draw the molecule
        molecule = Draw.MolToImage(compounds["mol"][i], size = (300, 300), kekulize = True)
        
        # a filter for taking pixel intensities to 0/1
        def filter(x):
            if x == 255:
                return(255)
            else:
                return(0)

        # convert the images into black and white
        molecule = molecule.convert("L").point(filter, mode = "1")
        
        # save the image to a png file
        molecule.save("molecule_drawings/%s.png" %compounds["SubstanceID"][i])
        
    except:
        errors.append(compounds["SubstanceID"][i])
    
# write out the substance IDs of the errors
with open("errors.csv", "w") as csvfile:
    writer = csv.writer(csvfile, delimiter = ",")
    writer.writerow(errors)

