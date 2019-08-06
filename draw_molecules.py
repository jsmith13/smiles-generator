### Draw line representations of the compounds, store in 100 x 100 png files.

# load libraries
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd

# load the list of compounds
compounds = pd.read_csv("compounds.csv")

# generate rdkit molecules from the SMILES strings
compounds["mol"] = compounds["SMILES"].apply(Chem.MolFromSmiles)

# remove the structures that rdkit cannot understand
compounds = compounds[compounds["mol"].isnull() == False]

# draw each molecule, write to file
# the molecule drawer occasionally outputs an error while still producing an image, the try/except clause stops this from terminating the script
for i in range(len(compounds["mol"])):
    try:
        Draw.MolToFile(compounds["mol"][i], fileName = "molecule_drawings/%s.png" %compounds["SubstanceID"][i], size = (300, 300))
    except:
        print(i)
    