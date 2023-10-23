# Import the necessary package
from rdkit.Chem import Draw, Descriptors, Crippen, QED, rdFreeSASA, AllChem, Lipinski

# Import the necessary package
from rdkit.Chem import QED

# Define your method
def list_descriptors():
    # Retrieve all callable functions from the package
    descriptor_functions = [name for name, func in vars(Descriptors).items() if callable(func)]
    for func_name in descriptor_functions:
        print(func_name)

# Usage:
list_descriptors()
