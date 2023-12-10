# Utility functions, e.g., RDKit calculations

# This file has not been utilized in the project. It is a template for future development.

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import QED
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski, AllChem
import pandas as pd
from typing import List, Dict, Union, Optional

def calculate_descriptors(smiles_string: str, selected_descriptors: Optional[List[str]] = None) -> Dict[str, Union[float, str]]:
    """
    Compute selected molecular descriptors for a given SMILES string.

    Parameters:
    - smiles_string (str): The SMILES representation of the molecule.
    - selected_descriptors (List[str], optional): List of descriptor names to compute. 
      If None, all available descriptors are computed. Defaults to None.

    Returns:
    - Dict[str, Union[float, str]]: A dictionary containing computed descriptors. The keys 
      are the descriptor names and the values are the computed values. If an error occurs, 
      a dictionary with an "error" key is returned.
    """
    mol = get_mol(smiles_string)

    # All possible descriptors with their respective calculation functions
    all_descriptors = {
        "molecular_weight": Descriptors.MolWt,
        "QED": QED.qed,
        "clogP": Crippen.MolLogP,
        "polar_surface_area": Descriptors.TPSA,
        "number_of_h-bond_donors": Lipinski.NumHDonors,
        "number_of_h-bond_acceptors": Lipinski.NumHAcceptors,
        "number_of_rotatable_bonds": Lipinski.NumRotatableBonds,
        "number_of_rings": Lipinski.RingCount,
        "number_of_aromatic_rings": Lipinski.NumAromaticRings,
        "number_of_heavy_atoms": Lipinski.HeavyAtomCount,
        "number_of_heteroatoms": Lipinski.NumHeteroatoms,
        "number_of_amide_bonds": Lipinski.NumAmideBonds,
        "number_of_rotatable_bonds": Lipinski.NumRotatableBonds,
        "number_of_strict_rotable_bonds": Lipinski.NumStrictRotatableBonds,
        "number_of_nitrogen_atoms": Lipinski.NumN,
        "number_of_carbon_atoms": Lipinski.NumC,
        "number_of_oxygen_atoms": Lipinski.NumO,
        # FreeSASA <not implemented yet>

    }

    # If no descriptors are selected, return all by default
    if not selected_descriptors:
        selected_descriptors = all_descriptors.keys()

    results = {}
    for desc in selected_descriptors:
        if desc in all_descriptors:
            results[desc] = all_descriptors[desc](mol)
    
    return results



def get_mol(smiles_string: str) -> Chem.Mol:
    """Get RDKit molecule object from SMILES string."""
    mol = Chem.MolFromSmiles(smiles_string)
    
    if mol is None:
        return {"error": "Invalid SMILES string."}
    
    #Add additional error handling here

    return mol


def calculated_descriptors_to_csv(descriptors: Dict[str, float], filename: str) -> None:
    """
    Writes the calculated molecular descriptors to a CSV file.

    Parameters:
    - descriptors (Dict[str, float]): A dictionary containing the calculated descriptors.
      The keys are the descriptor names and the values are the computed values.
    - filename (str): The name of the CSV file to which the descriptors should be written.

    Returns:
    - None: The function writes to a file and does not return a value.
    """
    # Convert the dictionary to a DataFrame
    df = pd.DataFrame([descriptors])

    # Write the DataFrame to a CSV file
    df.to_csv(filename, index=False)

# Example usage:
# descriptors = {"molecular_weight": 46.07, "QED": 0.73, "clogP": -0.64}
# calculated_descriptors_to_csv(descriptors, "output.csv")

    

if __name__ == "__main__":
    smiles = "CCO"  # Ethanol as an example
    print(calculate_descriptors(smiles))
