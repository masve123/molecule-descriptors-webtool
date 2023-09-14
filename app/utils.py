# Utility functions, e.g., RDKit calculations

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import QED
from rdkit.Chem import Crippen
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
        # You can add more descriptors here as needed
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
