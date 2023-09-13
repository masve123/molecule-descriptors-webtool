# Utility functions, e.g., RDKit calculations

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import QED
from rdkit.Chem import Crippen
import csv


def calculate_descriptors(smiles_string, selected_descriptors=None) -> dict:
    """Compute molecular descriptors for a given SMILES string.

    Parameters:
    - smiles_string (str): SMILES representation of the molecule.
    - selected_descriptors (list): List of descriptors to calculate. If None, all descriptors are calculated.

    Returns:
    - dict: Dictionary containing molecular descriptors."""
    
    mol = get_mol(smiles_string)

    # All possible descriptors with their respective calculation functions
    all_descriptors = {
        "molecular_weight": Descriptors.MolWt,
        "QED": QED.qed,
        "clogP": Crippen.MolLogP,
        # You can add more descriptors here
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


def calculated_descriptors_to_csv(descriptors: dict) -> str:
    with open('descriptors.csv', 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=descriptors.keys())
        writer.writeheader()
        writer.writerow(descriptors)
    

if __name__ == "__main__":
    smiles = "CCO"  # Ethanol as an example
    print(calculate_descriptors(smiles))
