from flask import Flask, render_template, request
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
app = Flask(__name__)
import utils # legg til funksjoner her og endre uthenting av molekyler

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@app.route('/identify_molecule', methods=['POST'])
def identify_molecule():
    smiles = request.form.get('inputField', '')
    
    molecule = utils.get_mol(smiles)
    
    if molecule is None:
        return render_template('index.html', molecule_name="Invalid SMILES")
    
    molecule_weight = utils.calculate_descriptors(molecule)
    return render_template('index.html', molecule_name=molecule_weight)

if __name__ == '__main__':
    app.run(debug=True)
