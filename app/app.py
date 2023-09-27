from flask import Flask, render_template, request
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

from io import BytesIO
import base64

#from app.routes import main  # Replace 'your_folder_name' with the actual folder name where routes.py is located

app = Flask(__name__)
#app.register_blueprint(main)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/identify_molecule', methods=['POST'])
def identify_molecule():
    smiles = request.form.get('inputField', '')
    selected_options = request.form.getlist('displayOptions')
    print(type(smiles))
    print(smiles)
    molecule = Chem.MolFromSmiles(smiles)
    
    descriptors = {}
    img_str = None  # Initialize img_str here
    if molecule is not None:
        if 'Image' in selected_options:
            img = Draw.MolToImage(molecule)
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
            descriptors['Image'] = img_str
        if 'MolecularWeight' in selected_options:
            descriptors['MolecularWeight'] = Descriptors.ExactMolWt(molecule)
        if 'PSA' in selected_options:
            descriptors['PSA'] = Descriptors.TPSA(molecule)
        # Add other descriptors here based on selected_options
    else:
        descriptors['Error'] = "Invalid SMILES"
    
    return render_template('index.html', descriptors=descriptors, image=img_str)

@app.route('/draw')
def draw():
    return render_template('draw.html')

@app.route('/docs')
def docs():
    return render_template('docs.html')

@app.route('/about')
def about():
    return render_template('about.html')

if __name__ == '__main__':
    app.run(debug=True)
