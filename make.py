# Create the complete HTML and Flask Python files with the requested modifications

html_content = """
<!DOCTYPE html>
<html lang="en">

<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <link rel="stylesheet" type="text/css" href="static/css/styles.css">
    <title>MolDescriptor</title>
</head>

<body>

    <div class="top-bar">
        <div class="logo"><a class="title-link" href="index.html">MolDescriptor</a></div>
        <div class="links">
            <a href="pages/draw.html" class="draw-link">DRAW</a>
            <a href="pages/docs.html" class="docs-link">DOCS</a>
            <a href="pages/about.html" class="about-link">OM</a>
        </div>
    </div>

    <div class="form-container">
        <h1 class="subheader">PyDescriptor</h1>
        <form action="/identify_molecule" method="post">
            <label for="inputField">Input a molecule using the 'SMILES' format:</label>
            <br>
            <input type="text" id="smiles" name="inputField" placeholder="n1ccccc1" required>
            <div class="checkbox-container">
                <label><input type="checkbox" name="displayOptions" value="MolecularWeight"> Molecular Weight</label>
                <label><input type="checkbox" name="displayOptions" value="PSA"> Polar Surface Area</label>
                <!-- Add other checkboxes here -->
            </div>
            <input type="submit" value="Get Result">
        </form>
        {% if descriptors %}
            <h2>Descriptors:</h2>
            <ul>
            {% for key, value in descriptors.items() %}
                <li>{{ key }}: {{ value }}</li>
            {% endfor %}
            </ul>
        {% endif %}
        <div id="resultArea"></div>
    </div>
</body>

</html>
"""

python_content = """
from flask import Flask, render_template, request
from rdkit import Chem
from rdkit.Chem import Descriptors

app = Flask(__name__)

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@app.route('/identify_molecule', methods=['POST'])
def identify_molecule():
    smiles = request.form.get('inputField', '')
    selected_options = request.form.getlist('displayOptions')
    
    molecule = Chem.MolFromSmiles(smiles)
    
    descriptors = {}
    if molecule is not None:
        if 'MolecularWeight' in selected_options:
            descriptors['MolecularWeight'] = Descriptors.ExactMolWt(molecule)
        if 'PSA' in selected_options:
            descriptors['PSA'] = Descriptors.TPSA(molecule)
        # Add other descriptors here based on selected_options
    else:
        descriptors['Error'] = "Invalid SMILES"
    
    return render_template('index.html', descriptors=descriptors)

if __name__ == '__main__':
    app.run(debug=True)
"""
import os
# Update the existing files in the project directory
project_folder = '/mnt/data/updated_project_with_checkbox'
os.makedirs(f'{project_folder}/app/templates', exist_ok=True)
os.makedirs(f'{project_folder}/app/static/css', exist_ok=True)

# Save the new HTML and Python files
with open(f'{project_folder}/app/templates/index.html', 'w') as f:
    f.write(html_content)

with open(f'{project_folder}/app/app.py', 'w') as f:
    f.write(python_content)

# Create a new ZIP file with the updated structure
updated_zip_with_checkbox_path = '/mnt/data/updated_project_with_checkbox.zip'
with ZipFile(updated_zip_with_checkbox_path, 'w') as zipf:
    for root, _, files in os.walk(project_folder):
        for file in files:
            zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), project_folder))

updated_zip_with_checkbox_path
