from flask import Flask, render_template, request, send_from_directory, render_template_string, Response
from rdkit import Chem, rdBase
from rdkit.Chem import Draw, Descriptors, Crippen, QED, rdFreeSASA, AllChem, Lipinski
from collections import OrderedDict
from io import StringIO
import csv
from rdkit.Chem.AtomPairs import Pairs, Sheridan, Torsions, Utils

from io import BytesIO
import base64

app = Flask(__name__, template_folder='templates', static_folder='static')
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # for 16 MB max-limit

@app.route('/')
def index():
    all_descriptors = get_all_descriptors()
    return render_template('index.html', all_descriptors=all_descriptors)



@app.route('/identify_molecule', methods=['POST'])
def identify_molecule():
    global_descriptors = get_all_descriptors()  # Renamed variable here
    selected_options = request.form.getlist('displayOptions')


    # Get the checkbox state for excluding invalid SMILES
    exclude_invalid = request.form.get('excludeInvalid') == 'true'

    file = request.files.get('csvFile')
    if file and file.filename != '':
        # Process the uploaded file
        if not file.filename.endswith('.csv'):
            # Invalid file type
            return render_template('index.html', error="Invalid file type! Please upload a .csv file.")

        smiles_list = []
        csv_file = csv.reader(file.stream.read().decode("utf-8").splitlines())
        for row in csv_file:
            smiles_list.append(row[0])
    else:
        # No file uploaded, proceed with the input field
        smiles_input = request.form.get('inputField', '')
        smiles_list = [s.strip() for s in smiles_input.split(',') if s.strip() != '']

    # Exclude invalid SMILES if checkbox is checked
    if exclude_invalid:
        smiles_list = [smiles for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]

    computed_descriptors = []  # Renamed variable here
    for smiles in smiles_list:
        descriptors, img_str = compute_descriptors(smiles, selected_options)
        computed_descriptors.append(descriptors)  # Renamed variable here

    #print(selected_options)
    return render_template('index.html', descriptors_list=computed_descriptors, all_descriptors=global_descriptors)
    #return render_template('index.html', descriptors_list=all_descriptors)





def generate_csv(data):
    output = StringIO()
    writer = csv.writer(output)

    # Write headers
    writer.writerow(data.keys())

    # Write data
    writer.writerow(data.values())

    return output.getvalue()

import inspect

# ... and other relevant modules

def get_all_descriptors():
    all_descriptors = {
        'chem': {name: func for name, func in inspect.getmembers(Descriptors, inspect.isfunction) if filter_method(name)},
        'lipinski': {name: func for name, func in inspect.getmembers(Lipinski, inspect.isfunction) if filter_method(name)},
        'crippen': {name: func for name, func in inspect.getmembers(Crippen, inspect.isfunction) if filter_method(name)},
        'qed': {name: func for name, func in inspect.getmembers(QED, inspect.isfunction) if filter_method(name)},
        'rdfreesasa': {name: func for name, func in inspect.getmembers(rdFreeSASA, inspect.isfunction) if filter_method(name)},
        'allchem': {name: func for name, func in inspect.getmembers(AllChem, inspect.isfunction) if filter_method(name)},
    }

    

    return all_descriptors

### Checks if the mtehod contains name of methods that dont work.
def filter_method(name):
        not_working_methods = ["rundoctest", "auto", 'namedtuple']
        for not_working in not_working_methods:
            if not_working.lower() in name.lower():
                return False
        return True


all_descriptors = get_all_descriptors()



def compute_descriptors(smiles, selected_options):
    molecule = Chem.MolFromSmiles(smiles)
    descriptors = {}
    img_str = None  # Initialize img_str here

    if molecule is not None:
        descriptors['SMILES'] = smiles
        for option in selected_options:
            if option in all_descriptors['chem']:
                descriptors[option] = all_descriptors['chem'][option](molecule)
            elif option in all_descriptors['lipinski']:
                descriptors[option] = all_descriptors['lipinski'][option](molecule)
            elif option in all_descriptors['crippen']:
                descriptors[option] = all_descriptors['crippen'][option](molecule)
            elif option in all_descriptors['qed']:
                descriptors[option] = all_descriptors['qed'][option](molecule)
            # ... continue this pattern for other descriptor categories

            if 'Image' in selected_options:
                img = Draw.MolToImage(molecule)
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
                descriptors['Image'] = img_str
            if 'FreeSASA' in selected_options:
                # 1. Generate 3D coordinates for the molecule
                AllChem.EmbedMolecule(molecule, AllChem.ETKDG())

                # 2. Classify atoms and get radii
                radii = rdFreeSASA.classifyAtoms(molecule)

                # 3. Define SASA options: Using ShrakeRupley algorithm and OONS classifier as an example
                sasa_opts = rdFreeSASA.SASAOpts(rdFreeSASA.SASAAlgorithm.ShrakeRupley, rdFreeSASA.SASAClassifier.OONS)

                # 4. Compute the SASA and store in descriptors dictionary
                descriptors['FreeSASA'] = rdFreeSASA.CalcSASA(molecule, radii, confIdx=-1, opts=sasa_opts)
    else:
        descriptors['Error'] = f"Invalid SMILES: {smiles}"

    return descriptors, img_str



def compute_descriptors2(smiles, selected_options):
    molecule = Chem.MolFromSmiles(smiles)
    descriptors = {}
    img_str = None  # Initialize img_str here

    if molecule is not None:
        descriptors['SMILES'] = smiles
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
        if 'clogP' in selected_options:
            descriptors['clogP'] = Crippen.MolLogP(molecule)
        if 'QED' in selected_options:
            descriptors['QED'] = QED.qed(molecule)
        if 'numberOfAtoms' in selected_options:
            atom_counts = get_atom_counts(molecule)
            descriptors.update(atom_counts)
                #descriptors['numberOfAtoms'] = molecule.GetNumAtoms()
        if 'NumHDonors' in selected_options:
            descriptors['Number of H Donors'] = Lipinski.NumHDonors(molecule)
        if 'NumHAcceptors' in selected_options:
            descriptors['Number of H Acceptors'] = Lipinski.NumHAcceptors(molecule)
        if 'NumRotatableBonds' in selected_options:
            descriptors['Number of Rotatable Bonds'] = Lipinski.NumRotatableBonds(molecule)

        if 'NumAliphaticCarbocycles' in selected_options:
            descriptors['Number of Aliphatic Carbocycles'] = Lipinski.NumAliphaticCarbocycles(molecule)
        if 'NumAliphaticHeterocycles' in selected_options:
            descriptors['Number of Aliphatic Heterocycles'] = Lipinski.NumAliphaticHeterocycles(molecule)
        if 'NumAliphaticRings' in selected_options:
            descriptors['Number of Aliphatic Rings'] = Lipinski.NumAliphaticRings(molecule)
        if 'NumAromaticCarbocycles' in selected_options:
            descriptors['Number of Aromatic Carbocycles'] = Lipinski.NumAromaticCarbocycles(molecule)
        if 'NumAromaticHeterocycles' in selected_options:
            descriptors['Number of Aromatic Heterocycles'] = Lipinski.NumAromaticHeterocycles(molecule)
        if 'NumAromaticRings' in selected_options:
            descriptors['Number of Aromatic Rings'] = Lipinski.NumAromaticRings(molecule)
        if 'NumSaturatedCarbocycles' in selected_options:
            descriptors['Number of Saturated Carbocycles'] = Lipinski.NumSaturatedCarbocycles(molecule)
        if 'NumSaturatedHeterocycles' in selected_options:
            descriptors['Number of Saturated Heterocycles'] = Lipinski.NumSaturatedHeterocycles(molecule)
        if 'NumSaturatedRings' in selected_options:
            descriptors['Number of Saturated Rings'] = Lipinski.NumSaturatedRings(molecule)
        if 'FractionCSP3' in selected_options:
            descriptors['Fraction of C atoms that are SP3 hybridized'] = Lipinski.FractionCSP3(molecule)
        if 'HeavyAtomCount' in selected_options:
            descriptors['Number of Heavy Atoms'] = Lipinski.HeavyAtomCount(molecule)
        if 'NHOHCount' in selected_options:
            descriptors['Number of NHs or OHs'] = Lipinski.NHOHCount(molecule)
        if 'NOCount' in selected_options:
            descriptors['Number of Nitrogens and Oxygens'] = Lipinski.NOCount(molecule)
        if 'NumHeteroatoms' in selected_options:
            descriptors['Number of Heteroatoms'] = Lipinski.NumHeteroatoms(molecule)
        if 'RingCount' in selected_options:
            descriptors['Ring Count'] = Lipinski.RingCount(molecule)

        # rdkit.Chem.AtomPairs.Sheridan module
        if 'AssignPattyTypes' in selected_options:
            descriptors['AssignPattyTypes'] = Sheridan.AssignPattyTypes(molecule)

        # rdkit.Chem.AtomPairs.Torsions module
        if 'GetTopologicalTorsionFingerprintAsIds' in selected_options:
            descriptors['GetTopologicalTorsionFingerprintAsIds'] = Torsions.GetTopologicalTorsionFingerprintAsIds(molecule)
        if 'FreeSASA' in selected_options:
            # 1. Generate 3D coordinates for the molecule
            AllChem.EmbedMolecule(molecule, AllChem.ETKDG())

            # 2. Classify atoms and get radii
            radii = rdFreeSASA.classifyAtoms(molecule)

            # 3. Define SASA options: Using ShrakeRupley algorithm and OONS classifier as an example
            sasa_opts = rdFreeSASA.SASAOpts(rdFreeSASA.SASAAlgorithm.ShrakeRupley, rdFreeSASA.SASAClassifier.OONS)

            # 4. Compute the SASA and store in descriptors dictionary
            descriptors['FreeSASA'] = rdFreeSASA.CalcSASA(molecule, radii, confIdx=-1, opts=sasa_opts)
    else:
        descriptors['Error'] = f"Invalid SMILES: {smiles}"

    return descriptors, img_str


def get_atom_counts(molecule):
    atom_counts = OrderedDict()
    # Convert implicit hydrogens to explicit ones
    mol_with_hydrogens = Chem.AddHs(molecule)
    for atom in mol_with_hydrogens.GetAtoms():
        symbol = atom.GetSymbol()
        key = "Number of "+symbol+" atoms"
        atom_counts[key] = atom_counts.get(key, 0) + 1
    
    atom_counts["Number of atoms total:"] = sum(atom_counts.values())
    
    return atom_counts



@app.route('/feedback')
def feedback():
    return render_template('feedback.html')

@app.route('/docs')
def docs():
    return render_template('docs.html')

@app.route('/about')
def about():
    return render_template('about.html')


import os

@app.route('/editor')
def editor():
    # Get the directory of the current script
    dir_path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(dir_path, 'frontend', 'editor.html')

    with open(file_path, "r") as f:
        content = f.read()
    return render_template_string(content)

@app.route('/editor/gui/<path:filename>')
def editor_resources(filename):
    return send_from_directory('frontend/gui', filename)


@app.route('/download_csv', methods=['POST'])
def download_csv():
    
    # Fetching the selected options again
    selected_options = request.form.getlist('displayOptions')

    # Get the checkbox state for excluding invalid SMILES
    exclude_invalid = request.form.get('excludeInvalid') == 'true'
    
    smiles_input = request.form.get('inputField', '')
    smiles_list = [s.strip() for s in smiles_input.split(',')]

    # Exclude invalid SMILES if checkbox is checked
    if exclude_invalid:
        smiles_list = [smiles for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]
    
    output = StringIO()
    writer = csv.writer(output)

    # Get the RDKit version
    rdkit_version = rdBase.rdkitVersion

    # Get all descriptors for all molecules
    all_descriptors = [compute_descriptors(smiles, selected_options)[0] for smiles in smiles_list]

    # Find the complete set of descriptor keys
    all_keys = set()
    for desc in all_descriptors:
        all_keys.update(desc.keys())

    # Remove the 'Image' descriptor key if it's present and 'SMILES' key
    all_keys.discard('Image')
    all_keys.discard('SMILES')

    # Convert set to a list to maintain order
    all_keys = sorted(list(all_keys))

    # Order the atom count columns (assuming they start with "Number of")
    atom_columns = sorted([key for key in all_keys if key.startswith("Number of")])
    
    # Create the final list of columns ensuring SMILES is first and atom columns are in order
    other_keys = [key for key in all_keys if key not in atom_columns and key != 'Error']
    ordered_keys = ['SMILES'] + atom_columns + other_keys + (['Error'] if 'Error' in all_keys else [])


    # Write the headers
    writer.writerow(ordered_keys)

    # Write data for each molecule
    for desc in all_descriptors:
        writer.writerow([desc.get(key, '') for key in ordered_keys])

    csv_data = output.getvalue()

    return Response(
        csv_data,
        mimetype="text/csv",
        headers={"Content-disposition": f"attachment; filename=data_RDKit_{rdkit_version}.csv"}
    )

if __name__ == '__main__':
    app.run(debug=True)
