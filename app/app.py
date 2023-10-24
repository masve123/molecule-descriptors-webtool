import json
from flask import Flask, render_template, request, send_from_directory, render_template_string, Response
import rdkit
from rdkit import Chem, rdBase
from rdkit.Chem import Draw, Descriptors, Crippen, QED, rdFreeSASA, AllChem, Lipinski
from collections import OrderedDict
from io import StringIO, BytesIO
import csv
import base64
from flask import jsonify
import os

app = Flask(__name__, template_folder='templates', static_folder='static')
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # for 16 MB max-limit

@app.route('/')
def index():
    return render_template('index.html')


@app.route('/identify_molecule', methods=['POST'])
def identify_molecule():
    try:
        #selected_options = request.form.getlist('displayOptions')
        selected_methods_json = request.form.get('selectedMethods')
        
        if selected_methods_json is None:
            print("Error: 'selectedMethods' not received by server.")
            # Handle the error as appropriate for your application
            return "An error occurred"
        
        selected_methods = json.loads(selected_methods_json)

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
            


        selected_options = request.form.getlist('displayOptions')

        all_descriptors = []
        for smiles in smiles_list:
            descriptors = {}
            
            molecule = Chem.MolFromSmiles(smiles)
            if molecule:
                # Check if 'Image' was selected and add to the descriptors dictionary
                if 'Image' in selected_options:
                    img_str = draw_to_html(smiles, selected_options)
                    if img_str:
                        descriptors['Image'] = img_str

                for method in selected_methods:
                    module_name, method_name = method.split('.')
                    module = getattr(rdkit.Chem, module_name, None)
                    if module:
                        func = getattr(module, method_name, None)
                        if func:
                            descriptors[method_name] = func(molecule)
            else:
                if exclude_invalid:
                    continue
                descriptors['Error'] = 'Invalid Molecule'
            
            all_descriptors.append(descriptors)

        return render_template('index.html', descriptors_list=all_descriptors)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return render_template('index.html', error="An unexpected error occurred.")

def get_all_methods(*modules):
    all_methods = {}
    for module in modules:
        method_names = [attr for attr in dir(module) if callable(getattr(module, attr)) and not attr.startswith("__")]
        all_methods[module.__name__.split('.')[-1]] = method_names
    return all_methods

rdkit_methods = get_all_methods(Chem, Descriptors, Crippen, QED, rdFreeSASA, AllChem, Lipinski, Draw)


@app.route('/api/methods', methods=['GET'])
def get_methods():
    return jsonify(rdkit_methods)


def generate_csv(data):
    output = StringIO()
    writer = csv.writer(output)

    # Write headers
    writer.writerow(data.keys())

    # Write data
    writer.writerow(data.values())

    return output.getvalue()


def draw_to_html(smiles, selected_options):
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
    else:
        descriptors['Error'] = f"Invalid SMILES: {smiles}"

    return img_str

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



@app.route('/draw')
def draw():
    return render_template('draw.html')

@app.route('/docs')
def docs():
    return render_template('docs.html')

@app.route('/about')
def about():
    return render_template('about.html')


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
