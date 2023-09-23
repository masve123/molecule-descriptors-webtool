from flask import Flask, render_template
from app.routes import main  # Replace 'your_folder_name' with the actual folder name where routes.py is located

app = Flask(__name__)
app.register_blueprint(main)

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)
