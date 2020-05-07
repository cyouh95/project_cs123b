from flask import Flask, url_for, render_template, redirect, request, send_file, flash
from utils import data_dir
from HBB import HBB
import settings
import zipfile
import io


app = Flask(__name__)
app.secret_key = settings.secret_key


@app.route('/')
def index():
    return render_template('index.j2')


@app.route('/results')
def results():
    return render_template('results.j2')


@app.route('/upload', methods=['POST'])
def upload():
    basic_file = request.files['basic_file']
    related_file = request.files['related_file']

    if not basic_file or not related_file:
        flash('Must upload both basic and related dataset.')
        return redirect(url_for('index'))

    basic_dataset = basic_file.stream.read().decode('utf-8')
    related_dataset = related_file.stream.read().decode('utf-8')

    hbb = HBB(basic_dataset, related_dataset)

    return redirect(url_for('results'))


@app.route('/download')
def download():
    files = ['basic_dataset_training.fasta', 'basic_dataset_testing.fasta',
             'related_dataset_training.fasta', 'related_dataset_testing.fasta']

    data = io.BytesIO()

    with zipfile.ZipFile(data, mode='w') as z:
        for file in files:
            z.write(f'{data_dir}/{file}', file)

    data.seek(0)

    return send_file(data, mimetype='application/zip', as_attachment=True, attachment_filename='output.zip')


if __name__ == '__main__':
    app.run(debug=True)