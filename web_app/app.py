import os
import uuid
import shutil
import zipfile
import markdown
from pathlib import Path
from flask import Flask, render_template, request, redirect, url_for, send_from_directory, flash

# Import local analysis functions
import protein_analysis

app = Flask(__name__)
app.secret_key = "protein_secret_key"

# Configuration
UPLOAD_FOLDER = Path(__file__).parent / 'static' / 'uploads'
RESULTS_FOLDER = Path(__file__).parent / 'static' / 'results'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['RESULTS_FOLDER'] = RESULTS_FOLDER

# Ensure directories exist
UPLOAD_FOLDER.mkdir(parents=True, exist_ok=True)
RESULTS_FOLDER.mkdir(parents=True, exist_ok=True)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])
def analyze():
    if 'pdb_file' not in request.files:
        flash('No file part')
        return redirect(request.url)
    
    file = request.files['pdb_file']
    if file.filename == '':
        flash('No selected file')
        return redirect(request.url)
    
    if file and file.filename.endswith('.pdb'):
        job_id = str(uuid.uuid4())
        job_dir = RESULTS_FOLDER / job_id
        job_dir.mkdir(parents=True, exist_ok=True)
        
        pdb_path = UPLOAD_FOLDER / f"{job_id}.pdb"
        file.save(pdb_path)
        
        try:
            # Configure style for matplotlib
            protein_analysis.configure_style()
            
            # Run analysis
            residues, atoms = protein_analysis.load_residues(pdb_path)
            results = protein_analysis.analyze_structure(residues, atoms)
            
            # Generate plots
            prefix = "protein"
            protein_analysis.plot_3d_backbone(
                results["ca_coords"],
                results["residue_scores"],
                job_dir / f"{prefix}_3d_backbone.png",
                dpi=150 # Lower DPI for web display speed
            )
            protein_analysis.plot_residue_score(
                results["residue_scores"],
                job_dir / f"{prefix}_residue_score.png",
                dpi=150
            )
            protein_analysis.plot_contact_map(
                results["distance_map"],
                job_dir / f"{prefix}_ca_distance_map.png",
                dpi=150
            )
            protein_analysis.plot_ramachandran(
                results["phi"],
                results["psi"],
                job_dir / f"{prefix}_ramachandran.png",
                dpi=150
            )
            protein_analysis.plot_overview_panel(
                results["residue_scores"],
                results["summary"]["aa_composition"],
                job_dir / f"{prefix}_overview.png",
                dpi=150
            )
            
            # Generate Markdown report
            report_path = job_dir / f"{prefix}_report.md"
            protein_analysis.generate_markdown_report(results, report_path, prefix)
            
            # Create a zip file for download
            zip_path = RESULTS_FOLDER / f"{job_id}.zip"
            with zipfile.ZipFile(zip_path, 'w') as zipf:
                for root, dirs, files in os.walk(job_dir):
                    for f in files:
                        zipf.write(os.path.join(root, f), f)
            
            return redirect(url_for('result', job_id=job_id))
            
        except Exception as e:
            flash(f"Error during analysis: {str(e)}")
            return redirect(url_for('index'))
    else:
        flash('Allowed file type is .pdb')
        return redirect(request.url)

@app.route('/result/<job_id>')
def result(job_id):
    job_dir = RESULTS_FOLDER / job_id
    if not job_dir.exists():
        flash('Result not found')
        return redirect(url_for('index'))
    
    # List files in job_dir to display
    files = os.listdir(job_dir)
    images = [f for f in files if f.endswith('.png')]
    report_file = next((f for f in files if f.endswith('_report.md')), None)
    
    report_html = ""
    if report_file:
        with open(job_dir / report_file, 'r') as f:
            content = f.read()
            # Convert markdown to HTML
            # Adjust image paths in markdown to point to static files
            content = content.replace('](' + "protein", '](' + url_for('static', filename=f'results/{job_id}/protein'))
            report_html = markdown.markdown(content, extensions=['extra'])
            
    return render_template('result.html', job_id=job_id, images=images, report_html=report_html)

@app.route('/download/<job_id>')
def download(job_id):
    zip_filename = f"{job_id}.zip"
    return send_from_directory(RESULTS_FOLDER, zip_filename, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True, port=5000)
