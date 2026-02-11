# Protein Structural Profiler

A comprehensive, reusable tool for automated single-protein PDB analysis and publication-quality visualization. Designed for structural biologists and bioinformaticians to quickly profile protein structures, assess stability, and generate integrated research reports.

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioinformatics](https://img.shields.io/badge/domain-Bioinformatics-green.svg)]()

## 🌟 Key Features

- **3D Backbone Visualization**: Generate interactive-style 3D plots of the CA backbone, color-coded by B-factor or confidence scores (e.g., pLDDT).
- **Comprehensive Geometric Analysis**:
  - **Ramachandran Plots**: Assessment of backbone dihedral angles ($\phi$, $\psi$) with favored region overlays.
  - **Contact Maps**: High-resolution CA-CA distance matrices to identify structural domains and folding patterns.
  - **Radius of Gyration ($R_g$)**: Quantify protein compactness with comparisons to theoretical globular models.
- **Biochemical Profiling**: Detailed amino acid composition analysis categorized by physicochemical properties (Hydrophobic, Polar, Acidic, Basic).
- **Automated Reporting**: Generates a professional Markdown report summarizing all structural insights, statistics, and visualizations.
- **Publication-Ready Graphics**: All plots are generated with high DPI (default 300) using a clean, scientific aesthetic.

## 📸 Visualization Showcase

| 3D Backbone (B-factor Color) | Ramachandran Plot |
|:---:|:---:|
| ![3D Backbone](results_52/52_3d_backbone.png) | ![Ramachandran](results_52/52_ramachandran.png) |

| Contact Map (CA-CA Distance) | Overview Panel |
|:---:|:---:|
| ![Contact Map](results_52/52_ca_distance_map.png) | ![Overview](results_52/52_overview.png) |

## 🚀 Quick Start

### 1. Prerequisites

Ensure you have Python 3.8+ installed. Install the required dependencies:

```bash
pip install numpy matplotlib
```

### 2. Basic Usage

Run the analysis script by providing a PDB file and an output directory:

```bash
python protein_analysis.py \
  --pdb path/to/your_protein.pdb \
  --chain A \
  --output-dir ./results_output \
  --prefix my_protein
```

### 3. Command Line Arguments

| Argument | Description | Default |
|:--- |:--- |:--- |
| `--pdb` | **(Required)** Path to the input PDB file. | N/A |
| `--chain` | Specific chain ID to analyze (e.g., A, B). | All chains |
| `--output-dir`| Directory to save the results. | `results` |
| `--prefix` | Prefix for all output filenames. | `protein` |
| `--dpi` | Resolution for output images. | `300` |

## 📊 Output Artifacts

The tool generates a structured output directory containing:

- `*_report.md`: A comprehensive summary report in Markdown format.
- `*_summary.json`: Raw statistical data (sequence, composition, $R_g$, B-factor stats).
- `*_3d_backbone.png`: 3D representation of the protein backbone.
- `*_ramachandran.png`: Dihedral angle distribution.
- `*_ca_distance_map.png`: Heatmap of inter-residue distances.
- `*_residue_score.png`: Per-residue B-factor/confidence profile.
- `*_overview.png`: Combined panel of residue scores and AA composition.

## 🛠️ Methodology

1.  **Parsing**: Robust PDB parsing supporting `ATOM` and `HETATM` records, handling multiple chains and insertion codes.
2.  **Geometry**: Calculation of Euclidean distances, dihedral angles, and the Radius of Gyration.
3.  **Statistics**: Sequence-based analysis including amino acid distribution and property mapping.
4.  **Visualization**: Rendering using `matplotlib` with a custom scientific style for high-impact figures.

## 📝 License

Distributed under the MIT License. See `LICENSE` for more information.

---

*Developed for rapid protein structural bioinformatics profiling.*
