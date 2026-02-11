#!/usr/bin/env python3
"""Generic single-protein PDB analysis with publication-style visualizations."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


@dataclass
class AtomRecord:
    atom_name: str
    resname: str
    chain: str
    resseq: int
    icode: str
    coord: np.ndarray
    bfactor: float
    occupancy: float
    element: str


@dataclass
class ResidueRecord:
    resname: str
    chain: str
    resseq: int
    icode: str
    atoms: Dict[str, AtomRecord]

    @property
    def residue_id(self) -> str:
        suffix = self.icode.strip()
        return f"{self.chain}:{self.resseq}{suffix}" if suffix else f"{self.chain}:{self.resseq}"


# ------------------------- Parsing -------------------------

def parse_atom_line(line: str) -> Optional[AtomRecord]:
    if not line.startswith(("ATOM", "HETATM")):
        return None

    altloc = line[16].strip()
    if altloc not in ("", "A"):
        return None

    atom_name = line[12:16].strip()
    resname = line[17:20].strip()
    chain = line[21].strip() or "_"
    resseq = int(line[22:26])
    icode = line[26].strip()
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    occupancy = float(line[54:60])
    bfactor = float(line[60:66])
    element = line[76:78].strip() if len(line) >= 78 else ""

    return AtomRecord(
        atom_name=atom_name,
        resname=resname,
        chain=chain,
        resseq=resseq,
        icode=icode,
        coord=np.array([x, y, z], dtype=float),
        bfactor=bfactor,
        occupancy=occupancy,
        element=element,
    )


def load_residues(pdb_path: Path, chain_filter: Optional[str] = None) -> Tuple[List[ResidueRecord], List[AtomRecord]]:
    residues: Dict[Tuple[str, int, str], ResidueRecord] = {}
    atoms: List[AtomRecord] = []

    with pdb_path.open("r", encoding="utf-8") as f:
        for line in f:
            rec = parse_atom_line(line)
            if rec is None:
                continue
            if chain_filter and rec.chain != chain_filter:
                continue

            atoms.append(rec)
            key = (rec.chain, rec.resseq, rec.icode)
            if key not in residues:
                residues[key] = ResidueRecord(
                    resname=rec.resname,
                    chain=rec.chain,
                    resseq=rec.resseq,
                    icode=rec.icode,
                    atoms={},
                )
            residues[key].atoms[rec.atom_name] = rec

    ordered = [residues[k] for k in sorted(residues.keys(), key=lambda x: (x[0], x[1], x[2]))]
    if not ordered:
        raise ValueError("No residues parsed. Check chain selection or PDB format.")
    return ordered, atoms


# ------------------------- Geometry -------------------------

def dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def compute_rg(atoms: List[AtomRecord]) -> float:
    coords = np.array([a.coord for a in atoms], dtype=float)
    center = np.mean(coords, axis=0)
    return float(np.sqrt(np.mean(np.sum((coords - center) ** 2, axis=1))))


# ------------------------- Analysis -------------------------

def analyze_structure(residues: List[ResidueRecord], atoms: List[AtomRecord]) -> Dict[str, object]:
    ca_residues: List[ResidueRecord] = [r for r in residues if "CA" in r.atoms]
    ca_coords = np.array([r.atoms["CA"].coord for r in ca_residues], dtype=float)
    residue_scores = np.array([r.atoms["CA"].bfactor for r in ca_residues], dtype=float)

    if len(ca_coords) == 0:
        raise ValueError("No CA atoms found; cannot perform residue-level analysis.")

    diff = ca_coords[:, None, :] - ca_coords[None, :, :]
    ca_distance_map = np.linalg.norm(diff, axis=2)

    sequence = "".join(AA3_TO_1.get(r.resname, "X") for r in ca_residues)

    aa_counts: Dict[str, int] = {}
    for aa in sequence:
        aa_counts[aa] = aa_counts.get(aa, 0) + 1

    phi: List[float] = []
    psi: List[float] = []
    for i in range(1, len(ca_residues) - 1):
        prev_res = ca_residues[i - 1]
        curr_res = ca_residues[i]
        next_res = ca_residues[i + 1]

        needed_prev = {"C"}
        needed_curr = {"N", "CA", "C"}
        needed_next = {"N"}

        if not needed_prev.issubset(prev_res.atoms) or not needed_curr.issubset(curr_res.atoms) or not needed_next.issubset(next_res.atoms):
            continue

        phi_val = dihedral(
            prev_res.atoms["C"].coord,
            curr_res.atoms["N"].coord,
            curr_res.atoms["CA"].coord,
            curr_res.atoms["C"].coord,
        )
        psi_val = dihedral(
            curr_res.atoms["N"].coord,
            curr_res.atoms["CA"].coord,
            curr_res.atoms["C"].coord,
            next_res.atoms["N"].coord,
        )
        phi.append(phi_val)
        psi.append(psi_val)

    chain_ids = sorted({r.chain for r in residues})

    summary = {
        "chain_ids": chain_ids,
        "residue_count": len(ca_residues),
        "atom_count": len(atoms),
        "sequence": sequence,
        "aa_composition": aa_counts,
        "radius_of_gyration": compute_rg(atoms),
        "score_stats": {
            "mean": float(np.mean(residue_scores)),
            "std": float(np.std(residue_scores)),
            "min": float(np.min(residue_scores)),
            "max": float(np.max(residue_scores)),
        },
    }

    return {
        "summary": summary,
        "ca_residues": ca_residues,
        "ca_coords": ca_coords,
        "residue_scores": residue_scores,
        "distance_map": ca_distance_map,
        "phi": np.array(phi, dtype=float),
        "psi": np.array(psi, dtype=float),
    }


# ------------------------- Plotting -------------------------

def configure_style() -> None:
    plt.style.use("seaborn-v0_8-paper")
    plt.rcParams.update(
        {
            "figure.figsize": (8.0, 6.0),
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman", "Liberation Serif"],
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "axes.grid": True,
            "grid.alpha": 0.3,
            "grid.linestyle": "--",
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 5,
            "ytick.major.size": 5,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def plot_3d_backbone(coords: np.ndarray, score: np.ndarray, out_file: Path, dpi: int) -> None:
    fig = plt.figure(figsize=(8.5, 6.5))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], color="#546E7A", linewidth=1.2, alpha=0.9)
    sc = ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=score, cmap="viridis", s=20, alpha=0.95)

    ax.set_title("Protein CA Backbone (color: residue score)")
    ax.set_xlabel("X (A)")
    ax.set_ylabel("Y (A)")
    ax.set_zlabel("Z (A)")
    cb = plt.colorbar(sc, ax=ax, pad=0.08)
    cb.set_label("CA B-factor / confidence")

    fig.savefig(out_file, dpi=dpi)
    plt.close(fig)


def plot_residue_score(score: np.ndarray, out_file: Path, dpi: int) -> None:
    x = np.arange(1, len(score) + 1)
    mean = np.mean(score)

    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.plot(x, score, color="#0D47A1", linewidth=1.8)
    ax.fill_between(x, score, mean, color="#90CAF9", alpha=0.3)
    ax.axhline(mean, color="#D81B60", linestyle="--", linewidth=1.3, label=f"Mean = {mean:.2f}")

    ax.set_title("Per-residue score profile")
    ax.set_xlabel("Residue index")
    ax.set_ylabel("CA B-factor / confidence")
    ax.legend(loc="best", frameon=False)

    fig.savefig(out_file, dpi=dpi)
    plt.close(fig)


def plot_contact_map(distance_map: np.ndarray, out_file: Path, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    # Using 'viridis_r' or 'GnBu_r' for a cleaner look
    im = ax.imshow(distance_map, cmap="viridis_r", origin="lower", extent=[1, distance_map.shape[0], 1, distance_map.shape[1]])

    ax.set_title("CA-CA Contact Map", fontweight="bold", pad=12)
    ax.set_xlabel("Residue Index")
    ax.set_ylabel("Residue Index")
    
    # Add minor ticks for better readability
    from matplotlib.ticker import MultipleLocator
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.yaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(10))

    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(r"Distance ($\AA$)")

    fig.savefig(out_file, dpi=dpi)
    plt.close(fig)


def plot_ramachandran(phi: np.ndarray, psi: np.ndarray, out_file: Path, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 5.8))
    
    # Standard favored regions (approximate)
    from matplotlib.patches import Rectangle
    regions = [
        # Beta sheets
        {"xy": (-150, 100), "w": 100, "h": 60, "label": "Beta"},
        # Alpha helix (Right)
        {"xy": (-100, -60), "w": 60, "h": 40, "label": "Alpha-R"},
        # Alpha helix (Left)
        {"xy": (40, 20), "w": 60, "h": 60, "label": "Alpha-L"},
    ]
    
    for r in regions:
        ax.add_patch(Rectangle(r["xy"], r["w"], r["h"], color="gray", alpha=0.1, zorder=0))
    
    ax.scatter(phi, psi, s=20, alpha=0.6, color="#1565C0", edgecolors="white", linewidth=0.5, zorder=2)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
    ax.set_yticks([-180, -120, -60, 0, 60, 120, 180])
    
    ax.axhline(0, color="0.6", linewidth=0.8, linestyle="--", zorder=1)
    ax.axvline(0, color="0.6", linewidth=0.8, linestyle="--", zorder=1)

    ax.set_title("Ramachandran Plot", fontweight="bold", pad=15)
    ax.set_xlabel(r"$\phi$ (degrees)")
    ax.set_ylabel(r"$\psi$ (degrees)")

    fig.savefig(out_file, dpi=dpi)
    plt.close(fig)


def plot_overview_panel(score: np.ndarray, aa_comp: Dict[str, int], out_file: Path, dpi: int) -> None:
    x = np.arange(1, len(score) + 1)
    aa_items = sorted(aa_comp.items(), key=lambda kv: kv[1], reverse=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Profile Plot
    axes[0].plot(x, score, color="#1976D2", linewidth=1.5)
    axes[0].fill_between(x, score, alpha=0.1, color="#1976D2")
    axes[0].set_title("Residue Score Profile", fontweight="bold")
    axes[0].set_xlabel("Residue Index")
    axes[0].set_ylabel("Score / B-factor")

    # Property-based coloring for AA composition
    # Hydrophobic: A,V,I,L,M,F,W,P -> #2196F3
    # Polar: S,T,C,Y,N,Q -> #4CAF50
    # Positive: K,R,H -> #F44336
    # Negative: D,E -> #FF9800
    prop_colors = {
        "A": "#2196F3", "V": "#2196F3", "I": "#2196F3", "L": "#2196F3", "M": "#2196F3", "F": "#2196F3", "W": "#2196F3", "P": "#2196F3",
        "S": "#4CAF50", "T": "#4CAF50", "C": "#4CAF50", "Y": "#4CAF50", "N": "#4CAF50", "Q": "#4CAF50",
        "K": "#F44336", "R": "#F44336", "H": "#F44336",
        "D": "#FF9800", "E": "#FF9800", "G": "#9E9E9E"
    }

    labels = [k for k, _ in aa_items]
    values = [v for _, v in aa_items]
    colors = [prop_colors.get(l, "#9E9E9E") for l in labels]
    
    axes[1].bar(labels, values, color=colors, alpha=0.85, edgecolor="gray", linewidth=0.5)
    axes[1].set_title("Amino Acid Composition", fontweight="bold")
    axes[1].set_xlabel("Residue Type")
    axes[1].set_ylabel("Frequency")

    # Add legend for properties
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color="#2196F3", lw=4, label="Hydrophobic"),
        Line2D([0], [0], color="#4CAF50", lw=4, label="Polar"),
        Line2D([0], [0], color="#F44336", lw=4, label="Positive"),
        Line2D([0], [0], color="#FF9800", lw=4, label="Negative"),
    ]
    axes[1].legend(handles=legend_elements, loc="upper right", fontsize=8, frameon=False)

    fig.tight_layout()
    fig.savefig(out_file, dpi=dpi)
    plt.close(fig)


def generate_markdown_report(results: Dict[str, object], out_file: Path, prefix: str) -> None:
    summary = results["summary"]
    stats = summary["score_stats"]
    seq = summary["sequence"]

    # Calculate some extra metrics
    cys_count = summary["aa_composition"].get("C", 0)
    cys_perc = (cys_count / len(seq)) * 100 if len(seq) > 0 else 0

    # Theoretical Rg (approximation for globular proteins)
    # Rg ~ 0.395 * N^(3/5)
    theoretical_rg = 0.395 * (len(seq) ** 0.6)

    md = f"""# Protein Structural Bioinformatics Analysis Report (ID: {prefix})

## 1. Executive Summary
This report presents a systematic analysis of the protein structure provided in the PDB file. The protein consists of {summary['residue_count']} residues across {len(summary['chain_ids'])} chain(s) ({', '.join(summary['chain_ids'])}). 

## 2. Sequence Composition & Biochemical Properties
- **Total Residues**: {summary['residue_count']} AA
- **Atom Count**: {summary['atom_count']}
- **Cysteine Content**: {cys_count} residues ({cys_perc:.2f}%)
  *The {'high' if cys_perc > 3 else 'moderate'} percentage of Cysteine suggests the {'potential' if cys_perc > 0 else 'absence of'} presence of disulfide bonds, which are critical for structural stability.*

![Amino Acid Composition]({prefix}_overview.png)

## 3. Compactness & Folding State
- **Radius of Gyration ($R_g$)**: **{summary['radius_of_gyration']:.2f} Å**
- **Theoretical $R_g$ (globular)**: ~{theoretical_rg:.2f} Å
- **Analysis**: The measured $R_g$ of {summary['radius_of_gyration']:.2f} Å indicates a {'highly compact' if summary['radius_of_gyration'] < theoretical_rg + 2 else 'slightly extended'} fold. 

![Contact Map]({prefix}_ca_distance_map.png)

## 4. Conformational Analysis (Ramachandran Plot)
The Ramachandran plot shows the distribution of dihedral angles ($\\phi$ and $\\psi$). Most residues are expected to fall within the favored $\\alpha$-helical and $\\beta$-sheet regions.

![Ramachandran Plot]({prefix}_ramachandran.png)

## 5. Flexibility & Stability Assessment (B-factor)
- **Mean B-factor**: {stats['mean']:.2f}
- **B-factor Range**: {stats['min']:.2f} - {stats['max']:.2f}
- **Analysis**: Lower B-factor values indicate more stable regions (typically the hydrophobic core), while higher values correspond to flexible loops or termini.

![B-factor Profile]({prefix}_residue_score.png)
![3D Backbone]({prefix}_3d_backbone.png)

## 6. Conclusion
The protein (ID: {prefix}) exhibits a well-defined fold with a clear stable core. {'The presence of ' + str(cys_count) + ' Cysteines suggests that disulfide bonding analysis could be a fruitful area for further research.' if cys_count > 1 else ''}

---
*Report generated automatically by protein_analysis.py*
"""
    out_file.write_text(md, encoding="utf-8")


# ------------------------- Main -------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Single-protein PDB analysis with reusable interface.")
    p.add_argument("--pdb", required=True, type=Path, help="Input PDB file path")
    p.add_argument("--chain", default=None, help="Chain ID to analyze, e.g. A. Default: all chains")
    p.add_argument("--output-dir", type=Path, default=Path("results"), help="Output directory")
    p.add_argument("--prefix", default="protein", help="Output file prefix")
    p.add_argument("--dpi", type=int, default=300, help="Figure DPI")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    configure_style()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    residues, atoms = load_residues(args.pdb, args.chain)
    results = analyze_structure(residues, atoms)

    summary_path = args.output_dir / f"{args.prefix}_summary.json"
    summary_path.write_text(json.dumps(results["summary"], indent=2, ensure_ascii=False), encoding="utf-8")

    plot_3d_backbone(
        results["ca_coords"],
        results["residue_scores"],
        args.output_dir / f"{args.prefix}_3d_backbone.png",
        args.dpi,
    )
    plot_residue_score(
        results["residue_scores"],
        args.output_dir / f"{args.prefix}_residue_score.png",
        args.dpi,
    )
    plot_contact_map(
        results["distance_map"],
        args.output_dir / f"{args.prefix}_ca_distance_map.png",
        args.dpi,
    )
    plot_ramachandran(
        results["phi"],
        results["psi"],
        args.output_dir / f"{args.prefix}_ramachandran.png",
        args.dpi,
    )
    plot_overview_panel(
        results["residue_scores"],
        results["summary"]["aa_composition"],
        args.output_dir / f"{args.prefix}_overview.png",
        args.dpi,
    )

    # Generate Markdown Report
    report_path = args.output_dir / f"{args.prefix}_report.md"
    generate_markdown_report(results, report_path, args.prefix)

    print(f"Analysis complete. Outputs written to: {args.output_dir}")
    print(f"Summary file: {summary_path}")
    print(f"Markdown Report: {report_path}")


if __name__ == "__main__":
    main()
