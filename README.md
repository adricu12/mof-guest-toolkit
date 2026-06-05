# <img src="docs/figures/main-logo.png" alt="mascot" height="80" style="vertical-align:middle"> MOF-Guest Toolkit

`mof-guest-toolkit` is a Python package of reusable utilities for computational
chemistry workflows, with an emphasis on MOF–guest systems, cheminformatics,
molecular data handling, and research automation.

The package originated from scripts developed during a Master's thesis in
Theoretical Chemistry at TU Dresden and has been reorganised into a structured
toolkit for research, coursework, and independent computational projects.

---

<img src="docs/figures/searching.png" alt="cheminformatics" height="150" align="right" hspace="16">

**Cheminformatics**

Resolve compound identifiers (CID, name, formula, SMILES, InChIKey) to
canonical SMILES and rich metadata. Compute RDKit descriptors for single
compounds or large batches. Generate 3D conformers from PubChem data or
locally with RDKit. Visualise molecules interactively in the browser.

- [Interactive 3D viewer](./examples.md/#example-1--interactive-3d-viewer)
- [Resolve a compound identifier](./examples.md/#example-2--resolve-a-compound-identifier)
- [Compute descriptors for one compound](./examples.md/#example-3--compute-descriptors-for-one-compound)
- [Batch descriptor computation](./examples.md/#example-4--batch-descriptor-computation)
- [Generate 3D structure files](./examples.md/#example-5--generate-3d-structure-files)
- [Convert between structure formats](./examples.md/#example-6--convert-between-structure-formats)
- [Low-level Python helpers](./examples.md/#example-7--low-level-python-helpers)

<br clear="both">

---

<img src="docs/figures/computer.png" alt="HPC" height="140" align="right" hspace="16">

**HPC & Simulation Setup**

Generate submission scripts for HPC clusters. Restart failed or incomplete
geometry optimisations. Build topology and input files for LAMMPS MD simulations.

<br clear="both">

---

<img src="docs/figures/thinking2.png" alt="analysis" height="140" align="right" hspace="16">

**Computational Output Analysis**

Parse geometry optimisation outputs from AMS and ORCA. Analyse normal modes,
imaginary frequencies, and gradient convergence. Process LAMMPS trajectory files
and run automated structure checks.

<br clear="both">

---

<img src="docs/figures/lab.png" alt="wet lab" height="120" align="right" hspace="16">

**Wet Lab Helpers**

Solution preparation calculator. UV–Vis calibration curves. Standardised PXRD
analysis and laboratory plots.

<br clear="both">

---

## Requirements

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda (recommended)
- Git
- An internet connection (for PubChem and EBI API calls)

> **Note:** The conda install path is recommended because it pulls RDKit from
> `conda-forge`, which provides the most reliable binary dependencies.
> A pure `pip install` is also supported (RDKit has been pip-installable since 2022).

---

## Installation

### Option A — conda (recommended)

```bash
git clone https://github.com/adricu12/mof-guest-toolkit.git
cd mof-guest-toolkit
conda env create -f environment.yml
conda activate mof-toolkit
```

### Option B — pip only

```bash
git clone https://github.com/adricu12/mof-guest-toolkit.git
cd mof-guest-toolkit
pip install -e .
```

Activate the conda environment in every new terminal session.
To deactivate: `conda deactivate`.

---

### Verify the installation

```bash
mol_explorer -input 3033
```

Expected output:

```
Descriptor                       Value
--------------------------------------------------------------
CID                              3033
IUPAC_Name                       2-(2,6-dichloroanilino)phenylacetic acid
Common_Name                      Diclofenac
SMILES                           OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl
Formula                          C14H11Cl2NO2
InChIKey                         DCOPUUMXTXDBNB-UHFFFAOYSA-N
```

---

## Updating

```bash
cd mof-guest-toolkit
conda activate mof-toolkit
git pull
conda env update -f environment.yml --prune
pip install -e .
```

---

## Command reference

### Cheminformatics CLI

| Command | Description |
|---|---|
| `interactive_explorer` | Local web viewer — rotatable 3D structure + descriptor table + EBI cross-refs |
| `mol_explorer -input <query> [-format <type>]` | Resolve a compound and print its identifiers |
| `mol_explorer -input <query> -descriptor <default\|full\|NAME…>` | Compute descriptors for one compound |
| `mol_explorer -input <query> -descriptor <…> -output <file>` | Save descriptors to CSV/TXT |
| `mol_explorer -batch <csv> -descriptor <…> -output <file>` | Batch-compute descriptors from a CSV |
| `mol_get_xyz -input <query> [-outputformat <xyz sdf pdb mol>]` | Generate 3D structure file(s) for one compound |
| `mol_get_xyz -batch <csv> [-outputformat <…>] -output <dir>` | Batch-generate 3D structure files from a CSV |
| `mol_file_translate -input <file> -output <file>` | Convert between xyz / sdf / mol / pdb formats |

All commands support `--help` for full usage details.

**Supported input types for `mol_explorer` and `mol_get_xyz`** (pass via `-format`):

| Value | Accepts |
|---|---|
| `cid` | PubChem CID integer |
| `name` | IUPAC or common name |
| `formula` | Molecular formula (e.g. `C9H8O4`, `c6h6`, `NaHCO3`) |
| `smiles` | SMILES string |
| `inchikey` | 27-character InChIKey |

### Python helpers

| Function | Module | Description |
|---|---|---|
| `resolve_compound_input(query, input_type)` | `rdkit_descriptors` | Resolve any identifier to a unified metadata dict |
| `get_rdkit_dict(query, ...)` | `rdkit_descriptors` | Full descriptor dict for one compound |
| `get_descriptor(smiles, names)` | `rdkit_descriptors` | Compute named RDKit descriptors from SMILES |
| `get_smiles(query, input_type)` | `rdkit_descriptors` | Resolve any identifier to canonical SMILES |
| `get_all_rdkit_descriptors(mol)` | `rdkit_descriptors` | All ~210 RDKit descriptors for an RDKit Mol |
| `display_table(dict)` | `rdkit_descriptors` | Print a descriptor dict as an aligned table |
| `get_3d_structure(query, ...)` | `molecule_manager` | Generate 3D structure file(s) for one compound |
| `get_smiles_from_coords(filepath)` | `molecule_manager` | Read canonical SMILES from a 3D structure file |
| `embed_3d(mol)` | `molecule_manager` | RDKit ETKDGv3 + MMFF94 conformer for an RDKit Mol |

---

## Project structure

```
mof-guest-toolkit/
├── mof_toolkit/
│   ├── __init__.py               ← public API re-exports
│   ├── rdkit_descriptors.py      ← resolver, descriptors, mol_explorer CLI, Flask viewer
│   ├── molecule_manager.py       ← 3D conformer generation, mol_get_xyz, mol_file_translate
│   ├── templates/
│   │   └── viewer.html           ← Flask/3Dmol.js interactive viewer template
│   ├── static/
│   │   └── figures/              ← mascot and UI images
│   ├── ccdc.py                   ← [coming soon]
│   ├── geometry_opt.py           ← [coming soon]
│   ├── normal_modes.py           ← [coming soon]
│   └── bonding.py                ← [coming soon]
├── docs/
│   └── figures/
├── tests/
│   └── data/
│       └── molecules_example.csv ← 6-compound demo CSV
├── pyproject.toml                ← package metadata and CLI entry points
├── environment.yml               ← conda environment definition
├── examples.md                   ← worked examples for all features
└── README.md
```

---

## Citation

If you use this package in your work, please cite:

> Ugarte, A. (2025). *Impact of Framework Topology on the Selective Separation
> of Pharmaceuticals and Cannabinoids in Metal-Organic Frameworks.*
> Master's Thesis, TU Dresden.
> Code: https://github.com/adricu12/mof-guest-toolkit

---

## License

MIT License. See `LICENSE` for details.