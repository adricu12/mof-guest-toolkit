# <img src="docs/figures/main-logo.png" alt="mascot" height="80" style="vertical-align:middle"> MOF-Guest-Toolkit

`mof-guest-toolkit` is a Python package of reusable utilities for computational chemistry workflows, with an emphasis on MOF–guest systems, cheminformatics, molecular data handling, and research automation.

The package originated from scripts developed during a Master's thesis in Theoretical Chemistry at TU Dresden and has been reorganized into a structured toolkit for research, coursework, and independent computational projects. It provides command-line tools and importable Python functions for molecular structure retrieval, descriptor calculation, 3D conformer generation, and workflow support.


<img src="docs/figures/searching.png" alt="cheminformatics" height="150" align="right" hspace="16">

**Cheminformatics**

Fetch and visualise molecular structures from PubChem. Compute RDKit descriptors for
single compounds or large batches. Generate 3D conformers from SMILES or PubChem data.

- [Quick interactive exploration of molecules](./examples.md/#example-1--interactive-3d-viewer)
- [Computing a specific RDKit descriptor for a molecule](./examples.md/#example-2--check-a-single-rdkit-descriptor)
- [Getting descriptors for a set of molecules](./examples.md/#example-4--batch-descriptor-computation)
- [Getting 3D coordinate files of molecules](./examples.md/#example-5--generate-3d-structure-files)

<br clear="both">

---

<img src="docs/figures/computer.png" alt="HPC" height="140" align="right" hspace="16">

**HPC & Simulation Setup**

Generate submission scripts for HPC clusters. Restart failed or incomplete geometry
optimisations. Build topology and input files for LAMMPS MD simulations of molecules
and MOFs.

- Generating HPC submission scripts
- Restarting failed or incomplete geometry optimisations
- Generating topology files for LAMMPS MD simulations

<br clear="both">

---

<img src="docs/figures/thinking2.png" alt="analysis" height="140" align="right" hspace="16">

**Computational Output Analysis**

Parse geometry optimisation outputs from AMS and ORCA. Analyse normal modes, imaginary
frequencies, and gradient convergence. Process LAMMPS trajectory files and run
automated structure checks.

- Parsing geometry optimisation outputs from AMS/ORCA
- Analysing normal modes and gradient convergence
- Analysing LAMMPS trajectory files
- Structure checker

<br clear="both">

---

<img src="docs/figures/lab.png" alt="wet lab" height="120" align="right" hspace="16">

**Wet Lab Helpers**

Solution preparation calculator for concentrations and dilutions. UV–Vis calibration
curves — input absorbance values and retrieve concentrations directly. Standardised
PXRD analysis and laboratory plots.

- Solution preparation calculator (concentrations, dilutions)
- UV–Vis calibration tools
- PXRD analysis and standardised laboratory plots

<br clear="both">

---

## Requirements

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda (recommended)
- Git
- An internet connection (for PubChem API calls)

> **Note:** The conda install path is recommended because it pulls RDKit from `conda-forge`,
> which provides the most reliable binary dependencies.
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

Activate the conda environment every time you open a new terminal session.
To deactivate: `conda deactivate`.

---

### Verify the installation

```bash
rdkit_check_descrpt 3033 TPSA
```

Expected output:

```
CID         : 3033
IUPAC_Name  : 2-(2,6-dichloroanilino)phenylacetic acid
Common_Name : Diclofenac
SMILES      : OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl
Function    : TPSA
Value       : 49.330000
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

| Command | Description |
|---|---|
| `pubchem_interactive` | Local web viewer — rotatable 3D structure + descriptor table |
| `rdkit_check_descrpt <cid\|name\|smiles> <descriptor>` | Evaluate one RDKit descriptor on a compound |
| `rdkit_default_descrpts <cid\|name\|smiles>` | Print all default descriptors for a compound |
| `rdkit_batch_fetcher <in.csv> <out.csv>` | Batch-compute descriptors from a CSV file |
| `fetch_xyz_batch <in.csv> <out_dir>` | Batch-generate 3D structure files from a CSV file |
| `smiles_to_3d <SMILES>` | Generate 3D structure file(s) locally from a SMILES string |

All commands support `--help` for full usage details and available options.

**Python-only helpers** (for use in scripts and notebooks):

| Function | Module | Description |
|---|---|---|
| `get_rdkit_dict(query)` | `rdkit_descriptors` | Returns a descriptor dictionary for one compound |
| `resolve_compound_input(query)` | `rdkit_descriptors` | Resolves CID, name, or SMILES to a unified dict |
| `get_3d_structure(query, ...)` | `molecule_manager` | Generates 3D structure file(s) for one compound |
| `embed_3d(mol)` | `molecule_manager` | RDKit ETKDGv3 + MMFF94 conformer for an RDKit Mol object |

---

## Project structure

```
mof-guest-toolkit/
├── mof_toolkit/
│   ├── __init__.py
│   ├── rdkit_descriptors.py  ← PubChem fetchers, RDKit descriptors, CLI tools
│   ├── molecule_manager.py   ← 3D conformer generation and structure file writing
│   ├── static/
│   │   └── figures/          
│   ├── ccdc.py               ← [coming soon]
│   ├── geometry_opt.py       ← [coming soon]
│   ├── normal_modes.py       ← [coming soon]
│   └── bonding.py            ← [coming soon]
├── docs/
│   └── figures/
├── tests/
│   └── data/
│       └── molecules_example.csv   ← 6-compound demo file
├── pyproject.toml            ← package metadata and CLI entry points
├── environment.yml           ← conda environment definition
├── examples.md
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
