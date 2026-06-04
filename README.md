# mof-guest-toolkit

A Python package of helper utilities for computational MOF–guest chemistry workflows.
Developed as part of the thesis *"Impact of Framework Topology on the Selective Separation
of Pharmaceuticals and Cannabinoids in Metal-Organic Frameworks"* (TU Dresden, 2025).

The package provides command-line tools and importable functions for:

- Fetching and visualising molecular structures from PubChem
- Batch-computing RDKit molecular descriptors
- Retrieving CIF files from the CCDC by refcode *(coming soon)*
- Parsing geometry optimisation outputs from AMS/ORCA *(coming soon)*
- Analysing normal modes and gradient convergence *(coming soon)*

---

## Requirements

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- Git
- An internet connection (for PubChem API calls)

> **Why conda?**  
> This package depends on [RDKit](https://www.rdkit.org/), which must be installed
> via `conda-forge`. A plain `pip install` will not work for RDKit.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/adricu12/mof-guest-toolkit.git
cd mof-guest-toolkit
```

### 2. Create the conda environment

This single command creates a dedicated environment called `mof-toolkit`,
installs all dependencies (including RDKit from conda-forge), and installs
this package itself in editable mode so CLI commands are immediately available.

```bash
conda env create -f environment.yml
```

### 3. Activate the environment

```bash
conda activate mof-toolkit
```

You need to activate the environment every time you open a new terminal session.
To deactivate: `conda deactivate`.

### 4. Verify the installation

```bash
pubchem_check 3033 rdMolDescriptors.CalcTPSA
```

Expected output:
```
Chem.rdMolDescriptors.CalcTPSA  (CID 3033):
  81.0799...
```

---

## Updating

If the code changes (e.g. after `git pull`), you do **not** need to reinstall —
the editable install (`-e .`) means changes are picked up automatically.

If new dependencies are added to `environment.yml`:

```bash
conda env update -f environment.yml --prune
```

---

## Usage

### `pubchem_interactive`

Opens an interactive Jupyter notebook in your browser where you can enter a CID
or molecule name and see the 3D structure and a descriptor table.

```bash
pubchem_interactive
```

> **WSL users:** Jupyter will not open a browser automatically.
> Copy the `localhost:8888/?token=...` URL printed in the terminal
> and paste it into your Windows browser.

---

### `pubchem_check`

Evaluate any single RDKit function on a compound and print the result.

```bash
pubchem_check <cid_or_name> <rdkit_function>
```

```bash
# by CID
pubchem_check 3033 rdMolDescriptors.CalcTPSA

# by name
pubchem_check aspirin rdMolDescriptors.CalcNumAromaticRings

# fragment counts
pubchem_check 644019 Fragments.fr_Ar_OH
```

Supported namespaces: `Chem`, `rdMolDescriptors`, `Fragments`, `Descriptors`, `GraphDescriptors`.

If the function requires extra arguments beyond `mol`, a clear error message is printed.

---

### `pubchem_batch_fetcher`

Compute descriptors for a list of compounds from a CSV file.

```bash
pubchem_batch_fetcher <input.csv> <output.csv>
```

**Minimum input CSV** (only `CID` is required):

```
CID
3033
644019
2244
```

**Full input CSV** (metadata columns are passed through to the output):

```
CID,Name,Abbreviation,Guest_Type
3033,Diclofenac,DIC,Pharmaceutical
644019,Cannabidiol,CBD,Cannabinoid
```

**With extra property columns** (use `::` to specify a custom RDKit function):

```
CID,Name,Guest_Type,Chi0v::Chem.rdMolDescriptors.CalcChi0v
3033,Diclofenac,Pharmaceutical
```

The output CSV contains all metadata columns plus all default descriptors
(MolecularWeight, HBA, HBD, RotatableBonds, TPSA, NumRings, NumAromaticRings,
and fragment counts) plus any extra computed columns.

A ready-to-use input file for the 37 pharmaceutical and cannabinoid guests
from the thesis is included at `molecules.csv`.

```bash
pubchem_batch_fetcher molecules.csv results.csv
```

---

### `get_xyz_cid`

Print the default descriptor table for one compound in the terminal.

```bash
get_xyz_cid <cid_or_name>

get_xyz_cid 3033
get_xyz_cid aspirin
```

---

### `fetch_xyz_batch`

Download PubChem 3D conformers as `.xyz` files for a list of compounds.

```bash
fetch_xyz_batch <input.csv> <output_dir>
```

The input CSV must have `CID` and `Abbreviation` columns.
Each compound is saved as `<Abbreviation>.xyz` in the output directory.

```bash
fetch_xyz_batch molecules.csv ./xyz_files/
```

---

## Using functions in your own scripts

All functions are importable directly:

```python
from mof_toolkit.pubchem import (
    compute_properties,
    fetch_smiles_from_cid,
    fetch_and_save_xyz,
    show_molecule,          # for use inside Jupyter notebooks
)

# compute default descriptors for Diclofenac
props = compute_properties(3033)
print(props["TPSA"])

# get SMILES
smiles = fetch_smiles_from_cid(644019)
print(smiles)
```

---

## Default descriptor set

| Property | RDKit function |
|---|---|
| MolecularWeight | `CalcExactMolWt` |
| NumRings | `CalcNumRings` |
| NumAromaticRings | `CalcNumAromaticRings` |
| HBA | `CalcNumHBA` |
| HBD | `CalcNumHBD` |
| RotatableBonds | `CalcNumRotatableBonds` |
| TPSA | `CalcTPSA` |
| fr_Al_OH | `Fragments.fr_Al_OH` |
| fr_Ar_OH | `Fragments.fr_Ar_OH` |
| fr_COO | `Fragments.fr_COO` |
| fr_C_O_noCOO | `Fragments.fr_C_O_noCOO` |
| fr_Ar_N | `Fragments.fr_Ar_N` |
| fr_NH2 / fr_NH1 / fr_NH0 | `Fragments.fr_NH*` |
| fr_ether | `Fragments.fr_ether` |
| fr_sulfonamd | `Fragments.fr_sulfonamd` |

---

## Project structure

```
mof-guest-toolkit/
├── mof_toolkit/
│   ├── __init__.py
│   ├── pubchem.py          ← PubChem fetchers, RDKit descriptors, CLI tools
│   ├── ccdc.py             ← CIF fetcher by CCDC refcode        [coming soon]
│   ├── geometry_opt.py     ← AMS/ORCA .out parser               [coming soon]
│   ├── normal_modes.py     ← frequency and imaginary mode tools  [coming soon]
│   └── bonding.py          ← connectivity and bond distance checks [coming soon]
├── molecules.csv           ← 37 thesis guest molecules, ready to use
├── pyproject.toml          ← package metadata and CLI entry points
├── environment.yml         ← conda environment definition
└── README.md
```

---

## Citation

If you use this package in your work, please cite:

> Ugarte, A. (2025). *Impact of Framework Topology on the Selective Separation
> of Pharmaceuticals and Cannabinoids in Metal-Organic Frameworks.*
> Master's Thesis, TU Dresden.
> Code: https://github.com/YOUR_USERNAME/mof-guest-toolkit

---

## License

MIT License. See `LICENSE` for details.