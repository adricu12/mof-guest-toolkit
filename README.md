# Host-Guest-toolkit

A Python package of helper utilities for computational chemistry workflows in the context of host-guest studies in MOFs.
Developed as part of the thesis *"Impact of Framework Topology on the Selective Separation
of Pharmaceuticals and Cannabinoids in Metal-Organic Frameworks"* (TU Dresden, 2025).

The package provides command-line tools and importable functions for:

- Fetching and visualising molecular structures from PubChem library
- Batch-computing RDKit molecular descriptors
- Parsing geometry optimisation outputs from AMS/ORCA *(coming soon)*
- Analysing normal modes and gradient convergence *(coming soon)*

---

## Requirements

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- Git
- An internet connection (for PubChem API calls)


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
pubchem_check_prop 3033 rdMolDescriptors.CalcTPSA
```

Expected output:
```
  Compound : 3033  (CID 3033)
  Function : rdMolDescriptors.CalcTPSA
  Value    : 49.33...
```

---

## Updating

If the code changes run the following commands:
```bash
cd mof-guest-toolkit
conda activate mof-toolkit

git pull
conda env update -f environment.yml --prune
pip install -e .
```
---

## Usage examples

All examples below use the small demo file at `tests/data/molecules_example.csv`:

```
CID,Name,Abbreviation,Guest_Type
3033,Diclofenac,DIC,Pharmaceutical
3672,Ibuprofen,IBU,Pharmaceutical
2554,Carbamazepine,CAR,Pharmaceutical
30219,Cannabichromene,CBC,Cannabinoid
3084339,Cannabichromenic acid,CBCA,Cannabinoid
644019,Cannabidiol,CBD,Cannabinoid
```

---

### Example 1 — interactive 3D viewer

Launches a local web app. Open the printed URL in your browser
(Windows browser for WSL users). Type a CID or name and press Enter.

```bash
pubchem_interactive
```

```
  PubChem viewer running at: http://localhost:5050
  WSL users: open that URL in your Windows browser.
  Press Ctrl+C to stop.
```

The viewer shows the 3D structure (rotatable, zoomable) and a descriptor table.

---

### Example 2 — check a single RDKit property

```bash
pubchem_check_prop <cid_or_name> <rdkit_function>
```

```bash
pubchem_check_prop 3033 rdMolDescriptors.CalcTPSA
```
```
  Compound : 3033  (CID 3033)
  Function : rdMolDescriptors.CalcTPSA
  Value    : 49.33...
```

```bash
pubchem_check_prop cannabidiol rdMolDescriptors.CalcNumAromaticRings
```
```
  Compound : cannabidiol  (CID 644019)
  Function : rdMolDescriptors.CalcNumAromaticRings
  Value    : 1
```

```bash
pubchem_check_prop 3033 Fragments.fr_Ar_OH
```
```
  Compound : 3033  (CID 3033)
  Function : Fragments.fr_Ar_OH
  Value    : 0
```

Supported namespaces: `Chem`, `rdMolDescriptors`, `Fragments`, `Descriptors`, `GraphDescriptors`.
If the function requires extra arguments beyond `mol`, a clear error message is printed.

---

### Example 3 — fetch properties for one compound in Python

`get_xyz_cid` is a Python helper (not a CLI command) that returns a property
dict for one compound. Use it in scripts and notebooks to build results lists.

```python
from mof_toolkit.pubchem import get_xyz_cid, display_table

# Single compound — by CID or name
props = get_xyz_cid(3033)
props = get_xyz_cid("aspirin")
props = get_xyz_cid("cannabidiol")

# Print as a table
display_table(props)
```

```
Property               Value
----------------------------------------------------
CID                    3033
Name                   2-(2,6-dichloroanilino)...
MolecularWeight        295.0185
HBA                    2
HBD                    2
RotatableBonds         3
TPSA                   49.3300
...
```

---

### Example 4 — loop over a list of CIDs

```python
from mof_toolkit.pubchem import get_xyz_cid, display_table
import csv

cids = [3033, 3672, 644019]
results = []

for cid in cids:
    props = get_xyz_cid(cid)
    if props:
        results.append(props)

# Print each as a table
for r in results:
    display_table(r)
    print()
```

**Save the results list to CSV:**

```python
import csv

fieldnames = list(results[0].keys())
with open("my_results.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)
```

> For large lists, use `pubchem_batch_fetcher` directly, it handles errors, skips bad CIDs gracefully, and is faster for many compounds.

---

### Example 5 — batch compute descriptors

The easiest way to process a list of compounds:

```bash
pubchem_batch_fetcher tests/data/molecules_example.csv results.csv
```

```
  [1/6] CID 3033 ... done
  [2/6] CID 3672 ... done
  [3/6] CID 2554 ... done
  [4/6] CID 30219 ... done
  [5/6] CID 3084339 ... done
  [6/6] CID 644019 ... done

Saved 6/6 rows → results.csv
```

The output carries through all metadata columns and appends all default descriptors:

```
CID,Name,Abbreviation,Guest_Type,MolecularWeight,NumRings,...,fr_sulfonamd
3033,Diclofenac,DIC,Pharmaceutical,295.0185,2,...
```

**With a custom extra property** — add a `PropName::RDKitFunction` column header:

```
CID,Name,Guest_Type,Chi0v::Chem.rdMolDescriptors.CalcChi0v
3033,Diclofenac,Pharmaceutical
644019,Cannabidiol,Cannabinoid
```

```bash
pubchem_batch_fetcher my_list.csv results.csv
```

---

### Example 6 — download XYZ files

```bash
fetch_xyz_batch tests/data/molecules_example.csv ./xyz_files/
```

```
  [1/6] DIC (CID 3033)
    Saved: DIC.xyz
  [2/6] IBU (CID 3672)
    Saved: IBU.xyz
  ...
```

Each `.xyz` file follows the standard format:
```
30
DIC  CID=3033  source: PubChem 3D conformer
C      1.234567    -0.123456     0.000000
...
```

---

## Command reference

| Command | Description |
|---|---|
| `pubchem_interactive` | Local web viewer — 3D structure + descriptor table in browser |
| `pubchem_check_prop <cid> <func>` | Evaluate one RDKit function on a compound |
| `pubchem_batch_fetcher <in.csv> <out.csv>` | Batch-compute descriptors for a list |
| `fetch_xyz_batch <in.csv> <out_dir>` | Batch-download 3D XYZ files |

**Python-only helper** (import in scripts/notebooks, not available as CLI):

| Function | Description |
|---|---|
| `get_xyz_cid(cid_or_name)` | Returns property dict for one compound |

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
├── tests/
│   └── data/
│       └── molecules_example.csv   ← 6-compound demo file
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
