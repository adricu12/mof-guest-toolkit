# mof-guest-toolkit

A Python package of helper utilities for computational chemistry workflows.
The package provides command-line tools and importable functions for:

- [Quick interactive exploration of molecules](#example-1--interactive-3d-viewer)
- [Compute specific rdkit descriptor for a molecule]()
- [Getting descriptors for a set of molecules]()
- [Getting coord files of molecules]()

- Parsing geometry optimisation outputs from AMS/ORCA *(coming soon)*
- Analysing normal modes and gradient convergence *(coming soon)*
- Analysing LAMMPS trajectory files *(coming soon)*
- Structure checker *(coming soon)*


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

You need to activate the conda environment every time you open a new terminal.
To deactivate: `conda deactivate`.

---

### Verify the installation

```bash
rdkit_check_prop 3033 TPSA
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
(Windows browser for WSL users). Type a CID, name, or SMILES and press Enter.

```bash
pubchem_interactive
```

```
  PubChem viewer running at: http://localhost:5050
  WSL users: open that URL in your Windows browser.
  Press Ctrl+C to stop.
```

The viewer shows the rotatable/zoomable 3D structure and a descriptor table.

---

### Example 2 — check a single RDKit property

```bash
rdkit_check_prop <cid|name|smiles> <property>
```

The compound can be given as a CID, a common/IUPAC name, or a SMILES string.
The property can be a `DEFAULT_PROPERTIES` key (e.g. `TPSA`, `HBA`) or any
dotted RDKit callable.

```bash
rdkit_check_prop 3033 TPSA
```
```
CID         : 3033
IUPAC_Name  : 2-(2,6-dichloroanilino)phenylacetic acid
Common_Name : Diclofenac
SMILES      : OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl
Function    : TPSA
Value       : 49.330000
```

```bash
rdkit_check_prop cannabidiol NumAromaticRings
rdkit_check_prop 3033 Fragments.fr_Ar_OH
rdkit_check_prop "CC(=O)O" rdMolDescriptors.CalcTPSA
```

Supported namespaces: `Chem`, `rdMolDescriptors`, `Fragments`, `Descriptors`, `GraphDescriptors`.

Run `rdkit_check_prop --help` for the full list of DEFAULT_PROPERTIES keys.

---

### Example 3 — print all default descriptors for one compound

```bash
rdkit_default_props <cid|name|smiles>
```

```bash
rdkit_default_props 2244
rdkit_default_props aspirin
rdkit_default_props "CC(=O)Oc1ccccc1C(=O)O"
```

Prints CID, IUPAC name, common name, SMILES, and all `DEFAULT_PROPERTIES` as
an aligned table.  If the SMILES is valid but not in PubChem, identifiers are
shown as N/A but all descriptors are still computed.

---

### Example 4 — fetch properties for one compound in Python

`get_rdkit_dict` accepts a CID, name, or SMILES and returns a property dict.
SMILES-only compounds (not in PubChem) return empty CID/name fields but full descriptors.

```python
from mof_toolkit.rdkit_properties import get_rdkit_dict, display_table

display_table(get_rdkit_dict(3033))
display_table(get_rdkit_dict("aspirin"))
display_table(get_rdkit_dict("CC(=O)Oc1ccccc1C(=O)O"))   # SMILES input

# SMILES not in PubChem — CID/names empty, all descriptors still computed
display_table(get_rdkit_dict("C1CC1"))
```

```
Property                 Value
--------------------------------------------------------
CID                      3033
IUPAC_Name               2-(2,6-dichloroanilino)...
Common_Name              Diclofenac
SMILES                   OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl
MolecularWeight          295.0185
HBA                      2
...
```

---

### Example 5 — loop over a list of compounds

```python
from mof_toolkit.rdkit_properties import get_rdkit_dict, display_table
import csv

queries = [3033, "ibuprofen", "CC(=O)Oc1ccccc1C(=O)O"]
results = [get_rdkit_dict(q) for q in queries if get_rdkit_dict(q)]

# Save to CSV
fieldnames = list(results[0].keys())
with open("my_results.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)
```

---

### Example 6 — batch compute descriptors

```bash
rdkit_batch_fetcher tests/data/molecules_example.csv results.csv
```

The input CSV must have **at least one** of `CID`, `Name`, or `SMILES` columns.
All three may be present simultaneously.

**Resolution rules (highest priority first): SMILES > CID > Name**

- Each field is validated independently. Invalid fields are flagged and skipped.
- If two valid fields point to different PubChem compounds, a `MISMATCH` flag is printed
  and the higher-priority input is used.
- The resolved identifiers (CID, IUPAC name, common name, SMILES) are always appended
  as `CID_resolved`, `Name_IUPAC_resolved`, `Name_common_resolved`, `SMILES_resolved`.

Example output:

```
  [1/6] Using CID as primary input (CID=3033)
    → done
  [2/6] Using CID as primary input (CID=3672)
    → done
  ...

Saved 6/6 rows → results.csv
```

**With custom extra properties** — add `PropName::RDKitFunction` column headers:

```
CID,Name,Guest_Type,Chi0v::Chem.rdMolDescriptors.CalcChi0v
3033,Diclofenac,Pharmaceutical,
644019,Cannabidiol,Cannabinoid,
```

```bash
rdkit_batch_fetcher my_list.csv results.csv
```

**Mixed-input example** — each row can use a different identifier:

```
CID,Name,SMILES,Guest_Type
3033,,,Pharmaceutical
,Ibuprofen,,Pharmaceutical
,,CC(=O)Oc1ccccc1C(=O)O,Other
3033,Aspirin,,Mismatch_example
```

Row 4 triggers a `MISMATCH` flag (CID 3033 = Diclofenac ≠ Aspirin) and CID takes priority.

Run `rdkit_batch_fetcher --help` for full usage.

---

### Example 7 — generate 3D structure files (single compound)

`get_3d_structure` is the recommended Python helper for single-compound use.
It first tries to download a PubChem 3D conformer; if unavailable it falls back
to RDKit ETKDGv3 + MMFF94.

```python
from mof_toolkit.molecule_manager import get_3d_structure

# By CID — downloads PubChem conformer when available
get_3d_structure(3033, formats=["xyz", "sdf"], output_dir="./structures/")

# By name
get_3d_structure("aspirin", formats=["xyz"], output_dir="./structures/")

# By SMILES — always uses RDKit (no CID known)
get_3d_structure("CC(=O)Oc1ccccc1C(=O)O", formats=["sdf", "pdb"],
                 output_stem="./structures/aspirin")

# Force RDKit even when a CID is available
get_3d_structure(3033, formats=["xyz"], output_dir="./out/", source="rdkit")
```

**Output file naming** (when `output_stem` is not given):

| Available info | Filename |
|---|---|
| CID + common name | `3033_Diclofenac.xyz` |
| CID only | `3033.xyz` |
| common name only | `Diclofenac.xyz` |
| SMILES only, not in PubChem | molecular formula, e.g. `C3H6.xyz` |

**Loop over a list of CIDs:**

```python
cids = [3033, 3672, 644019]
for cid in cids:
    get_3d_structure(cid, formats=["xyz"], output_dir="./structures/")
```

**Loop over a mixed list (CIDs, names, SMILES):**

```python
queries = [3033, "ibuprofen", "CC(=O)Oc1ccccc1C(=O)O"]
for q in queries:
    get_3d_structure(q, formats=["xyz", "sdf"], output_dir="./structures/")
```

**Loop with custom output names:**

```python
entries = [("CC(=O)O", "acetic_acid"), ("c1ccccc1", "benzene")]
for smiles, name in entries:
    get_3d_structure(smiles, formats=["xyz", "sdf"],
                     output_stem=f"./structures/{name}")
```

---

### Example 8 — batch 3D structure files from CSV

```bash
fetch_xyz_batch tests/data/molecules_example.csv ./structures/
fetch_xyz_batch tests/data/molecules_example.csv ./structures/ --format xyz sdf
fetch_xyz_batch tests/data/molecules_example.csv ./structures/ --format xyz pdb --code-col Abbreviation
```

The input CSV uses the **same format** as `rdkit_batch_fetcher` (at least one of `CID`, `Name`,
`SMILES`; same resolution rules apply).

**File naming priority:**
1. `--code-col` value (e.g. `Abbreviation` column → `DIC.xyz`)
2. `CID_CommonName` (e.g. `3033_Diclofenac.xyz`)
3. `CID` only (e.g. `3033.xyz`)
4. Common name only (e.g. `Diclofenac.xyz`)
5. `missing_id01.xyz`, `missing_id02.xyz`, ... (SMILES valid but not in PubChem and no name)

Each compound is saved as `<stem>.<format>` in the output directory.

Run `fetch_xyz_batch --help` for full usage.

---

### Example 9 — generate 3D structure from SMILES (CLI, no internet)

```bash
smiles_to_3d "CC(=O)Oc1ccccc1C(=O)O"
smiles_to_3d "CC(=O)Oc1ccccc1C(=O)O" --format xyz sdf pdb mol
smiles_to_3d "CC(=O)Oc1ccccc1C(=O)O" --output aspirin --format sdf pdb
```

3D coordinates are generated **locally** with RDKit (ETKDGv3 + MMFF94).
No internet connection needed. Supported formats: `xyz`, `sdf`, `pdb`, `mol`.

Run `smiles_to_3d --help` for full usage.

---

## Command reference

| Command | Description |
|---|---|
| `pubchem_interactive` | Local web viewer — 3D structure + descriptor table in browser |
| `rdkit_check_prop <cid|name|smiles> <prop>` | Evaluate one RDKit property on a compound |
| `rdkit_default_props <cid|name|smiles>` | Print all default descriptors for a compound |
| `rdkit_batch_fetcher <in.csv> <out.csv>` | Batch-compute descriptors from CSV (CID/Name/SMILES) |
| `fetch_xyz_batch <in.csv> <out_dir>` | Batch-generate 3D structure files from CSV |
| `smiles_to_3d <SMILES>` | Generate 3D structure file(s) from a SMILES string (local) |

All commands accept `--help` for full usage details.

**Python-only helpers** (import in scripts/notebooks):

| Function | Module | Description |
|---|---|---|
| `get_rdkit_dict(query)` | `rdkit_properties` | Property dict for one compound (CID/name/SMILES) |
| `resolve_compound_input(query)` | `rdkit_properties` | Resolves CID/name/SMILES → unified dict |
| `get_3d_structure(query, ...)` | `molecule_manager` | Generate 3D file(s) for one compound |
| `embed_3d(mol)` | `molecule_manager` | RDKit ETKDGv3 + MMFF94 conformer for an RDKit Mol |

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

All outputs also include `CID`, `IUPAC_Name`, `Common_Name`, and `SMILES`.

---

## Project structure

```
mof-guest-toolkit/
├── mof_toolkit/
│   ├── __init__.py
│   ├── rdkit_properties.py  ← PubChem fetchers, RDKit descriptors, batch CLI tools
│   ├── molecule_manager.py  ← 3D conformer generation and structure file writing
│   ├── ccdc.py              ← CIF fetcher by CCDC refcode        [coming soon]
│   ├── geometry_opt.py      ← AMS/ORCA .out parser               [coming soon]
│   ├── normal_modes.py      ← frequency and imaginary mode tools  [coming soon]
│   └── bonding.py           ← connectivity and bond distance checks [coming soon]
├── tests/
│   └── data/
│       └── molecules_example.csv   ← 6-compound demo file
├── pyproject.toml           ← package metadata and CLI entry points
├── environment.yml          ← conda environment definition
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
