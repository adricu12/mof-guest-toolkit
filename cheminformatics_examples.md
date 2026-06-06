# Cheminformatics examples

## Index

1. [Interactive 3D viewer](#example-1--interactive-3d-viewer)
2. [Resolve a compound identifier](#example-2--resolve-a-compound-identifier)
3. [Compute descriptors for one compound](#example-3--compute-descriptors-for-one-compound)
4. [Batch descriptor computation](#example-4--batch-descriptor-computation)
5. [Generate 3D structure files](#example-5--generate-3d-structure-files)
6. [Convert between structure formats](#example-6--convert-between-structure-formats)
7. [Get SMILES](#example-7--get-smiles)

---

### Example 1 — interactive 3D viewer

Launch the browser-based Molecule Explorer to visualise any compound in 3D
and inspect its descriptors.

```bash
interactive_explorer
```

```
  Molecule Explorer running at: http://localhost:5050
  WSL users: open that URL in your Windows browser.
  Press Ctrl+C to stop.
```

The viewer includes:

- **Input type selector** — choose CID, Name, Formula, SMILES, or InChIKey.
- **Rotatable/zoomable 3D structure** rendered with 3Dmol.js.
- **Descriptor table** showing the [default](#default-descriptor-set) set.
- **EBI Unichem cross-references** for every compound that has an InChIKey,
  lists which databases (ChEMBL, DrugBank, KEGG, HMDB, …) contain the compound.
- **Surprise me 🎲** — picks a random valid PubChem compound (biased toward
  common, well-characterised compounds).
- **Multi-match warning** — when a name or formula maps to several PubChem
  entries, the viewer shows the first result and lists the other CIDs.

---

### Example 2 — resolve a compound identifier

`mol_explorer -input` resolves any supported identifier and prints the
canonical CID, SMILES, name, molecular formula, and InChIKey.

**Supported input types** (pass with `-format`; default: auto-detect):

| Flag value | Accepts |
|---|---|
| `cid`      | PubChem CID integer (e.g. `3033`) |
| `name`     | IUPAC or common name (e.g. `aspirin`, `cannabidiol`) |
| `formula`  | Molecular formula (e.g. `C9H8O4`, `c6h6`, `NaHCO3`) |
| `smiles`   | SMILES string (e.g. `CC(=O)Oc1ccccc1C(=O)O`) |
| `inchikey` | 27-character InChIKey (e.g. `BSYNRYMUTXBXSQ-UHFFFAOYSA-N`) |

```bash
mol_explorer -input 3033
mol_explorer -input aspirin          -format name
mol_explorer -input C9H8O4           -format formula
mol_explorer -input "CC(=O)Oc1ccccc1C(=O)O" -format smiles
mol_explorer -input BSYNRYMUTXBXSQ-UHFFFAOYSA-N -format inchikey
```

Example output for `mol_explorer -input 3033`:

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

**Multi-match warning** — when a formula or name maps to several PubChem entries,
the first (lowest CID) is used and a warning is printed:

```
Warning: 'C6H6' matched 14 PubChem entries. Using CID 241. Other CIDs: 8092, 12282, …
```

#### **Python API** version — use **`resolve_compound_input(query, input_type)`** to produce full resolution with all metadata:

```python
from mof_toolkit.rdkit_descriptors import resolve_compound_input

res = resolve_compound_input("aspirin")
print(res["cid"])         # 2244
print(res["smiles"])      # 'CC(=O)Oc1ccccc1C(=O)O'  (canonical)
print(res["formula"])     # 'C9H8O4'
print(res["inchikey"])    # 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
print(res["all_cids"])    # [2244]  — only one match for 'aspirin'

# Formula → multiple matches
res = resolve_compound_input("C6H6", "formula")
print(res["cid"])         # 241  (benzene — lowest CID = most canonical)
print(res["all_cids"])    # [241, 8092, 12282, …]
```

---

### Example 3 — compute descriptors for one compound

**[Default descriptor set](#default-descriptor-set)** (printed when no `-descriptor` flag is given):

```bash
mol_explorer -input aspirin -descriptor default
mol_explorer -input 3033    -descriptor default -output diclofenac.csv
```

**Specific RDKit descriptors** — use the exact names from [`Descriptors.descList`](#full-descriptor-set):

```bash
mol_explorer -input aspirin -descriptor TPSA MolLogP NumHDonors
mol_explorer -input "CC(=O)O" -descriptor TPSA MolLogP -output acetic_acid.csv
```

**Full ~210 descriptor set**:

```bash
mol_explorer -input aspirin    -descriptor full
mol_explorer -input cannabidiol -descriptor full -output cbd_all.csv
```

#### **Python API `get_rdkit_dict`** computes descriptor from any of the supporter input type and **`display_table`** display a table of the computed descriptors from any of the supported input_type. 

```python
from mof_toolkit.rdkit_descriptors import get_rdkit_dict, display_table

# Default set
display_table(get_rdkit_dict(3033))
display_table(get_rdkit_dict("aspirin"))
display_table(get_rdkit_dict("CC(=O)Oc1ccccc1C(=O)O"))  # SMILES

# Specific descriptors
display_table(get_rdkit_dict("aspirin", descriptors=["TPSA", "MolLogP"]))

# Full set
display_table(get_rdkit_dict("aspirin", full=True))

# SMILES not in PubChem — CID/names are empty; all descriptors still computed
display_table(get_rdkit_dict("C1CC1"))
```

#### When dealing directly with **SMILES, Python API `get_descriptor(smiles, descriptor_names)` and `get_all_rdkit_descriptors(mol)`** computes any RDKit descriptors:

```python
from mof_toolkit.rdkit_descriptors import get_descriptor

# Single descriptor
get_descriptor("CC(=O)O", "TPSA")
# → {'TPSA': 37.3}

# Multiple descriptors
get_descriptor("CC(=O)O", ["TPSA", "MolLogP", "NumHDonors"])
# → {'TPSA': 37.3, 'MolLogP': -0.17, 'NumHDonors': 1}

# The full ~210-name list
from mof_toolkit.rdkit_descriptors import get_all_rdkit_descriptors
from rdkit import Chem
mol = Chem.MolFromSmiles("CC(=O)O")
all_desc = get_all_rdkit_descriptors(mol)
print(sorted(all_desc.keys()))   # prints all ~210 available names
```


---

### Example 4 — batch descriptor computation

This tool process a CSV input with a long list of (guest) molecules, with **at least one** of the columns:
`CID`, `Name`, `SMILES`, `InChIKey`.
All four may be present simultaneously.

**Resolution priority (highest first): SMILES > CID > InChIKey > Name**

- Each field is validated independently; invalid values are flagged and ignored.
- If a SMILES column is present it is authoritative — other identifier columns
  are updated from the SMILES lookup and a note is printed.
- If two valid fields point to different PubChem compounds, a `MISMATCH` warning
  is printed and the higher-priority field is used.
- When a query maps to multiple PubChem entries the first (lowest CID) is used
  and a warning lists up to 10 other CIDs; more entries are noted with "… and more!".

**Example input CSV** ([`molecules_example.csv`](./tests/data/molecules_example.csv)):

```
CID,Name,Abbreviation,Guest_Type
3033,Diclofenac,DIC,Pharmaceutical
3672,Ibuprofen,IBU,Pharmaceutical
2554,Carbamazepine,CAR,Pharmaceutical
30219,Cannabichromene,CBC,Cannabinoid
3084339,Cannabichromenic acid,CBCA,Cannabinoid
644019,Cannabidiol,CBD,Cannabinoid
```

**[Default descriptors:](#default-descriptor-set)**

```bash
mol_explorer -batch ./tests/data/molecules_example.csv -output results.csv
```

**[Full ~210 descriptors:](#full-descriptor-set)**

```bash
mol_explorer -batch ./tests/data/molecules_example.csv -descriptor full -output results_full.csv
```

**Specific descriptors:**

```bash
mol_explorer -batch ./tests/data/molecules_example.csv -descriptor TPSA MolLogP NumHDonors -output results_selected.csv
```

Resolved columns appended to every output row:
`CID_resolved`, `IUPAC_Name_resolved`, `Common_Name_resolved`, `SMILES_resolved`,
`Formula_resolved`, `InChIKey_resolved`.

For small set of molecules, a **Python loop** can be used as follows:

```python
from mof_toolkit.rdkit_descriptors import get_rdkit_dict
import csv

queries = [3033, "ibuprofen", "CC(=O)Oc1ccccc1C(=O)O"]
results = [r for q in queries if (r := get_rdkit_dict(q)) is not None]

with open("my_results.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
    writer.writeheader()
    writer.writerows(results)
```

---

### Example 5 — generate 3D structure files

`mol_get_xyz` generates 3D coordinate files from any supported identifier.
PubChem 3D conformers are used when a CID is available; otherwise RDKit
ETKDGv3 + MMFF94 is used locally (no internet required).

**Single compound:**

```bash
mol_get_xyz -input 3033
mol_get_xyz -input aspirin -outputformat xyz sdf
mol_get_xyz -input "CC(=O)Oc1ccccc1C(=O)O" -outputformat sdf pdb -output aspirin
mol_get_xyz -input DCOPUUMXTXDBNB-UHFFFAOYSA-N -outputformat xyz  # InChIKey
```

**Force local RDKit generation** (ignore PubChem 3D conformer):

```bash
mol_get_xyz -input 3033 -outputformat xyz --source rdkit
```

**Batch from CSV:**

```bash
mol_get_xyz -batch example.csv -output ./structures/
mol_get_xyz -batch example.csv -outputformat xyz sdf -output ./structures/

# example using the csv file
mol_get_xyz -batch ./tests/data/molecules_example.csv -outputformat xyz -output ./structures/
```

**Custom output naming with `-namecol`** — use any column(s) from the CSV as the file stem.
One or more column names can be given (case-insensitive); multiple values are joined with `_`.

```bash
# Use the Name column → Diclofenac.xyz, Ibuprofen.xyz, …
mol_get_xyz -batch ./tests/data/molecules_example.csv -outputformat xyz -output ./structures/ -namecol name

# Combine CID and Name → 3033_Diclofenac.xyz, 3672_Ibuprofen.xyz, …
mol_get_xyz -batch ./tests/data/molecules_example.csv -outputformat xyz -output ./structures/ -namecol cid name

# Use the Abbreviation column → DIC.xyz, IBU.xyz, …
mol_get_xyz -batch ./tests/data/molecules_example.csv -outputformat xyz -output ./structures/ -namecol abbreviation
```

When `-namecol` is omitted the default auto-naming (CID + PubChem resolved name) is used.

**Auto-generated filename rules** (default, when `-namecol` is not specified):

| Available info      | Filename example        |
|---------------------|-------------------------|
| CID + common name   | `3033_Diclofenac.xyz`   |
| CID only            | `3033.xyz`              |
| Common name only    | `Diclofenac.xyz`        |
| SMILES only (no PubChem match) | molecular formula, e.g. `C3H6.xyz` |

#### **Python API `get_3d_structure()`** — generate or download a 3D structure file (requires internet for PubChem lookup).

```python
from mof_toolkit.molecule_manager import get_3d_structure

# By CID — PubChem 3D conformer preferred
get_3d_structure(3033, formats=["xyz", "sdf"], output_dir="./structures/")

# By name
get_3d_structure("aspirin", formats=["xyz"], output_dir="./structures/")

# By SMILES — RDKit always used (no CID)
get_3d_structure("CC(=O)Oc1ccccc1C(=O)O", formats=["sdf", "pdb"],
                 output_stem="./structures/aspirin")

# Force RDKit
get_3d_structure(3033, formats=["xyz"], output_dir="./out/", source="rdkit")

# Loop
for cid in [3033, 3672, 644019]:
    get_3d_structure(cid, formats=["xyz"], output_dir="./structures/")
```

**Supported output formats:** `xyz`, `sdf`, `mol`, `pdb`

#### **Python API `embed_3d(mol)`** — generate a 3D conformer locally from an RDKit Mol object (no internet required).

```python
from rdkit import Chem
from mof_toolkit.molecule_manager import embed_3d

mol_2d = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
mol_3d = embed_3d(mol_2d)   # ETKDGv3 + MMFF94; Hs added automatically
print(mol_3d.GetNumConformers())  # 1
```

---

### Example 6 — convert between structure formats

`mol_file_translate` converts any supported 3D structure file to another format.
The input and output formats are inferred from the file extensions.

```bash
mol_file_translate -input aspirin.sdf  -output aspirin.xyz
mol_file_translate -input aspirin.xyz  -output aspirin.pdb
mol_file_translate -input molecule.pdb -output molecule.mol
mol_file_translate -input compound.mol -output compound.sdf
```
---
### Example 7 — get SMILES  

#### **Python API `get_smiles(query, input_type)`** resolves any identifier to canonical SMILES

```python
from mof_toolkit.rdkit_descriptors import get_smiles

get_smiles(3033)                        # CID   → 'OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl'
get_smiles("aspirin")                   # name  → 'CC(=O)Oc1ccccc1C(=O)O'
get_smiles("C9H8O4", "formula")         # formula → 'CC(=O)Oc1ccccc1C(=O)O'
get_smiles("DCOPUUMXTXDBNB-UHFFFAOYSA-N", "inchikey")  # InChIKey
```

#### **Python API `get_smiles_from_coords`** reads SMILES from a structure file:

```python
from mof_toolkit.molecule_manager import get_smiles_from_coords

smiles = get_smiles_from_coords("aspirin.sdf")   # → 'CC(=O)Oc1ccccc1C(=O)O'
smiles = get_smiles_from_coords("molecule.xyz")  # → canonical SMILES (bonds perceived from 3D coords)
smiles = get_smiles_from_coords("compound.pdb")  # → canonical SMILES
```

For XYZ files, bonds and bond orders are reconstructed from the 3D coordinates
using a distance + valence algorithm (adapted from [xyz2mol](https://github.com/jensengroup/xyz2mol)):

1. **Connectivity**: two atoms are bonded if their distance falls within
   `(r_cov_i + r_cov_j) × 1.3 Å` (covalent radii from [Alvarez 2008](https://doi.org/10.1039/B801115J)).
2. **Bond orders**: starting from all-single bonds, pairs are iteratively
   upgraded to double/triple wherever both atoms still have unsatisfied valence,
   greedy by highest combined remaining valence first.
3. **Aromaticity**: RDKit's sanitiser perceives it automatically from the
   resulting Kekulé form; no special ring handling is needed.

RDKit's built-in `rdDetermineBonds` (RDKit ≥ 2022.03)
is tried first and is the preferred path.


---

## Default descriptor set

| Key | RDKit function | Description |
|---|---|---|
| `MolecularFormula`  | `CalcMolFormula`          | Molecular formula string |
| `ExactMolWt`        | `CalcExactMolWt`          | Exact molecular weight (Da) |
| `NumHeavyAtoms`     | `mol.GetNumHeavyAtoms()`  | Non-hydrogen atom count |
| `RingCount`         | `CalcNumRings`            | Total ring count |
| `NumAromaticRings`  | `CalcNumAromaticRings`    | Aromatic ring count |
| `NumHBA`            | `CalcNumHBA`              | Hydrogen bond acceptors |
| `NumHBD`            | `CalcNumHBD`              | Hydrogen bond donors |
| `NumRotatableBonds` | `CalcNumRotatableBonds`   | Rotatable bond count |
| `TPSA`              | `CalcTPSA`                | Topological polar surface area (Å²) |

All outputs also include `CID`, `IUPAC_Name`, `Common_Name`, `SMILES`,
`Formula`, and `InChIKey`.

## Full descriptor set
To see all ~210 available descriptor names:

```python
from rdkit.Chem import Descriptors
print([name for name, _ in Descriptors.descList])
```