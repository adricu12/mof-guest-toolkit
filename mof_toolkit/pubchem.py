"""
pubchem.py
----------
Fetch chemical information and compute molecular descriptors
from PubChem CIDs or compound names.

CLI commands
------------
  pubchem_interactive              — open an interactive Jupyter notebook
  pubchem_check <cid> <func>       — evaluate one RDKit function on a compound
  pubchem_batch_fetcher <in> <out> — batch-compute descriptors from a CSV
  get_xyz_cid <cid_or_name>        — print descriptor table in the terminal
"""

import csv
import json
import os
import subprocess
import tempfile

import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import Fragments, rdMolDescriptors


# ---------------------------------------------------------------------------
# Descriptor registry
# ---------------------------------------------------------------------------

DEFAULT_PROPERTIES = {
    "MolecularWeight":  rdMolDescriptors.CalcExactMolWt,
    "NumRings":         rdMolDescriptors.CalcNumRings,
    "NumAromaticRings": rdMolDescriptors.CalcNumAromaticRings,
    "HBA":              rdMolDescriptors.CalcNumHBA,
    "HBD":              rdMolDescriptors.CalcNumHBD,
    "RotatableBonds":   rdMolDescriptors.CalcNumRotatableBonds,
    "TPSA":             rdMolDescriptors.CalcTPSA,
    "fr_Al_OH":         Fragments.fr_Al_OH,
    "fr_Ar_OH":         Fragments.fr_Ar_OH,
    "fr_COO":           Fragments.fr_COO,
    "fr_C_O_noCOO":     Fragments.fr_C_O_noCOO,
    "fr_Ar_N":          Fragments.fr_Ar_N,
    "fr_NH2":           Fragments.fr_NH2,
    "fr_NH1":           Fragments.fr_NH1,
    "fr_NH0":           Fragments.fr_NH0,
    "fr_ether":         Fragments.fr_ether,
    "fr_sulfonamd":     Fragments.fr_sulfonamd,
}


# ---------------------------------------------------------------------------
# Core fetchers
# ---------------------------------------------------------------------------

def fetch_cid_from_name(name: str) -> int:
    """Resolve a compound name to a PubChem CID."""
    compounds = pcp.get_compounds(name, "name")
    if not compounds:
        raise ValueError(f"No compound found for name: '{name}'")
    return compounds[0].cid


def fetch_smiles_from_cid(cid: int) -> str:
    """Fetch canonical SMILES for a given PubChem CID."""
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
        "/property/CanonicalSMILES/TXT"
    )
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    return response.text.strip()


def fetch_pubchem_metadata(cid: int) -> dict:
    """Fetch name, formula, MW, and SMILES for a CID from the PubChem REST API."""
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
        "/property/IUPACName,MolecularWeight,MolecularFormula,CanonicalSMILES/JSON"
    )
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    p = response.json()["PropertyTable"]["Properties"][0]
    return {
        "CID": cid,
        "Name": p.get("IUPACName", ""),
        "MolecularFormula": p.get("MolecularFormula", ""),
        "MolecularWeight_API": p.get("MolecularWeight", ""),
        "SMILES": p.get("CanonicalSMILES", ""),
    }


def smiles_to_mol(smiles: str) -> Chem.Mol:
    """Convert SMILES string to RDKit Mol, raising on failure."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: '{smiles}'")
    return mol


def resolve_to_cid(input_strg: str) -> int:
    """Accept a numeric CID string or a compound name and return the CID."""
    if input_strg.strip().isdigit():
        return int(input_strg.strip())
    return fetch_cid_from_name(input_strg)


# ---------------------------------------------------------------------------
# RDKit function resolver (used by pubchem_check and batch_fetcher)
# ---------------------------------------------------------------------------

def _resolve_rdkit_func(func_str: str):
    """
    Resolve a dotted RDKit expression string to a callable via eval.

    Supported namespaces: Chem, rdMolDescriptors, Fragments, Descriptors,
    GraphDescriptors.

    Raises ValueError if the string cannot be resolved, TypeError if the
    result is not callable.
    """
    from rdkit.Chem import Descriptors, GraphDescriptors  # noqa: PLC0415

    eval_ctx = {
        "Chem": Chem,
        "rdMolDescriptors": rdMolDescriptors,
        "Fragments": Fragments,
        "Descriptors": Descriptors,
        "GraphDescriptors": GraphDescriptors,
    }
    try:
        func = eval(func_str, {"__builtins__": {}}, eval_ctx)  # noqa: S307
    except Exception as exc:
        raise ValueError(f"Cannot resolve '{func_str}': {exc}") from exc
    if not callable(func):
        raise TypeError(
            f"'{func_str}' evaluates to {type(func).__name__}, not a callable"
        )
    return func


# ---------------------------------------------------------------------------
# Descriptor computation
# ---------------------------------------------------------------------------

def compute_properties(cid: int, properties: dict = None) -> dict:
    """Compute RDKit descriptors for a CID. Defaults to DEFAULT_PROPERTIES."""
    if properties is None:
        properties = DEFAULT_PROPERTIES
    smiles = fetch_smiles_from_cid(cid)
    mol = smiles_to_mol(smiles)
    return {name: func(mol) for name, func in properties.items()}


# ---------------------------------------------------------------------------
# Terminal display
# ---------------------------------------------------------------------------

def display_table(property_dict: dict):
    """Print a property dict as an aligned table."""
    print(f"\n{'Property':<22} {'Value':<30}")
    print("-" * 52)
    for prop, value in property_dict.items():
        if isinstance(value, float):
            value = f"{value:.4f}"
        print(f"{prop:<22} {str(value):<30}")


def interact_with_pubchem(input_strg: str, properties: dict = None):
    """Resolve a CID/name, compute descriptors, and print a table."""
    try:
        cid = resolve_to_cid(input_strg)
        print(f"Resolved to CID: {cid}")
        result = compute_properties(cid, properties)
        display_table(result)
        return result
    except Exception as e:
        print(f"Error: {e}")
        return None


# ---------------------------------------------------------------------------
# Jupyter display helper (call this from inside a notebook)
# ---------------------------------------------------------------------------

def show_molecule(cid: int):
    """Render the 3D structure and a short property table inside a Jupyter cell."""
    import py3Dmol  # noqa: PLC0415
    import pandas as pd  # noqa: PLC0415
    from IPython.display import display  # noqa: PLC0415

    url_3d = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
        "/SDF?record_type=3d"
    )
    resp = requests.get(url_3d, timeout=10)
    if resp.status_code != 200:
        print(f"No 3D structure found for CID={cid} (HTTP {resp.status_code})")
    else:
        viewer = py3Dmol.view(width=500, height=350)
        viewer.addModel(resp.text, "sdf")
        viewer.setStyle({"stick": {}, "sphere": {"scale": 0.25}})
        viewer.zoomTo()
        viewer.show()

    meta = fetch_pubchem_metadata(cid)
    mol = smiles_to_mol(meta["SMILES"])
    props = {
        "CID": cid,
        "Name": meta["Name"],
        "Molecular Weight": f"{rdMolDescriptors.CalcExactMolWt(mol):.4f}",
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "Rotatable Bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
    }
    display(
        pd.DataFrame(props.items(), columns=["Property", "Value"]).set_index("Property")
    )


# ---------------------------------------------------------------------------
# XYZ file writer
# ---------------------------------------------------------------------------

def fetch_and_save_xyz(cid, abbreviation, output_dir):
    """Download the PubChem 3D conformer and write it as <abbreviation>.xyz."""
    compounds = pcp.get_compounds(cid, "cid", record_type="3d")
    if not compounds:
        print(f"  No 3D conformer for CID={cid} ({abbreviation}) — skipping")
        return
    atoms = compounds[0].atoms
    lines = [
        str(len(atoms)),
        f"{abbreviation}  CID={cid}  source: PubChem 3D conformer",
    ]
    for atom in atoms:
        lines.append(
            f"{atom.element:<3}  {atom.x:12.6f}  {atom.y:12.6f}  {atom.z:12.6f}"
        )
    filepath = os.path.join(output_dir, f"{abbreviation}.xyz")
    with open(filepath, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  Saved: {os.path.basename(filepath)}")


# ---------------------------------------------------------------------------
# Notebook cell source strings (self-contained; no dependency on mof_toolkit)
# ---------------------------------------------------------------------------

_NB_IMPORTS = """\
# PubChem Interactive Viewer
import requests
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import ipywidgets as widgets
from IPython.display import display, clear_output
import py3Dmol
import pandas as pd
print("Libraries loaded — enter a CID or name below and click Look up.")
"""

_NB_HELPERS = """\
def _resolve_cid(query):
    if query.strip().isdigit():
        return int(query.strip())
    hits = pcp.get_compounds(query, 'name')
    if not hits:
        raise ValueError(f'No compound found: {query!r}')
    return hits[0].cid

def _fetch_smiles(cid):
    r = requests.get(
        f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}'
        '/property/CanonicalSMILES/TXT',
        timeout=10)
    r.raise_for_status()
    return r.text.strip()

def _fetch_name(cid):
    r = requests.get(
        f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}'
        '/property/IUPACName/JSON',
        timeout=10)
    r.raise_for_status()
    return r.json()['PropertyTable']['Properties'][0].get('IUPACName', 'N/A')

def show_compound(query):
    cid = _resolve_cid(query)
    smiles = _fetch_smiles(cid)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'RDKit could not parse SMILES for CID {cid}')
    name = _fetch_name(cid)

    r3d = requests.get(
        f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d',
        timeout=10)
    if r3d.status_code == 200:
        v = py3Dmol.view(width=500, height=350)
        v.addModel(r3d.text, 'sdf')
        v.setStyle({'stick': {}, 'sphere': {'scale': 0.25}})
        v.zoomTo()
        v.show()
    else:
        print(f'No 3D structure available for CID {cid} (HTTP {r3d.status_code})')

    props = {
        'CID': cid,
        'Name': name,
        'Molecular Weight': f'{rdMolDescriptors.CalcExactMolWt(mol):.4f}',
        'HBA': rdMolDescriptors.CalcNumHBA(mol),
        'HBD': rdMolDescriptors.CalcNumHBD(mol),
        'Rotatable Bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
    }
    display(
        pd.DataFrame(props.items(), columns=['Property', 'Value']).set_index('Property')
    )
"""

_NB_WIDGET = """\
inp = widgets.Text(
    value='',
    placeholder='CID or molecule name  (e.g.  aspirin  or  2244)',
    layout=widgets.Layout(width='440px'),
)
btn = widgets.Button(description='Look up', button_style='primary')
out = widgets.Output()

def on_click(_):
    with out:
        clear_output(wait=True)
        q = inp.value.strip()
        if not q:
            print('Please enter a CID or molecule name.')
            return
        try:
            show_compound(q)
        except Exception as e:
            print(f'Error: {e}')

btn.on_click(on_click)
display(widgets.VBox([inp, btn, out]))
"""


def _nb_code_cell(source: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source,
    }


# ---------------------------------------------------------------------------
# CLI: pubchem_interactive
# ---------------------------------------------------------------------------

def pubchem_interactive_cli():
    """pubchem_interactive — generate and open an interactive Jupyter notebook."""
    notebook = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python"},
        },
        "cells": [
            _nb_code_cell(_NB_IMPORTS),
            _nb_code_cell(_NB_HELPERS),
            _nb_code_cell(_NB_WIDGET),
        ],
    }

    nb_path = os.path.join(tempfile.gettempdir(), "pubchem_interactive.ipynb")
    with open(nb_path, "w") as fh:
        json.dump(notebook, fh, indent=1)

    is_wsl = False
    try:
        with open("/proc/version") as fv:
            is_wsl = "microsoft" in fv.read().lower()
    except OSError:
        pass

    print(f"Notebook written to: {nb_path}")
    if is_wsl:
        print("WSL detected — Jupyter will not open a browser automatically.")
        print("Copy the URL printed below into your Windows browser.\n")
        subprocess.run(["jupyter", "notebook", nb_path, "--no-browser"], check=False)
    else:
        subprocess.run(["jupyter", "notebook", nb_path], check=False)


# ---------------------------------------------------------------------------
# CLI: pubchem_check
# ---------------------------------------------------------------------------

def pubchem_check_cli():
    """pubchem_check <cid_or_name> <rdkit_function_string>

    Evaluate a single RDKit function on a compound and print the result.

    Examples
    --------
      pubchem_check 3033 Chem.rdMolDescriptors.BCUT2D
      pubchem_check aspirin rdMolDescriptors.CalcTPSA
      pubchem_check 2244 Fragments.fr_COO
    """
    import argparse  # noqa: PLC0415

    parser = argparse.ArgumentParser(
        description="Evaluate one RDKit function on a PubChem compound",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Supported namespaces: Chem, rdMolDescriptors, Fragments, "
            "Descriptors, GraphDescriptors\n\n"
            "Examples:\n"
            "  pubchem_check 3033 Chem.rdMolDescriptors.BCUT2D\n"
            "  pubchem_check aspirin rdMolDescriptors.CalcTPSA"
        ),
    )
    parser.add_argument("compound", help="PubChem CID (integer) or compound name")
    parser.add_argument(
        "function",
        help="RDKit callable, e.g. Chem.rdMolDescriptors.BCUT2D",
    )
    args = parser.parse_args()

    try:
        cid = resolve_to_cid(args.compound)
    except Exception as e:
        print(f"Error resolving compound '{args.compound}': {e}")
        return

    try:
        smiles = fetch_smiles_from_cid(cid)
        mol = smiles_to_mol(smiles)
    except Exception as e:
        print(f"Error fetching molecule for CID {cid}: {e}")
        return

    try:
        func = _resolve_rdkit_func(args.function)
    except (ValueError, TypeError) as e:
        print(f"Error: {e}")
        print(
            "Hint: check the function path. Valid namespaces: "
            "Chem, rdMolDescriptors, Fragments, Descriptors, GraphDescriptors."
        )
        return

    try:
        result = func(mol)
        print(f"{args.function}  (CID {cid}):\n  {result}")
    except TypeError as e:
        print(f"Error calling {args.function}(mol): {e}")
        print(
            "Hint: this function may require extra arguments beyond mol "
            "(e.g. integer parameters). Check the RDKit docs."
        )


# ---------------------------------------------------------------------------
# CLI: pubchem_batch_fetcher
# ---------------------------------------------------------------------------

def pubchem_batch_fetcher_cli():
    """pubchem_batch_fetcher <input.csv> <output.csv>

    Input CSV format
    ----------------
    Required column  : CID
    Optional metadata: Name, Guest_type, or any other column — passed through.
                       If Name is blank, the IUPAC name is fetched from PubChem.
    Extra properties : column header format  PropName::rdkit.func.Path
                       e.g.  BCUT2D::Chem.rdMolDescriptors.BCUT2D

    Output CSV
    ----------
    All input metadata columns + all DEFAULT_PROPERTIES + any extra computed columns.
    """
    import argparse  # noqa: PLC0415

    parser = argparse.ArgumentParser(
        description="Batch-compute molecular descriptors from a CSV of PubChem CIDs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Extra property columns use '::' in the header:\n"
            "  BCUT2D::Chem.rdMolDescriptors.BCUT2D\n"
            "  NumHeteroatoms::rdMolDescriptors.CalcNumHeteroatoms"
        ),
    )
    parser.add_argument("input_csv", help="Input CSV file (must have a CID column)")
    parser.add_argument("output_csv", help="Output CSV file path")
    args = parser.parse_args()

    with open(args.input_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        headers = reader.fieldnames or []
        rows = list(reader)

    # Split headers into passthrough cols and extra computed cols (PropName::func)
    extra_props: dict[str, str] = {}
    passthrough_cols: list[str] = []
    for h in headers:
        if "::" in h:
            prop_name, func_str = h.split("::", 1)
            extra_props[prop_name.strip()] = func_str.strip()
        elif h.strip().upper() != "CID":
            passthrough_cols.append(h)

    # Resolve extra property functions once up front
    resolved_extra: dict[str, object] = {}
    for prop_name, func_str in extra_props.items():
        try:
            resolved_extra[prop_name] = _resolve_rdkit_func(func_str)
        except (ValueError, TypeError) as e:
            print(f"Warning: skipping extra property '{prop_name}': {e}")

    output_rows: list[dict] = []
    total = len(rows)
    for i, row in enumerate(rows, 1):
        cid_raw = (row.get("CID") or row.get("cid") or "").strip()
        if not cid_raw:
            print(f"  Row {i}/{total}: missing CID — skipping")
            continue

        print(f"  [{i}/{total}] {cid_raw} ...", end=" ", flush=True)
        try:
            cid = resolve_to_cid(cid_raw)
            smiles = fetch_smiles_from_cid(cid)
            mol = smiles_to_mol(smiles)
        except Exception as e:
            print(f"ERROR ({e})")
            continue

        out_row: dict = {"CID": cid}

        # Passthrough metadata; auto-fetch Name from PubChem when blank
        for col in passthrough_cols:
            val = row.get(col, "")
            if col.strip().lower() == "name" and not val.strip():
                try:
                    val = fetch_pubchem_metadata(cid).get("Name", "")
                except Exception:
                    val = ""
            out_row[col] = val

        # Default descriptor set
        for prop_name, func in DEFAULT_PROPERTIES.items():
            try:
                out_row[prop_name] = func(mol)
            except Exception as e:
                out_row[prop_name] = f"ERROR: {e}"

        # Extra user-specified descriptors
        for prop_name, func in resolved_extra.items():
            try:
                out_row[prop_name] = func(mol)
            except TypeError as e:
                out_row[prop_name] = f"ERROR (extra args needed): {e}"
            except Exception as e:
                out_row[prop_name] = f"ERROR: {e}"

        output_rows.append(out_row)
        print("done")

    if not output_rows:
        print("No rows processed — output file not written.")
        return

    fieldnames = list(output_rows[0].keys())
    with open(args.output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"\nSaved {len(output_rows)} rows → {args.output_csv}")


# ---------------------------------------------------------------------------
# CLI: get_xyz_cid (original terminal entry point)
# ---------------------------------------------------------------------------

def get_xyz_cid_cli():
    """get_xyz_cid <cid_or_name> — print descriptor table in the terminal."""
    import argparse  # noqa: PLC0415

    parser = argparse.ArgumentParser(
        description="Fetch and display molecular descriptors from PubChem"
    )
    parser.add_argument(
        "input", type=str, help="PubChem CID (number) or compound name"
    )
    args = parser.parse_args()
    interact_with_pubchem(args.input)
