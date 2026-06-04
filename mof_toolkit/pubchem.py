"""
pubchem.py
----------
Fetch chemical information and compute molecular descriptors
from PubChem CIDs, compound names, or SMILES strings.

CLI commands
------------
  pubchem_interactive                    — open a local web viewer in the browser
  pubchem_check_prop <cid|name|smiles> <func|key>
                                         — evaluate one RDKit function on a compound
  rdkit_default_props <cid|name|smiles>  — print all default descriptors + SMILES
  pubchem_batch_fetcher <in.csv> <out.csv>
                                         — batch-compute descriptors from a CSV
  fetch_xyz_batch <in.csv> <out_dir>     — batch-download XYZ files from PubChem

Python-only helpers (use in scripts/notebooks, not CLI)
--------------------------------------------------------
  get_xyz_cid(cid_or_name)  — returns property dict for one compound
  resolve_compound_input(q) — resolve CID/name/SMILES → {cid, name, smiles, mol}
"""

import csv
import json
import os
import subprocess
import sys
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


def fetch_smiles_from_cid(cid: int) -> str | None:
    """
    Fetch canonical SMILES for a given PubChem CID.
    Returns None on any failure so batch runs never crash.
    """
    try:
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
            "/property/CanonicalSMILES/TXT"
        )
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            return response.text.strip()
        print(f"  Warning: PubChem returned status {response.status_code} for CID {cid}")
    except requests.exceptions.Timeout:
        print(f"  Warning: request timed out for CID {cid}")
    except requests.exceptions.ConnectionError:
        print("  Warning: no network connection")
    return None


def fetch_pubchem_metadata(cid: int) -> dict:
    """Fetch IUPAC name, formula, MW, and SMILES for a CID from PubChem."""
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
        "/property/IUPACName,MolecularWeight,MolecularFormula,CanonicalSMILES/JSON"
    )
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    p = response.json()["PropertyTable"]["Properties"][0]
    return {
        "CID":                 cid,
        "Name":                p.get("IUPACName", ""),
        # "MolecularFormula":    p.get("MolecularFormula", ""),
        # "MolecularWeight_API": p.get("MolecularWeight", ""),
        # "SMILES":              p.get("CanonicalSMILES", ""),
    }


def _lookup_by_smiles(smiles: str) -> tuple[int | None, str | None]:
    """
    Query PubChem for a CID matching the given SMILES.
    Returns (cid, iupac_name) or (None, None) if not found.
    """
    try:
        r = requests.get(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/TXT",
            params={"smiles": smiles},
            timeout=10,
        )
        if r.status_code != 200:
            return None, None
        cid = int(r.text.strip().splitlines()[0])
        meta = fetch_pubchem_metadata(cid)
        return cid, meta.get("Name", "")
    except Exception:
        return None, None


def smiles_to_mol(smiles: str) -> Chem.Mol | None:
    """Convert SMILES to RDKit Mol. Returns None on failure."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  Warning: RDKit could not parse SMILES: '{smiles}'")
    return mol


def resolve_to_cid(input_strg: str) -> int:
    """Accept a numeric CID string or a compound name and return the CID."""
    if str(input_strg).strip().isdigit():
        return int(str(input_strg).strip())
    return fetch_cid_from_name(str(input_strg))


def resolve_compound_input(query: str) -> dict:
    """
    Resolve a CID, compound name, or SMILES string into a unified dict.

    Returns
    -------
    dict with keys:
      cid    (int | None)      — None when the compound is not in PubChem
      name   (str)             — IUPAC name, or '' if unavailable
      smiles (str | None)      — canonical SMILES
      mol    (Chem.Mol | None) — parsed RDKit molecule

    When the input is a valid SMILES but PubChem has no record, a warning
    is printed and cid/name are returned as None/''.
    """
    q = query.strip()

    # 1. Numeric string → CID
    if q.isdigit():
        cid = int(q)
        smiles = fetch_smiles_from_cid(cid)
        try:
            name = fetch_pubchem_metadata(cid).get("Name", "")
        except Exception:
            name = ""
        mol = smiles_to_mol(smiles) if smiles else None
        return {"cid": cid, "name": name, "smiles": smiles, "mol": mol}

    # 2. Valid SMILES?
    mol_test = Chem.MolFromSmiles(q)
    if mol_test is not None:
        canonical = Chem.MolToSmiles(mol_test)
        cid, name = _lookup_by_smiles(q)
        if cid is None:
            print(
                f"  Warning: SMILES '{q}' is valid but was not found in PubChem "
                "— CID and Name will be empty."
            )
        return {"cid": cid, "name": name or "", "smiles": canonical, "mol": mol_test}

    # 3. Compound name → PubChem lookup
    cid = fetch_cid_from_name(q)   # raises ValueError if not found
    smiles = fetch_smiles_from_cid(cid)
    try:
        name = fetch_pubchem_metadata(cid).get("Name", "")
    except Exception:
        name = q
    mol = smiles_to_mol(smiles) if smiles else None
    return {"cid": cid, "name": name, "smiles": smiles, "mol": mol}


# ---------------------------------------------------------------------------
# RDKit function resolver
# ---------------------------------------------------------------------------

def _resolve_rdkit_func(func_str: str):
    """
    Resolve a dotted RDKit expression string to a callable via eval.
    Supported namespaces: Chem, rdMolDescriptors, Fragments,
    Descriptors, GraphDescriptors.
    Raises ValueError if unresolvable, TypeError if not callable.
    """
    from rdkit.Chem import Descriptors, GraphDescriptors

    eval_ctx = {
        "Chem":             Chem,
        "rdMolDescriptors": rdMolDescriptors,
        "Fragments":        Fragments,
        "Descriptors":      Descriptors,
        "GraphDescriptors": GraphDescriptors,
    }
    try:
        func = eval(func_str, {"__builtins__": {}}, eval_ctx)
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

def compute_properties(cid: int, properties: dict = None) -> dict | None:
    """
    Compute RDKit descriptors for a CID.
    Returns None if SMILES cannot be fetched or parsed.
    """
    if properties is None:
        properties = DEFAULT_PROPERTIES
    smiles = fetch_smiles_from_cid(cid)
    if smiles is None:
        return None
    mol = smiles_to_mol(smiles)
    if mol is None:
        return None
    return {name: func(mol) for name, func in properties.items()}


# ---------------------------------------------------------------------------
# get_xyz_cid — Python helper, not a CLI command
# ---------------------------------------------------------------------------

def get_xyz_cid(cid_or_name, properties: dict = None) -> dict | None:
    """
    Return a property dict for a single compound.

    This is a Python helper intended for use in scripts and notebooks —
    it is NOT registered as a CLI command. Use it in loops to build
    lists of results, then display or save them as needed.

    Parameters
    ----------
    cid_or_name : int or str
        PubChem CID (integer or numeric string) or compound name.
    properties : dict, optional
        Custom {name: callable(mol)} mapping. Defaults to DEFAULT_PROPERTIES.

    Returns
    -------
    dict with 'CID', 'Name' and all computed descriptor keys,
    or None if the compound could not be fetched.

    Examples
    --------
    >>> from mof_toolkit.pubchem import get_xyz_cid
    >>> props = get_xyz_cid(3033)
    >>> props = get_xyz_cid("aspirin")
    >>> props = get_xyz_cid("cannabidiol")
    """
    try:
        cid = resolve_to_cid(str(cid_or_name))
    except Exception as e:
        print(f"  Error resolving '{cid_or_name}': {e}")
        return None

    # Fetch IUPAC name (best-effort, don't crash if it fails)
    try:
        name = fetch_pubchem_metadata(cid).get("Name", "")
    except Exception:
        name = ""

    props = compute_properties(cid, properties)
    if props is None:
        return None

    return {"CID": cid, "Name": name, **props}


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


# ---------------------------------------------------------------------------
# XYZ file writer
# ---------------------------------------------------------------------------

def fetch_and_save_xyz(cid: int, abbreviation: str, output_dir: str):
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
# Interactive viewer — Flask + 3Dmol.js
# ---------------------------------------------------------------------------

_VIEWER_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>PubChem Viewer</title>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
  <style>
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: sans-serif; background: #f5f5f5; padding: 24px; }
    h2 { margin-bottom: 16px; color: #333; }
    .row { display: flex; gap: 12px; margin-bottom: 20px; align-items: center; }
    input { padding: 8px 12px; font-size: 15px; border: 1px solid #ccc;
            border-radius: 6px; width: 320px; }
    button { padding: 8px 20px; font-size: 15px; background: #4a6fa5;
             color: white; border: none; border-radius: 6px; cursor: pointer; }
    button:hover { background: #3a5a8a; }
    #error { color: #c0392b; margin-bottom: 12px; min-height: 20px; }
    #viewer { width: 100%; height: 420px; border-radius: 8px;
              border: 1px solid #ddd; background: white; position: relative; }
    table { border-collapse: collapse; margin-top: 20px; background: white;
            border-radius: 8px; overflow: hidden; box-shadow: 0 1px 4px #0001; }
    th, td { padding: 8px 18px; text-align: left; font-size: 14px; }
    th { background: #4a6fa5; color: white; }
    tr:nth-child(even) { background: #f0f4fa; }
    #spinner { display: none; margin-left: 12px; color: #888; font-size: 14px; }
  </style>
</head>
<body>
  <h2>PubChem Interactive Viewer</h2>
  <div class="row">
    <input id="query" type="text" placeholder="CID (e.g. 3033) or name (e.g. aspirin)"
           onkeydown="if(event.key==='Enter') lookup()">
    <button onclick="lookup()">Look up</button>
    <span id="spinner">Loading...</span>
  </div>
  <div id="error"></div>
  <div id="viewer"></div>
  <div id="table"></div>

  <script>
    let viewer = null;

    function lookup() {
      const q = document.getElementById('query').value.trim();
      if (!q) return;
      document.getElementById('error').textContent = '';
      document.getElementById('spinner').style.display = 'inline';
      document.getElementById('table').innerHTML = '';

      fetch('/lookup?q=' + encodeURIComponent(q))
        .then(r => r.json())
        .then(data => {
          document.getElementById('spinner').style.display = 'none';
          if (data.error) {
            document.getElementById('error').textContent = 'Error: ' + data.error;
            return;
          }
          render3D(data.sdf);
          renderTable(data.props);
        })
        .catch(e => {
          document.getElementById('spinner').style.display = 'none';
          document.getElementById('error').textContent = 'Request failed: ' + e;
        });
    }

    function render3D(sdf) {
      const el = document.getElementById('viewer');
      el.innerHTML = '';
      viewer = $3Dmol.createViewer(el, { backgroundColor: 'white' });
      viewer.addModel(sdf, 'sdf');
      viewer.setStyle({}, { stick: {}, sphere: { scale: 0.25 } });
      viewer.zoomTo();
      viewer.render();
    }

    function renderTable(props) {
      let html = '<table><tr><th>Property</th><th>Value</th></tr>';
      for (const [k, v] of Object.entries(props)) {
        const val = typeof v === 'number' && !Number.isInteger(v)
                    ? v.toFixed(4) : v;
        html += `<tr><td>${k}</td><td>${val}</td></tr>`;
      }
      html += '</table>';
      document.getElementById('table').innerHTML = html;
    }
  </script>
</body>
</html>
"""


def pubchem_interactive_cli():
    """
    pubchem_interactive — start a local web viewer.

    Launches a small Flask server at http://localhost:5050
    Open that URL in your browser (Windows browser for WSL users).
    Type a CID or compound name, press Enter or click Look up.
    Press Ctrl+C in the terminal to stop the server.
    """
    try:
        from flask import Flask, request, jsonify
    except ImportError:
        print("Error: Flask is not installed.")
        print("Run: pip install flask")
        sys.exit(1)

    app = Flask(__name__)
    app.logger.disabled = True
    import logging
    log = logging.getLogger("werkzeug")
    log.setLevel(logging.ERROR)   # suppress per-request logs

    @app.route("/")
    def index():
        from flask import Response
        return Response(_VIEWER_HTML, mimetype="text/html")

    @app.route("/lookup")
    def lookup():
        q = request.args.get("q", "").strip()
        if not q:
            return jsonify({"error": "empty query"})

        # Resolve CID
        try:
            cid = resolve_to_cid(q)
        except Exception as e:
            return jsonify({"error": str(e)})

        # Fetch 3D SDF
        url_3d = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
            "/SDF?record_type=3d"
        )
        resp = requests.get(url_3d, timeout=10)
        if resp.status_code != 200:
            return jsonify({"error": f"No 3D structure for CID {cid} "
                                     f"(HTTP {resp.status_code})"})

        # Compute properties
        props_raw = get_xyz_cid(cid)
        if props_raw is None:
            return jsonify({"error": f"Could not compute properties for CID {cid}"})

        # Round floats for JSON
        props_clean = {}
        for k, v in props_raw.items():
            props_clean[k] = round(v, 4) if isinstance(v, float) else v

        return jsonify({"sdf": resp.text, "props": props_clean})

    port = 5050
    print(f"\n  PubChem viewer running at: http://localhost:{port}")
    print("  WSL users: open that URL in your Windows browser.")
    print("  Press Ctrl+C to stop.\n")
    app.run(host="0.0.0.0", port=port, debug=False)


# ---------------------------------------------------------------------------
# CLI: pubchem_check_prop
# ---------------------------------------------------------------------------

def pubchem_check_prop_cli():
    """pubchem_check_prop <cid|name|smiles> <func|key>

    Evaluate a single RDKit property on a compound and print a labelled result.
    The second argument can be either:
      • A key from DEFAULT_PROPERTIES  (e.g. NumRings, TPSA, HBA)
      • Any dotted RDKit callable      (e.g. rdMolDescriptors.CalcTPSA)

    The compound can be given as a CID, an IUPAC/common name, or a SMILES string.

    Examples
    --------
      pubchem_check_prop 3033 TPSA
      pubchem_check_prop aspirin NumRings
      pubchem_check_prop "CC(=O)O" rdMolDescriptors.CalcTPSA
      pubchem_check_prop 2244 Fragments.fr_COO
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Evaluate one RDKit property on a compound",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "DEFAULT_PROPERTIES keys: " + ", ".join(DEFAULT_PROPERTIES) + "\n\n"
            "Supported RDKit namespaces: Chem, rdMolDescriptors, Fragments,\n"
            "Descriptors, GraphDescriptors\n\n"
            "Examples:\n"
            "  pubchem_check_prop 3033 TPSA\n"
            "  pubchem_check_prop aspirin rdMolDescriptors.CalcNumAromaticRings\n"
            '  pubchem_check_prop "CC(=O)O" HBA'
        ),
    )
    parser.add_argument(
        "compound", help="PubChem CID, compound name, or SMILES string"
    )
    parser.add_argument(
        "function",
        help="DEFAULT_PROPERTIES key (e.g. TPSA) or RDKit callable string",
    )
    args = parser.parse_args()

    # Resolve compound
    try:
        resolved = resolve_compound_input(args.compound)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    mol = resolved["mol"]
    if mol is None:
        print("Error: could not obtain a valid molecule.")
        sys.exit(1)

    # Resolve function — accept DEFAULT_PROPERTIES key or dotted RDKit path
    func_label = args.function
    if args.function in DEFAULT_PROPERTIES:
        func = DEFAULT_PROPERTIES[args.function]
    else:
        try:
            func = _resolve_rdkit_func(args.function)
        except (ValueError, TypeError) as e:
            print(f"Error: {e}")
            print(
                "Hint: use a DEFAULT_PROPERTIES key "
                f"({', '.join(DEFAULT_PROPERTIES)}) "
                "or a dotted RDKit path like rdMolDescriptors.CalcTPSA."
            )
            sys.exit(1)

    try:
        result = func(mol)
    except TypeError as e:
        print(f"Error calling {func_label}(mol): {e}")
        print(
            "Hint: this function may require extra arguments beyond mol. "
            "Check the RDKit docs for its full signature."
        )
        sys.exit(1)

    cid_display = str(resolved["cid"]) if resolved["cid"] else "N/A (not in PubChem)"
    val_display = f"{result:.6f}" if isinstance(result, float) else str(result)

    print(f"\nName     : {resolved['name'] or 'N/A'}")
    print(f"CID      : {cid_display}")
    print(f"SMILES   : {resolved['smiles'] or 'N/A'}")
    print(f"Function : {func_label}")
    print(f"Value    : {val_display}\n")


# ---------------------------------------------------------------------------
# CLI: rdkit_default_props
# ---------------------------------------------------------------------------

def rdkit_default_props_cli():
    """rdkit_default_props <cid|name|smiles>

    Print CID, Name, SMILES, and all DEFAULT_PROPERTIES for a compound.
    Accepts a PubChem CID, common/IUPAC name, or a SMILES string.
    If the SMILES is valid but not in PubChem, CID and Name are shown as N/A.

    Examples
    --------
      rdkit_default_props 2244
      rdkit_default_props aspirin
      rdkit_default_props "CC(=O)Oc1ccccc1C(=O)O"
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Print default RDKit descriptors for a compound",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  rdkit_default_props 2244\n"
            "  rdkit_default_props aspirin\n"
            '  rdkit_default_props "CC(=O)Oc1ccccc1C(=O)O"'
        ),
    )
    parser.add_argument(
        "compound", help="PubChem CID, compound name, or SMILES string"
    )
    args = parser.parse_args()

    try:
        resolved = resolve_compound_input(args.compound)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    mol = resolved["mol"]
    if mol is None:
        print("Error: could not obtain a valid molecule.")
        sys.exit(1)

    cid_display = str(resolved["cid"]) if resolved["cid"] else "N/A (not in PubChem)"
    row: dict = {
        "CID":    cid_display,
        "Name":   resolved["name"] or "N/A",
        "SMILES": resolved["smiles"] or "N/A",
    }
    for prop_name, func in DEFAULT_PROPERTIES.items():
        try:
            val = func(mol)
            row[prop_name] = f"{val:.4f}" if isinstance(val, float) else val
        except Exception as e:
            row[prop_name] = f"ERROR: {e}"

    display_table(row)


# ---------------------------------------------------------------------------
# CLI: pubchem_batch_fetcher
# ---------------------------------------------------------------------------

def pubchem_batch_fetcher_cli():
    """pubchem_batch_fetcher <input.csv> <output.csv>

    Input CSV format
    ----------------
    Required column  : CID
    Optional metadata: Name, Guest_Type, or any other column — passed through.
                       If Name is blank, IUPAC name is fetched from PubChem.
    Extra properties : column header format  PropName::rdkit.func.Path
                       e.g.  Chi0v::Chem.rdMolDescriptors.CalcChi0v

    Output CSV
    ----------
    All input metadata columns + DEFAULT_PROPERTIES + any extra computed columns.
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Batch-compute molecular descriptors from a CSV of PubChem CIDs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Extra property columns use '::' in the header:\n"
            "  Chi0v::Chem.rdMolDescriptors.CalcChi0v\n"
            "  NumHeteroatoms::rdMolDescriptors.CalcNumHeteroatoms"
        ),
    )
    parser.add_argument("input_csv",  help="Input CSV file (must have a CID column)")
    parser.add_argument("output_csv", help="Output CSV file path")
    args = parser.parse_args()

    if not os.path.isfile(args.input_csv):
        print(f"Error: input file not found: {args.input_csv}")
        sys.exit(1)

    with open(args.input_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        headers = reader.fieldnames or []
        rows    = list(reader)

    has_cid    = "CID"    in headers
    has_smiles = "SMILES" in headers

    if not has_cid and not has_smiles:
        print("Error: input CSV must have a 'CID' column (or a 'SMILES' column).")
        sys.exit(1)

    # Split headers: system columns, extra computed cols (PropName::func), passthrough
    _system = {"CID", "SMILES"}
    extra_props:      dict[str, str]  = {}
    passthrough_cols: list[str]       = []
    for h in headers:
        if "::" in h:
            prop_name, func_str = h.split("::", 1)
            extra_props[prop_name.strip()] = func_str.strip()
        elif h not in _system:
            passthrough_cols.append(h)

    resolved_extra: dict[str, object] = {}
    for prop_name, func_str in extra_props.items():
        try:
            resolved_extra[prop_name] = _resolve_rdkit_func(func_str)
        except (ValueError, TypeError) as e:
            print(f"Warning: skipping extra property '{prop_name}': {e}")

    output_rows: list[dict] = []
    total = len(rows)

    for i, row in enumerate(rows, 1):
        cid_raw    = (row.get("CID",    "") or "").strip()
        smiles_raw = (row.get("SMILES", "") or "").strip()

        if not cid_raw and not smiles_raw:
            print(f"  [{i}/{total}] Skipping row — no CID or SMILES")
            continue

        print(f"  [{i}/{total}] {cid_raw or smiles_raw[:30]} ...", end=" ", flush=True)

        cid   = None
        name  = ""
        mol   = None
        smiles = None

        if cid_raw:
            # Primary: resolve from CID
            try:
                cid = resolve_to_cid(cid_raw)
            except Exception as e:
                print(f"ERROR resolving CID ({e})")
                continue
            smiles = smiles_raw or fetch_smiles_from_cid(cid)
            if smiles is None:
                print("ERROR (could not fetch SMILES)")
                continue
            mol = smiles_to_mol(smiles)
            if mol is None:
                print("ERROR (could not parse SMILES)")
                continue
        else:
            # Fallback: resolve from SMILES
            mol = Chem.MolFromSmiles(smiles_raw)
            if mol is None:
                print("ERROR (invalid SMILES)")
                continue
            smiles = Chem.MolToSmiles(mol)
            cid, name = _lookup_by_smiles(smiles_raw)
            if cid is None:
                print("(not in PubChem — CID/Name will be empty)", end=" ")

        # Auto-fetch name from PubChem when we have a CID and no name yet
        if cid and not name:
            try:
                name = fetch_pubchem_metadata(cid).get("Name", "")
            except Exception:
                name = ""

        out_row: dict = {
            "CID":    cid if cid else "",
            "SMILES": smiles,
        }

        # Passthrough metadata columns
        for col in passthrough_cols:
            val = row.get(col, "")
            if col.strip().lower() == "name" and not val.strip():
                val = name
            out_row[col] = val

        # Default descriptors
        for prop_name, func in DEFAULT_PROPERTIES.items():
            try:
                out_row[prop_name] = func(mol)
            except Exception as e:
                out_row[prop_name] = f"ERROR: {e}"

        # Extra user-defined descriptors
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

    print(f"\nSaved {len(output_rows)}/{total} rows → {args.output_csv}")


# ---------------------------------------------------------------------------
# CLI: fetch_xyz_batch
# ---------------------------------------------------------------------------

def fetch_xyz_batch_cli():
    """fetch_xyz_batch <input.csv> <output_dir>

    Batch-download PubChem 3D conformers as .xyz files.
    Input CSV must have CID and Abbreviation columns.
    Each compound is saved as <Abbreviation>.xyz in output_dir.

    Example
    -------
      fetch_xyz_batch tests/data/molecules_example.csv ./xyz_files/
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Batch-download PubChem 3D conformers as XYZ files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n  fetch_xyz_batch molecules.csv ./xyz_files/",
    )
    parser.add_argument("input_csv",  help="CSV with CID and Abbreviation columns")
    parser.add_argument("output_dir", help="Directory to write .xyz files into")
    args = parser.parse_args()

    if not os.path.isfile(args.input_csv):
        print(f"Error: input file not found: {args.input_csv}")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.input_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        headers = reader.fieldnames or []
        rows = list(reader)

    if "CID" not in headers:
        if "SMILES" in headers:
            print(
                "Error: this command requires a CID column to download PubChem 3D conformers.\n"
                "For SMILES-only 3D structure generation, use:\n"
                "  smiles_to_3d <SMILES> [--format xyz sdf pdb mol]"
            )
        else:
            print("Error: input CSV must have a 'CID' column.")
        sys.exit(1)

    total = len(rows)
    for i, row in enumerate(rows, 1):
        cid_raw = row.get("CID", "").strip()
        abbr    = row.get("Abbreviation", cid_raw).strip()
        if not cid_raw:
            print(f"  [{i}/{total}] Skipping row with empty CID")
            continue
        print(f"  [{i}/{total}] {abbr} (CID {cid_raw})")
        try:
            fetch_and_save_xyz(int(cid_raw), abbr, args.output_dir)
        except Exception as e:
            print(f"    ERROR: {e}")
