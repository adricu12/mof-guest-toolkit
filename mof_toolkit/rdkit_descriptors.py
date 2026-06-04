"""
rdkit_descriptors.py
-------------------
Fetch chemical information and compute molecular descriptors
from PubChem CIDs, compound names, or SMILES strings.

CLI commands
------------
  pubchem_interactive                      — open a local web viewer in the browser
  rdkit_check_descript <cid|name|smiles> <func|key>
                                           — evaluate one RDKit function on a compound
  rdkit_default_descrpts <cid|name|smiles>    — print all default descriptors + identifiers
  rdkit_batch_fetcher <in.csv> <out.csv>   — batch-compute descriptors from a CSV
  fetch_xyz_batch <in.csv> <out_dir>       — batch-generate 3D structure files

Python-only helpers (use in scripts/notebooks, not CLI)
--------------------------------------------------------
  get_rdkit_dict(query, properties=None)
      Returns a descriptor dict for one compound (CID, name, or SMILES).
  resolve_compound_input(query)
      Resolves CID / name / SMILES → {cid, iupac_name, common_name, smiles, mol}.
"""

import csv
import os
import sys

import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import Fragments, rdMolDescriptors


# ---------------------------------------------------------------------------
# Descriptor registry
# ---------------------------------------------------------------------------

DEFAULT_DESCRIPTORS = {
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
    """
    Fetch IUPAC name, common name (first synonym), canonical SMILES,
    molecular formula, and molecular weight for a CID from PubChem.

    Returns a dict with keys:
      CID, IUPAC_Name, Common_Name, SMILES, MolecularFormula, MolecularWeight_PubChem
    """
    # Fetch properties
    prop_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
        "/property/IUPACName,MolecularWeight,MolecularFormula,CanonicalSMILES/JSON"
    )
    response = requests.get(prop_url, timeout=10)
    response.raise_for_status()
    p = response.json()["PropertyTable"]["Properties"][0]

    # Fetch first synonym as common name
    common_name = ""
    try:
        syn_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
            "/synonyms/JSON"
        )
        syn_resp = requests.get(syn_url, timeout=10)
        if syn_resp.status_code == 200:
            synonyms = syn_resp.json()["InformationList"]["Information"][0].get(
                "Synonym", []
            )
            # First synonym that is not just digits or a registry number
            for s in synonyms:
                if not s.isdigit() and "-" not in s[:5]:
                    common_name = s
                    break
            if not common_name and synonyms:
                common_name = synonyms[0]
    except Exception:
        pass

    return {
        "CID":                    cid,
        "IUPAC_Name":             p.get("IUPACName", ""),
        "Common_Name":            common_name,
        "SMILES":                 p.get("CanonicalSMILES", ""),
        "MolecularFormula":       p.get("MolecularFormula", ""),
        "MolecularWeight_PubChem": p.get("MolecularWeight", ""),
    }


def _lookup_by_smiles(smiles: str) -> tuple[int | None, str | None, str | None]:
    """
    Query PubChem for a CID matching the given SMILES.
    Returns (cid, iupac_name, common_name) or (None, None, None) if not found.
    """
    try:
        r = requests.get(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/TXT",
            params={"smiles": smiles},
            timeout=10,
        )
        if r.status_code != 200:
            return None, None, None
        cid = int(r.text.strip().splitlines()[0])
        meta = fetch_pubchem_metadata(cid)
        return cid, meta.get("IUPAC_Name", ""), meta.get("Common_Name", "")
    except Exception:
        return None, None, None


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
      cid         (int | None)      — None when the compound is not in PubChem
      iupac_name  (str)             — IUPAC name from PubChem, or ''
      common_name (str)             — first synonym from PubChem, or ''
      smiles      (str | None)      — canonical SMILES
      mol         (Chem.Mol | None) — parsed RDKit molecule

    When the input is a valid SMILES but PubChem has no record, a warning
    is printed and cid/iupac_name/common_name are returned as None/''/'' — but
    the mol and smiles are still populated so RDKit descriptors can be computed.
    """
    q = query.strip()

    # 1. Numeric string → CID
    if q.isdigit():
        cid = int(q)
        smiles = fetch_smiles_from_cid(cid)
        try:
            meta = fetch_pubchem_metadata(cid)
            iupac_name  = meta.get("IUPAC_Name", "")
            common_name = meta.get("Common_Name", "")
        except Exception:
            iupac_name = common_name = ""
        mol = smiles_to_mol(smiles) if smiles else None
        return {
            "cid": cid, "iupac_name": iupac_name, "common_name": common_name,
            "smiles": smiles, "mol": mol,
        }

    # 2. Valid SMILES?
    mol_test = Chem.MolFromSmiles(q)
    if mol_test is not None:
        canonical = Chem.MolToSmiles(mol_test)
        cid, iupac_name, common_name = _lookup_by_smiles(q)
        if cid is None:
            print(
                f"  Warning: SMILES '{q}' is valid but was not found in PubChem "
                "— CID, IUPAC_Name, and Common_Name will be empty."
            )
        return {
            "cid": cid, "iupac_name": iupac_name or "",
            "common_name": common_name or "", "smiles": canonical, "mol": mol_test,
        }

    # 3. Compound name → PubChem lookup
    cid = fetch_cid_from_name(q)   # raises ValueError if not found
    smiles = fetch_smiles_from_cid(cid)
    try:
        meta = fetch_pubchem_metadata(cid)
        iupac_name  = meta.get("IUPAC_Name", "")
        common_name = meta.get("Common_Name", "")
    except Exception:
        iupac_name = ""
        common_name = q
    mol = smiles_to_mol(smiles) if smiles else None
    return {
        "cid": cid, "iupac_name": iupac_name, "common_name": common_name,
        "smiles": smiles, "mol": mol,
    }


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

def _compute_descriptors(mol: Chem.Mol, descriptors: dict) -> dict:
    """Compute all descriptors in the given dict for the given mol."""
    result = {}
    for desc_name, func in descriptors.items():
        try:
            result[desc_name] = func(mol)
        except Exception as e:
            result[desc_name] = f"ERROR: {e}"
    return result


# ---------------------------------------------------------------------------
# get_rdkit_dict — Python helper, not a CLI command
# ---------------------------------------------------------------------------

def get_rdkit_dict(query, descriptors: dict = None) -> dict | None:
    """
    Return a descriptor dict for a single compound.

    Accepts a PubChem CID (int or numeric string), a compound name, or a
    SMILES string. When a valid SMILES is given that is not in PubChem,
    CID / IUPAC_Name / Common_Name are left empty but all RDKit descriptors
    are still computed from the SMILES.

    Parameters
    ----------
    query : int or str
        PubChem CID, compound name, or SMILES string.
    descriptors : dict, optional
        Custom {name: callable(mol)} mapping. Defaults to DEFAULT_DESCRIPTORS.

    Returns
    -------
    dict with keys CID, IUPAC_Name, Common_Name, SMILES, and all descriptor
    keys; or None if no valid molecule could be obtained at all.

    Examples
    --------
    >>> from mof_toolkit.rdkit_descriptors import get_rdkit_dict, display_table
    >>> display_table(get_rdkit_dict(3033))
    >>> display_table(get_rdkit_dict("aspirin"))
    >>> display_table(get_rdkit_dict("cannabidiol"))
    >>> display_table(get_rdkit_dict("CC(=O)Oc1ccccc1C(=O)O"))   # SMILES input
    >>> # SMILES not in PubChem — CID/names empty, descriptors still computed:
    >>> display_table(get_rdkit_dict("C1CC1"))
    """
    if descriptors is None:
        descriptors = DEFAULT_DESCRIPTORS

    try:
        resolved = resolve_compound_input(str(query))
    except Exception as e:
        print(f"  Error resolving '{query}': {e}")
        return None

    mol = resolved["mol"]
    if mol is None:
        print(f"  Error: could not obtain a valid molecule for '{query}'.")
        return None

    descriptor_vals = _compute_descriptors(mol, descriptors)

    return {
        "CID":         resolved["cid"] if resolved["cid"] else "",
        "IUPAC_Name":  resolved["iupac_name"],
        "Common_Name": resolved["common_name"],
        "SMILES":      resolved["smiles"] or "",
        **descriptor_vals,
    }


# ---------------------------------------------------------------------------
# Terminal display
# ---------------------------------------------------------------------------

def display_table(descriptor_dict: dict):
    """Print a descriptor dict as an aligned table."""
    if descriptor_dict is None:
        print("  (no data)")
        return
    print(f"\n{'Descriptor':<24} {'Value'}")
    print("-" * 56)
    for prop, value in descriptor_dict.items():
        if isinstance(value, float):
            value = f"{value:.4f}"
        print(f"{prop:<24} {value}")


# ---------------------------------------------------------------------------
# XYZ file writer (single compound, kept for fetch_and_save_xyz)
# ---------------------------------------------------------------------------

def fetch_and_save_xyz(cid: int, file_stem: str, output_dir: str) -> bool:
    """
    Download the PubChem 3D conformer for a CID and write it as
    <file_stem>.xyz in output_dir.

    Returns True on success, False if no 3D conformer is available.
    """
    compounds = pcp.get_compounds(cid, "cid", record_type="3d")
    if not compounds:
        print(f"  No 3D conformer available for CID={cid} — skipping")
        return False
    atoms = compounds[0].atoms
    lines = [
        str(len(atoms)),
        f"{file_stem}  CID={cid}  source: PubChem 3D conformer",
    ]
    for atom in atoms:
        lines.append(
            f"{atom.element:<3}  {atom.x:12.6f}  {atom.y:12.6f}  {atom.z:12.6f}"
        )
    filepath = os.path.join(output_dir, f"{file_stem}.xyz")
    with open(filepath, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"    Saved: {os.path.basename(filepath)}")
    return True


# ---------------------------------------------------------------------------
# Batch row resolver — shared by rdkit_batch_fetcher and fetch_xyz_batch
# ---------------------------------------------------------------------------

def _resolve_batch_row(
    row: dict,
    row_index: int,
    total: int,
) -> dict | None:
    """
    Resolve one CSV row that may contain CID, Name, and/or SMILES columns.

    Resolution priority (when multiple are present and valid):
      SMILES > CID > Name

    Validation and mismatch rules
    ------------------------------
    - Any present value is validated individually first.
    - Invalid values are flagged and ignored.
    - When two or more valid inputs yield different PubChem CIDs, a mismatch
      flag is printed and the higher-priority input is used.
    - At least one valid input must remain after filtering; otherwise the row
      is skipped.

    Returns
    -------
    dict with keys:
      mol, smiles, cid, iupac_name, common_name,
      input_cid, input_name, input_smiles   (original values from row)
    or None if the row cannot be resolved.
    """
    label = f"[{row_index}/{total}]"

    raw_cid    = (row.get("CID",    "") or "").strip()
    raw_name   = (row.get("Name",   "") or "").strip()
    raw_smiles = (row.get("SMILES", "") or "").strip()

    if not raw_cid and not raw_name and not raw_smiles:
        print(f"  {label} SKIP — no CID, Name, or SMILES found in row")
        return None

    # ------------------------------------------------------------------
    # Step 1: validate each field individually
    # ------------------------------------------------------------------
    valid: dict[str, dict] = {}   # key → resolved dict from resolve_compound_input

    if raw_cid:
        try:
            resolved_cid = resolve_compound_input(raw_cid)
            if resolved_cid["mol"] is not None:
                valid["CID"] = resolved_cid
            else:
                print(f"  {label} FLAG — CID '{raw_cid}' resolved but gave no valid molecule; ignoring")
        except Exception as e:
            print(f"  {label} FLAG — CID '{raw_cid}' is invalid ({e}); ignoring")

    if raw_name:
        try:
            resolved_name = resolve_compound_input(raw_name)
            if resolved_name["mol"] is not None:
                valid["Name"] = resolved_name
            else:
                print(f"  {label} FLAG — Name '{raw_name}' resolved but gave no valid molecule; ignoring")
        except Exception as e:
            print(f"  {label} FLAG — Name '{raw_name}' is invalid ({e}); ignoring")

    if raw_smiles:
        mol_test = Chem.MolFromSmiles(raw_smiles)
        if mol_test is None:
            print(f"  {label} FLAG — SMILES '{raw_smiles}' is not a valid SMILES; ignoring")
        else:
            try:
                resolved_smiles = resolve_compound_input(raw_smiles)
                valid["SMILES"] = resolved_smiles
            except Exception as e:
                print(f"  {label} FLAG — SMILES '{raw_smiles}' could not be resolved ({e}); ignoring")

    if not valid:
        print(f"  {label} SKIP — no valid input remained after validation")
        return None

    # ------------------------------------------------------------------
    # Step 2: mismatch check — compare PubChem CIDs of valid inputs
    # ------------------------------------------------------------------
    cids_found = {
        k: v["cid"] for k, v in valid.items() if v["cid"] is not None
    }
    unique_cids = set(cids_found.values())
    if len(unique_cids) > 1:
        details = ", ".join(f"{k}→CID {c}" for k, c in cids_found.items())
        print(f"  {label} MISMATCH — inputs point to different compounds: {details}")

    # ------------------------------------------------------------------
    # Step 3: pick the authoritative input by priority: SMILES > CID > Name
    # ------------------------------------------------------------------
    chosen_key = None
    for preferred in ("SMILES", "CID", "Name"):
        if preferred in valid:
            chosen_key = preferred
            break

    chosen = valid[chosen_key]
    print(
        f"  {label} Using {chosen_key} as primary input "
        f"(CID={chosen['cid'] or 'N/A'})"
    )

    return {
        "mol":         chosen["mol"],
        "smiles":      chosen["smiles"] or "",
        "cid":         chosen["cid"],
        "iupac_name":  chosen["iupac_name"],
        "common_name": chosen["common_name"],
        "input_cid":   raw_cid,
        "input_name":  raw_name,
        "input_smiles": raw_smiles,
        "chosen_key":  chosen_key,
    }


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
    body { font-family: sans-serif; background: #f0f2f5; }

    header {
      background: #1e2a3a;
      padding: 14px 24px;
      display: flex;
      align-items: center;
      gap: 14px;
    }
    .header-mascot {
      height: 52px;
      width: auto;
      border-radius: 6px;
      opacity: 0.92;
    }
    .header-text h2 {
      font-size: 16px; font-weight: 600;
      color: #e8edf2; letter-spacing: 0.3px;
    }
    .header-text .subtitle {
      font-size: 11px; color: #7a9ab5; margin-top: 2px; letter-spacing: 0.3px;
    }

    .main {
      padding: 20px 24px;
      display: flex; gap: 24px;
      align-items: flex-start; flex-wrap: wrap;
    }
    .left-panel { display: flex; flex-direction: column; gap: 12px; flex-shrink: 0; }

    .search-row { display: flex; gap: 0; align-items: stretch; }
    .search-icon-wrap {
      background: #fff;
      border: 1px solid #ced4da; border-right: none;
      border-radius: 6px 0 0 6px;
      padding: 4px 8px;
      display: flex; align-items: center;
    }
    .search-icon-wrap img { height: 26px; width: auto; }
    input {
      padding: 8px 12px; font-size: 14px;
      border: 1px solid #ced4da;
      border-left: none; border-right: none;
      width: 230px; background: #fff; color: #1e2a3a;
      outline: none;
    }
    input:focus { border-color: #4a6fa5; }
    button {
      padding: 8px 16px; font-size: 14px;
      background: #2c4a6e; color: #fff;
      border: none; border-radius: 0 6px 6px 0;
      cursor: pointer; transition: background 0.15s;
      white-space: nowrap;
    }
    button:hover { background: #1e3550; }

    #status-area {
      min-height: 90px;
      display: flex; align-items: center; gap: 12px;
    }
    #status-img { height: 80px; width: auto; display: none; }
    #status-text { font-size: 13px; color: #4a6080; font-style: italic; }
    #error-text { font-size: 13px; color: #c0392b; display: none; }

    #viewer {
      width: 350px; height: 350px;
      border-radius: 8px; border: 1px solid #d0d7e0;
      background: #ffffff; position: relative;
      box-shadow: 0 1px 4px rgba(0,0,0,0.08);
      display: none;
    }

    .right-panel { flex: 1; min-width: 260px; }

    #compound-header {
      display: none;
      background: #1e2a3a; color: #e8edf2;
      border-radius: 8px 8px 0 0;
      padding: 12px 16px;
    }
    #compound-header .cid {
      font-size: 11px; color: #7a9ab5;
      text-transform: uppercase; letter-spacing: 0.5px;
    }
    #compound-header .cname { font-size: 15px; font-weight: 600; margin: 2px 0; }
    #compound-header .iupac { font-size: 11px; color: #9ab4cc; }

    table {
      width: 100%; border-collapse: collapse;
      background: #fff; border-radius: 0 0 8px 8px;
      overflow: hidden; border: 1px solid #d0d7e0;
      border-top: none; font-size: 13px;
    }
    th {
      background: #2c4a6e; color: #c8d8ea;
      padding: 8px 14px; text-align: left;
      font-weight: 500; font-size: 11px;
      text-transform: uppercase; letter-spacing: 0.4px;
    }
    td { padding: 7px 14px; color: #1e2a3a; border-bottom: 1px solid #eef0f3; }
    td:first-child { color: #4a6080; font-weight: 500; width: 46%; }
    tr:last-child td { border-bottom: none; }
    tr:nth-child(even) td { background: #f7f9fc; }

    @keyframes spin { from { transform: rotate(0deg); } to { transform: rotate(360deg); } }
    .spinning { animation: spin 1.2s linear infinite; transform-origin: center; }
  </style>
</head>
<body>
  <header>
    <img class="header-mascot" src="data:image/jpeg;base64,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" alt="mascot">
    <div class="header-text">
      <div class="subtitle">mof-guest-toolkit · PubChem Interactive Viewer</div>
      <h2>Molecular Descriptor Explorer</h2>
    </div>
  </header>

  <div class="main">
    <div class="left-panel">

      <div class="search-row">
        <div class="search-icon-wrap">
          <img src="data:image/jpeg;base64,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" alt="search">
        </div>
        <input id="query" type="text"
               placeholder="CID, name, or SMILES"
               onkeydown="if(event.key===\'Enter\') lookup()">
        <button onclick="lookup()">Look up</button>
      </div>

      <div id="status-area">
        <img id="status-img" src="" alt="status">
        <div>
          <div id="status-text">Enter a CID, compound name, or SMILES above.</div>
          <div id="error-text"></div>
        </div>
      </div>

      <div id="viewer"></div>
    </div>

    <div class="right-panel">
      <div id="compound-header"></div>
      <div id="table"></div>
    </div>
  </div>

  <script>
    const IMG_LOADING = "data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEBkAGQAAD/2wBDAAMCAgMCAgMDAwMEAwMEBQgFBQQEBQoHBwYIDAoMDAsKCwsNDhIQDQ4RDgsLEBYQERMUFRUVDA8XGBYUGBIUFRT/2wBDAQMEBAUEBQkFBQkUDQsNFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBT/wgARCAE8AagDASIAAhEBAxEB/8QAHQABAAMBAQEBAQEAAAAAAAAAAAYHCAUEAwkCAf/EABsBAQACAwEBAAAAAAAAAAAAAAAEBQIDBgEH/9oADAMBAAIQAxAAAAHVIAAAAAAAAAAAAAAAAAAAABF8dsoGWqN1Hz+5V91egtOFAAAAAAAAAAAAAAAAAzfWfYuCh+tdS0Mb6Un8lNBP5IePz32D3wD55z0g0WmH55cGQ6X6b9NMZj+mm034zfpDoPjob6wAAAAAAAed56A9AAAA4ee9AYnqfoG+Kbp95nbF2ef0WPGPP6G2BB6fvTP9b2+pBZ8OyXrTNcDrdGMd2R5noQWHHgZT8+oOpXdnieaaIynE6DZAvPlgAAADzfXMcuN6IVa172UDIWqKl+eOV5+PzUbXzYdzrutuygYq2jlvQ+vOQkCqrKpNMZf1BPhBWWEPpfQkDg9TnON2RW9H9U6m1ML6sncpagvPloh+O6r+pUexK3tvp4/Zz7Lh6LovaHHq++yX3LsuzXLicsLj5wGWkB4/Y89D3w4/Y8zD3AADNkB4tx9bzeiOd0cjc9dTr+86ayuauLfKudext3spyXYtwzmWuMRb6262a9KYl0b9KWN8/pWTg17AMv3RMPRFvqnz/tjz65kXmHD7Eqg5+P8AuXhWdvJJIWvBBlqef0PPfn9D3wAAADn/AD4eV4PU9zXn5/2ZA6zWiNrv5fJBlqcbsx7PDNetct6ksoL8/P0DwLKj/PeWd9KadmMtm462L57UOVt54+mRepsGnLjrZ9c1TF9Kyo0vFHcAAAAKLvTP8TofnoTz+jZBDfWAAET5eqfYDn9DZCD3EAD55D1Zz4fRU/0L0+Z9BM50BHpD88sct6oxPrO3q+zjKU2ltwsXs8/oUlvi7aOfKtvqXWGTvPqRlM4JYWPq6Z7tZwyZ+ZBAmAAAAAAAAcvK/QsSp+hV/wCzxxeLezjTEDnlt89CVRAEboeNdaIq+l7YidFD+HalZxbyYKLGyvH/AH5J/Jf1LuDJ91fE6A1j97jl6Ivcj78x3XVMts6+2hS2wFe07HtKXtPMRRXAAEakub4TcVmxlM92FKsl4/ZEkA9AAc/oVvrmUHsTO+iInQx+QE3lgy1io9cyxM7w/QFb2lJ3xYiVRK/sDzyKbA8k15IKzu6jW4s+GDZDUjd1ay4vHuLMek9uHpZ10VH3Z8ubuvfAjSEFnWObGF0daw6Y4ehBmPh94PhJrO9Mvx+p+gbQgNQ3JZ8hQ9WfoG6Dj8RXFcda5Y2d1MZbFgzPQIMwBT9wU/GuvnclP3A9CTSDO+my+fD7GkIHV+P2FpwYe+AAAOf0OHlhzpbi7Ys6JnSvdst2qlbqK+aGncItljDq5g2yLuo6IorkeP3z2ZvvzLcDrPHMJhy4fRTSm41pXOLm+6LIzTvh+XXGBbm7P5X/ALYEUl8ffZ4pLYBT9wU/GulwU/cAEmlr+h+pelV38wFrwBHPXhI+lZ/Sp4HW60FlxQAAFU0hsWqbassr2Yx0zq2zVHvHClS1UNMTol15/nWkJkbjSko7YMcmeND4itq22pXXdn65lKasxx8aD6BOeH49OQbeh7k8+f8AfWWLnu2+j1fBeLheqLzavcaDzjmLoM9Cn7gjemyr+5IfMPMg31eT9IezqRL8JdBW+X90V3XdnkP0ef6UX1eYWB7NAXHzULb5+AABBad04mxMttSN+miLWkKJJCNIETywljOddWcHWmRuzpCRoxk0PVtpXR7Rma77hzqummjKr4b6tYHjyn2PM7kpe4Lc2Q8oWXLKH0Wf05f+3fqsvnaGd6/k0ezFf2BY8YGccAAD54/1xG4fRweyJA3V1P3Ae6g2wgAADI3Lt6vZrFXoy82axd9DTlaRewsvKY6moephlVtlehWTw1bcbwDUGcOx5aSaAyVsGLKsbn9CB8p1NF6wzfpCJ0Tx/TP8iovBj/ZmmxpOi92RPTYZ0vPx5702vYlnc+nmV0ejGd2TOcuB8/pO5UAAAAAAAAAAeP3z2KEr2yga94+Qe7v1aE8VU9XH2Q8Lo9vLypudpDl5YV7Tdjp0XRHY4/Y5i/V/YFfxrSv70zfF67stiVn0J5Ycdn+3ZG82kXlG6tyvPLk9ETos71/fFJwOruSD+ODjSlL8vzZrwXvykAAAAAAAABH5Ayxp2yuw26w0bgAAAAFf2BH9U2gJ3XeqIHW8fsRuu5nN2RSejKH02UT5dd6oru16mdNlrLh8XyjREbjXkfrezOh7rg9ySRN5kJFQAAAAAAAAAAAAAAAAABjfZGN9kVnbVPx7wSqRl+zKri3s5u8ncsG2AAAAAAAAAAAAAAAAAAAAAABlOcfSH0v0y0LIhcLnctTewKjvDXMCw5AAAAAAAAAAAAAAAAAAAAAAAVxrlSPOfLhdL9P1Rn+J2phI0x0DoPkAe4gAAAAAAAAAAAAAAAAAAAAAPP6KX1T85+M5b72uio7MkUuqHH7HS/Eg9wAAAAAAAAAAAAAAAAAAAAAAV3YjCTX+b9oQuH0lV6IrOzNsAJVCAAAAAAAAAAAB/8QAKhAAAgICAQMCBwADAQAAAAAABAUDBgIHAAEQIBQVExYwNTZAUBESFyX/2gAIAQEAAQUC/nWFxiiUat69Zp/6ANnVMe+zXXxy9Vi5Yifw7vcyfXx5ta5JUbrg/wDpSR4yx2LXEo/Am7FFKJsxpBHJJlLJq6THqk/h3BKWtc01vNbgbMglqzWq2PGxru4pkB0flbK4rYCdo5MopK+0tEB/7j2TKJJrh56FnywbGFC6a+lcHzdoB4hYrotwZV3Vs+eLfuGfJRrSPPgVBmRFHL4X6x5NmKnX7Rl0sFfJrpuszOkD79LMiKOX6Vh+wDz5iz2C5Hv+lVqs1jJHgwFg72STGOvaxjxzsPfaIGMLFVb2aYDX6Wc8zwwVBRk8dKYXa4WEmpWj6RBMQkNg2SVKUORZS4ld4bpCFbKFwBwwyFeNYLofYSY6I9ljDMYU5zHJjNH2CvzE6z97JZhq1AdtA6bjGxs23TisLoyYqlQ6ULw2a85rFT8ADhZcIAz/AGX1y6TkSlSihznSLmrGrm1p/jY13lMHATJ9GSTGGPYNp6sCaXQ8YY+bRWYyrtUl/wC4DRlCnANsXS/R0yk5V+bm0I8cLDXfx/l5ZdFta1et6Et+56sRpgVrhMRy3VSKs96NaujoXvaLBhXVlfTT2x1HHjFHw8CBmIDWFS7l9/E9Yqeo4H0/dw/c/LZ7r4xdNWYtbFw5gMsHvFsBaoNUQSYwbUZ/4j13VezbZgIU1psPzKeijyhR82kz+KbQFvVdWvG+ZkObXVKxhWwrDr4NniUCxqrOr2DCxLOM2MKgBiaXbnlfRw15d4TjxFRRx4xR/RLYCgcFaBHZq3GU158jJyH7bWyQ0NpzYDrJm95rUXIetXIrJragxcQQ7oVkHV+BC5HGdjJ/mazxx4wx+Pp4vUdpx4iolSQJJHJJjFHbrHLZWdKq/wAvhfoWSsjWWB7XDK7LxQPZ3fSvgFrV3dzPmIn1uJ6iz9pJMpZBh8yyAhcQQ6/8V1cebQjyzrvKCt6MrLy+OcVSDWCnrOz+nsay9R8NdVjGXp+jJJjFHYGRFusGercui7X/AENTsfCxfj+rfyDsYLkCZrBJ8cvlA/LeNFsLcCzoYq6brFX6VNy7vPfHlUS+wpfp2SksHNpHgwFg8z7WpWyfPqLgLAZnB5s1+LUBHXxK8Nz4ePWTwsX4/q+THCw8MMhXjFzkW+xJlmCZXymGQgWblhew11coCJttgggwFgvD32RJrtL7m7/TZsYVAD+3H2aRTrhkfi61/CkXatnzxb+Z7MRZGy2cAP1M2Y0n6fPr3ld2PP1n7yR4zR59CqfY1tjXNRrpc8nsmv6vknGLYigdrfQJyTR27utZSkNrWbU6vHWguXNzlYX9RSewo/09gPsmTas1yOsp59lt5ov/AGbgXVarDXBvFxYQEUbjZZpMgFecWWQDVo+ODIqq1mNtfmrLKOPKXOSwc6vSOvSNyTjmvYdDulnqY1liK1y7HkrOt8R+vNqQZ4uKeX62seF4e+yJNbIvXM/IqyrAT/qMC/QAVBd1e2Tny8q8mTcNPE82OYdxPTWlgkT0VWrj5eWRKuv8W1xk26VigQqJDE83qMExWWUSDLkccK8aCeIqLvtiTHrJrIjOaueF3e+9u6mk9hS+VwoZeR6K3sqxlHtNb1jg2GjmiDMhPG+hsKfCGraqgwyn8rPf4VEggTW4n12hhp/CceIqKCvrBcvHaZ+UK3Vv4/3sGu2DNyjVYpFXe6tsk9e14l90d+E8vwIF2z5M2ILAZlBx1SVjydvrhmB1zwyizQWMquTj7XiylWtBW43ltD7Bq6PHok8blesppKvQ53HBBIQBvp7BWZsq5rixxLZSCYhIT9oT4NPK7vPe3lTSewpPC3v5a6qjhVWzOYBzUyFuzjx+im5qm+XZmnCcRNtaLiYUzQ2nOxyMCx/HaH2DV/2Dwv1vzylpNJ9Z9Q86JaFX7MJZI+WHWxMROFScSE0+g+1y+N+tPs4et6/1NP8AG+2afo2Z1nLFdVNgZBdGFKTOY2uuGIOK21t0MqvZa8viOxtj7PzagMUJutivUVjx2h9g1f8AYO92eeyJaXV/mA3waH4q11FtZL2fyeLejdTX2slTfhmQsBvOw2EauBLQTLe8XgxLAe5hkK8Za0FbjbCgzhtOsI58E1o130l4jsh9YJr9uBsGLJQG4iuVMxrkIps4OdT2D6GLbEeXWPXQGYNb8dofYNX/AGDvswvpO+qIMQFd7NbEuS5wOQCpdlsegqPXZWQ9n87nTMXsdYs5NVNDMhYDeNh2GEsxGFaW5jXq8NXAvDarLr8TU3L/AFvNwFTrb8uSiFwnjX6xL3PLBSD0nVBsQtb0z9qvyqyUZairvEuxGMTFedKzvPjtD7Bq6THqk77HF9PZa99g7W2nYWTjJQYnl4PPmLOVc3ReFPxs0p3lc6Zi9jQWc+qzgXNOwj+YlXCranDjZbNWC9GlzcWDlf1uUdkCvGWD+NoP9ysGsWcQJR+xEwWdkcAuDUdkOQSioertfVl+S5JYKKC8yZVxvV5Rr3MVBHT1T4KwUQ1HBX2nsbZVc1TfLwsCOGwrqzW4q0F3KVBHSeF6aDr0PM48ouutcFxUX0LDTwbDjPqthjL/AMta8/5a14LqkfGRXX16XzZPAFHRSgYGkCiSGyERRQ9UyE1+QHrddCseVVhWZkuzZ4OLW4biJ26xSYH05E76uKa0r8iXYTFZzOSs3bo210zA6q7i3QcU7FWH9I5MZY/pSSYwx2p5JZnNX19EDxyiDfQIqLKhsH6Le1LEnVntSXrNI8fWCQHXro3lZr+NbXX+tNHbD5AfcwJNTy0u9TtC5I8ZcG2uFp+TRIzqpbWayOANcdRcFfG1MVN8nWvmKziO7MUnITK7d+NtYTR9RWLioF1y5h2D6RguJwiCqg139Iu4JguHbWi6cZXBw+6rNctjuLNcqQeQDxCxdrngxjc8TVaR7HW6CMgL4wL9ABVCmri1djCsQRB2wpC1W3EcwbNTxfAFVmnR1+9HJMYik17As9NKQZ1jYmY3ICIiov1GlOsU5I49lEiLY2ADg93eDRQNLC0k+X7C6IA1gyI4BrBaPwBWIrj8b3bM1E5xWJpHKDmfKg5epMoqrqwX/c/hE+Ao9i2LCetTjlPixBIQBuVuv41wO7pVkywxSeg5WtjYkdbBr8ZhEKxcVAtHscM7kcmMsf7JDleJKVdEgcnz+h5/1JVyTbGPSRhstoZFXrwcpNs/tFwBSssVhqwYGMftffxPVHJjIBpOXiuyPRGGtMxh6+QaQr7XBJgRaKdXS5IbBrXPDJLaGVbIGfobnC61lPBxc5aVcms3EayZfqGGQrxjtqxdOHbCcm8iAeWTguuXZEguqjcpP+TcF1UFjGJrZMNwWpJw43tIWu+kmrGXSSPVjLrInWYJ1na+/idByY4CN7WwYM65Yi3ufdW9waM+YDxRy8udgZ13iCXFkG41oaNIK7cV/Oe+YtYtfx5Z2r9R+p98Uq9ZLxOApwVn1r7+J6/LnjiV65H9MoVxJl/HY7AgKe/OFcdCsISCNpdm7qf/ANk+nY2N6pJaWo51Eruwi2ePaK3rHdLhHYsKHWVrpQqri5LJ/EsP2DWk+ENi46sISCPXHxS+zijnN2sceUslNpuKOPt8vKuSa/R5x/8AL1XAddKA5/49Rm9utfLTVDrMzpVZJUuOXy1+6EUGo+kx/nP/AIqa4d77ZOqYCh1f3gv+fskXEeyNLYP8lU+yH2SUifAWAkgu3Pli6FSB/P2mJ/uBT6yvsiqvVbrXCdmuvgCayS/AE/oWNd1bI9ds/Q2DK/g5Onx0tksQAMS0L+ff4ipK6uu7hd1Im9RPHJlFnrUH1Ng/oEQYFQGQ4jF9teWBep/o7BtPoYOwoc50gOtWpPFC/wBpWfz3NGXuzvkJFy0V/Ous9WA/5I/p3BdCxr2uxcYKx+1//8QAPhEAAAQDBAQLBwMEAwAAAAAAAQIDBAAFERIhMUEQE1FhBhQiMnGBkaGxwfAVICNAQtHhMDVSJDM08RZy4v/aAAgBAwEBPwH5UiSigCYhREA0EacUkKiprjKU7KhT7/J20ZbLiKoJgdK63ty7euJxLEkyFfMv7Ru713aKCIV0sJk4lx7SI9IZDDI7Gepioojygxu84mbD2g2FuBrMP2BmBwIY4Grs/X4PtCvWDhE2foIa8H0WSYupmbDL77egO+H7kHbgypS2QyAIqNKRLRE8ldFPgGjg4ch1zs1QqVQPC/7w74P8RIosupQoc3aOzo0/8iFs3K3ZJ2aZjfDWY+22CrY/92nbsHtx95dYES2hgiCq/LVNSEROitqTDUIVWBO7PKCNBW5a43w0+GodLIIOawUTDDMTmJbOOOiWOJSkQONpCJ9uIdlfIYZKNlUAUahyB3UjhS3cKIlUTHkBiHnoRRO4UKkmFRGJsoSWMiypIamG83rf4dMJ2LYazm50xpEum7JusRszQ5wgFRG+/t7KxMuEJEV1GopAoQLsc+wcIWMQ6gmTLZDZjTSAiUahoMQxKWgpX3C1cOLWVYcLaglqG6gruAOaEigu4MccoeOTJCBSQwqcxzmh8ehAIGcELYKBdLqbHWaJskwskKF+8fV9IYcIF2oapf4hN+PrcMORROuYWwCBRwrDYifB5pxpYKrHwDZ6z7IVVOscVFBqI6AEQGoRj7srljOZEsa6ypsp4bYncnI4sKawCAUKX90KlAhxKU1Q26D80YYB8GJj9MMErJLY5xL/AO2I74fpGtazKGadhIN8FHjDmuRfekXFQd6x2NAKFb9sTF8eYODLGwy3B7hWDs4WiImEOgYMUSDZMFB9yVPkJaU64ltK4BsDf684UnbhygdB1ygN3DpMFQpDFQAKKZsQhQ3HFATJzQgKByQhJQzOpDluipnqgfwCHSuqToGIw1R1KdM/0mTZCTtSu3BbSx+aHh+Ryw6SFm4qa10qVNMMaU8w84nTpJ49Mqjh400spO8f3pFu2jh66I9myiXf5qts2wPxf3hDNJBwHwGIATaen2MPlD7g+wVUtgbV7oTIws2ru2HqjQS/CC/dCrVJblDCSREQoWExsuzhpL/Uua5F9wXYpnEqpboByiP1Rj7jNHjDlNHaIBHCZcVH2qyIAd98HUOpzxroas1nqmqQLUY4pLZGFp2OsV2evPqCH0+dvOSA2C7AhqoRFciihagA4Qvwtu/p0u37B94dPFnqmtXNUdDkBFEwFhocDIhuhNYitQLlBECJmE5cR0OltUnvGGqWqToOOhsiZwsVIgVEYdSVm6Gwj8FX+I59G3pCu++H0heNg+KnaLtC/wDIQeXkHmjSNQ4QCiRqhDdfXBfiGmS/uCPTE6/cFunRL2CsxWBFPrHYEPZklKU+IS3H6jes/DLcIiYbRsfdNWyNIar64tDYwZgmYahdCaRUi2S6DGAoVGEgF0rrTYBowjgwmiu95d/JEQ6boWfNF1jy2Y4ANxvuOQ5V7YV9rSTlkNrEu38h4QvMJRMEjqOU9WoADeGfreHXDRyBOSccYTCjw1NMl/cEemJ1+4LdOgqgyeTlMnz1c/W7Dp0cWWAgKWBsjnS6GMsIrJ1nBy8rEB3F9D7rhuNdcjjCLoioX3DGsIGcKPEiYXwCSzoaq3F2QUoFCgaHx7KVnbEqf+z9WuUbwhQrDhLVZiei2ZRz9bQ64cvSSREibU1bIiAlHEeu7qGlBgzaVT0o6odWqOW379V+6JxwdfMr7HIDZhDZi9fOStUwtDj3Vx7umFm6rc1hYolHfolq5GztNZTABiaLJuHiiqQ8kR0PZmDtmi2s3kz7tEjnLpE5GlLZRw2h1+uqBoPJGJxM0mpjskW5Q6QDPZT3VGqSo1EI4ijBEE0+aGgxikCphg75MLi3xZUdmC2WgBCjc5BEoFEYly5WyoHLUD5DCc6ZzBMEpqS/+Qeqh1dkBweQ1pFU1gFIdo39VMe6PazuUuTN7QnIH8vLds3QkvK5iprEjCgt2V8h7hGHSzlqWxM0QVT/AJAHiH+oNJmMwC3LFb/4j6r4w4bqtVNUsWg+7KHzaXW3By2lPp84cTJ05X4wc/KDDd0Q7eLPTgouNRpT3ReGAb0xj2gQB5RRjj1eaQRi27VwCzBWVobSxqwRMifNDQ5IUig0GsMtWY1BLfDRIF3CaQ/UIB2jHCU5eO6ogAAFAA84l8uXfqACZalreMTeWoSzkWhE43hsAu/aMS/hCLNlqRCpi4dH4gibOcgDlD4K1e/z6r4VcgYAYzwlNhw++XqoBD3g84bhrG/xCbsfXRGH6eEHeolwvjjah+YnGtdjgT12xR4bMAgW7nn274Omuc9uxfBeaES3/NQ/7F8Yeyg8ymiwlMFAp080KXQUV5U5A31l6/AYdTNd6kCbihhAcc+iFUjomsKBQbu8KwvPVlGpWpQwAL8/XfviXLLuWAmfACifSFc61qIBd274Zt6cuTubv4GvD7h6viZMSvGqi7lKwoQK1AQv9bw6P0zABgoMERTT5oe+yPq3SR9hg8Ynztdm8ORAbNugiOY5Y7Lo4oqu3Vfqj/6ERviTqsxqm+KWyF+HKrsuxiZgwlqh3yxbShsAHcFLvMYZTdNrrCqNymKYa02bsBu3R7blaqepURMUuwML8cBCEHUiZqcYSAwmDCHszcPFDjaECm+mo0+U4T/FFB0H1l/PnATdwDUWd1ilO+vfElaJNERmrvAOb68O2Hz1V+uKyv8AoNnzCRmy8sbOnVaJDlfgNL92HhBmyMxmIIs7iDTqDPHOOEL0qioMkbiJ+P4w7fmGqPGVyIiNLQ0rEvkyzdusyXEBIfAQ2+qQWXtpGko7SC8C55j+bsIMYTmExsR+YYtVXjgqKOPhvgAEAoMTUzNRAUHalkB33wcClMIFGofMN3KzU1tA1kYk85couilWOJymuvGvWET1YVpgqI5Xdn63/8QAMxEAAgEDAgQDBwQCAwEAAAAAAQIDAAQREiEQEzFBBRQgIjAyM0BRYRVCgfAjUiRxkdH/2gAIAQIBAT8B+lLAdeBfmXQUdB9HhppSGOD2q3mJPKk6j0ywrMMNUgltTgNUMvKfVUUolGce/u3McqsKe6aQ6IaiTlrpJ4TbXCEcLsEKJF6io7vmEKo34+U1uXkOaeLy0ocdPVbQG4fT2qS6htzohTP5qdY57fzCDBqG3aX2v29zUl8sB5dsBir7EkUc3c1GhkYIO9X4jSQIgxjhMk7H/GdqkDq2H61YsgbB68GYINRqAGaQzt/FHONqmt5HUvI3SobQsofVg0oIG/pBB6ehgLW1Kd8Va2/mJNParqNba1Ma9zUzG2tUjXq1WFokwLyV4kFjVI1rw2PVKXPapH5jlz34pAFcyHrUtor+0uxpNQUa+tOTdvoX4RQAUYHuJppITnTtVtcFcjGc0DkZPBPiFeJHM/8AFeFdX/ivEptbiMdq8T2kVfxXhsy6OV3rxCXmTkfban/4tpp/c3qudfLwneoYxEmkegyxjYtWc9PRPE02F7ULZEYMm3FTpYGvEYyziVehqJfIxGV/iPSmLMdTd6miW/xJG2/2rC+Hwn/c1ZQ86XfoKvZhNMSOg91I7XDlEPsjrR5GNKLk1bI0cYVuMlxHF1Nc6eb5a4FSFk+KTf8AFRXcoGMZomXNRCTO9Q3k0A0jpU07znU5qUarFG+3Fh5S0x+5vQLJJY9UDZNG0nX9hojGx9EjaELVZLiLV96AA6cHkWMZatc1z8Gy1Haxx79TTgspApbD/Y0kaxjC8LQhZ1LVexsLg7dalgeHGvvUlzJKgRug4WMHOl36DeryfnykjpwdgiljSXMibt7QqC9XOY2waj8UdRhxmjc2twczJg1dW3lztup43Hymq2+UvCWUQrqNRwtOebNXT0rjUNXSry25DZX4TSeJSqMHepZnnbU/BVLnSvWpitlDyV+I9eAGdhV9qSPFLG6qJof/ACl5FzsRhqiguFdUjOQavrRn9qMdKkYGwXJ7/wD3jcfKarf5S8MeYuMHovDWucZqWYi4VB6bW6XTyJ/hq4snhOV3Fcpz0WorCaTc7D80ZobIaYd2+9Mxc6m68PDY9cur7V4hbPcTMoXOaaK58O2nX2fvSRm5Yl+9LJcWZz1Arw7xWFwRI3tGpJrSGFrgbdv7/elKyuMqeEyl4yoqFSkYVuEcOiRn+/C5t0YGTpwt4WfEjP6YryaEaVO1fqM9SXMsuztwVWc4UUnh0pGX2rXHYoeW2pjUV1HIoZmA/Ga8RjlljON0o28kR1QGvNtggr7VchJ0D9D+KZZohhvaWkVHOYW0mhcSxbTLSurjK+meN5sKOlLCiLoApI1jGF9IsEYZWUV+mORlWFfp2n45AK0WMPxHVTeIaRpgXAqSV5fjOeFo7tGMrivEOYkeQ21OdKk1Zg8vUe9SyrENzUEzTb42qW05kmqiZLf2G9paCY/yWx/io7tG2bY+8AzsKTw6d+oxXkIk2klrk2I6yH+/xWfD1+5oXVp8vRsailt4k0a8in3Y4NTfLb/qo7gQwL/e9HTOn4pIVjbK0GDDIpbVQ+s1MqpLiLY1I3a4T+ahl5bhUbIPu1YqQwqSeWX429cgyhFWsayRgtviuYquIlq4EnWPrUPNmAjXYCpLcvjDYxXlpwdQbemS5kGlqjhWMDbf6Sy21J9jXl018zvVy5kbkR1HGIl0j6hg6zOid61tDDqk61aR4HMbqfqHbQpapbhWdZF6iua9yRGfvQGNvqJXWNSzcIOYG1IM0Om/1DIrjDCri3RkyoxirZdMQ99//8QATRAAAgECAwQECQcJBgYCAwAAAQIDBBEAEiETMUFRBSJhcRAUICMygZGh0TBCUrGywfAzNEBQcnSSwuEVJENTgvE1YpOio9JUg2Ok4v/aAAgBAQAGPwL9XT1ZtnAtGp+c3D8csdKzSM7yHJcljrfNe/P9YjYV0JYtkCOcjE9gOvhi6MjPVh85J+0Rp7vtdmK+ouMkjrGBxuov/MP1I9D0fUPDDF1JWQZWL31136fHCuPGujy5v1gUD27OO/34NPUhIK4ahV9GQdnw/A+RZHUOjCxVhcEYlqejG2sQu3ix9Mdi/S4/1wUp6iamZGOaLhm3G6nTAWVIKmwPXZbMTw3ae7DO7F3Y3LMbknFSmYZxUElb6gZV+B/UlXNNEdhNMzpKPROYk2vz7MV9J0osdTGgUXy2LXzb/YN2FEcjmM+cgm3EevmPhja2EdRGcssYPHmOw/HyC9PPHUIDbNEwYX8uWqq2FG8YzGqUa7rdb6XDt5eFXRijqbhlNiDijhnirJKUyBX29OT1SdSWtf3/AKb0g6MUdaeQhlNiDlOGoZWtDU+jc6B/67v4fA0XR1quoB/KEeaHPv8Aq134qa2snkkopBZdqd73+aOA37vu08AihjSGNdyItgMVYfRoVM6HkVHwuPXirhB828GcjtDC32j5FRF5w0qvlePMLunzT32IPwxHNEc0cih1PMHEcTSIskl8iFtWtvt5LUSKY6ekdl1PpvuJ+H9cB5FFFFf/AB/Stx6vxtjYzdaNtY5hucfjhiSFpConhICcGYa/Vm/Q44mkRZJL5ELatbfb5PpL92k+ycRzRHLJGwdTyIw0TkQUhP5CPjyueP1absXN4qKM+cl/lHb9WI4YhljjUIo5AeR0kXYKPF3FyeJFhiUsoJWnYqSNxuvkUlWtgZ0KsAvFePv92HpKaUCMm6llzGPnl78f25VzSOVJWPPcmQ2tck8OH+3k+MJSQLUXJ2ojGa536+CWlmA6w6rEXyNwbFL415gxSjM4FwYzoSOel/kzLPKkMQ3vI2UDDR9FMIaZTpMUu7+3cPf9WBLBJ0rNE2542kIOGjnkeqUNZ4aq5Yc9d4On9MQ1lPfZSjTMLHkfBJUVEgihjF2Y4anpTJDSOciQR+lJw1tvvfduwrigNmF+tIgPsvjK2eKSNlM0AfRxyO8bj6r4WSNg6MLqym4I8MCJIBQzVAjEDRj0Cbd9+O/f5CNMryyS32cacbczw4e3BFLTQ0ylbXfrsDzG4e7GWqrJJEIsUHVU+oaeCnpTKIds4TOwvbCUtKmWNeJ3seZ8mLouJv8A8k1j/Cv3/wAOJukHAvUHJHprlG/XtP2fBJUVEgihjF2Y4aHopSmv5zIPqX2b/ZhpZpHmlbe7tcnBSngkqHAvliUsbYfZZ6aXTaQyr6XHUH8a4NSsJgKuY2Qm+vf6/LieWCOV4jdGdQSh7OXyTO7BEUXZmNgBj+z6ScNRIAXaM3Ejb943gaevArelIQ8jDzdNILhRzYc+zh37sU9cqEyxPs2ZR8w8/X9rFdTZfycgkzX35hb+X34mrKi+yiGuUXJ4DD9E0wNDKXzxGRrrMBfRrbufHdiWprDDNUnSIx3OQcd/P8b/AARFVCl6dSxA3m7D7sdGfu0f2R4KrUZ5xsFBG/Nv918T1jAEUyWXXXM39A3t8jLV00c4sQM66rffY8MLs0mpbf5Um/8AiviEpVPNt2bIhjtlUczftHDwijnv47Amrb9ou7N38/IabqNUNpDE3zj/AEw21L7NmMlRMoGl/j+N2FRFCIosFUWAHgemqU2kL2zLcjjfhgbChhDBs4ZhnYHvOuK7/R9tcTV7gXqDkj01yjfr2n7Pyn9n+MJ45lzbLj+Ozy4ejEPUh85J+0d3u+1ikikQyQqTI4tcWGuvZew9fgM9VMkEQ4ud/YOZwYej6sSFplEiZSCV1PEcwMdJTFfNOyIrcyL3+0MUfR6kXJ27i2vJf5vZim6bkn+ns4VHetyfbp3eB4qWB61kbKXzBYz3HW/sxHU7DxfLEI8ufNxJ+/HR6OpR1p4wysLEHKPBS0Ct1Yl2j2f5x3XHYPtYgLAh6gmcgnnu9wHlLQQRvI0SrGkYO8kZiezf/wBuGXNtKmWxlfh3D242lEEoaj/lHUfTdbh3jt34UsHp5o283Mvot3HiNfjhZuotSuk0Sn0T/XwTVdRfZRi5yi5xnSItNMcscKm+Ucvv9uFpYSX1zO5+c3Ps8kxTRpNG29HW4OFRFCIosFUWAHyS+M1MNPm9HauFv7cFKergqHAvlikDG2Ia1Lttqqw2u8Kxy+4Hy5pkiklnqHLCJbyMOz1D6sVVRUU0lNGIdn55CpJJB0/h+rwSQBg1PSExpZba6ZvePd4A5IInmaRbcPm/y4qxGJJMr7BEOpuNLAd9/biCmQkpCgjUtvsBbHSDoASU2fW5MQp+vwQUyEB5nEYLbrk28JYyFRV1IRWZdVUnKunYLYVEUIiiyqosAPK2+zTb5cm0y9bLyvy8JimjSaNt6OtwcOlFAIQ5u2pJPrOGd2CIouWY2AGNhTF3o1bLDGq6u3O3E8v98GSoVPHpfSI+Yv0b/j3foKLMzxSRX2cicL8xx4ezCJVBMr+hIjXDbr9vHwB6aprdiTbbPUMq/Xr6sLDW1hrZr3znh2X46318iumiOWWOB3U8iFwkmbLsI2ktb0vm/wA3hZ3Yu7G5Zjck4igiF5JWCKO04gpkJKQoI1Lb7AWxSy9RZZKnxhuWhzn6vBEVUkJUKWIG4WYff4IM1ilONuQSRu3e8jwTJoZqoGFFPI+kfZ7yMS17p5qBcqNr6Z5eq/tHyn9lU7lZHGadlPzfoev6u/H9q1UYYX/uwJ4je1vq/wBv0JndgiKLlmNgBjZ02eaPNs6aIcufr36/djqVgavuD1haK3Ec/X7sVfQ9ZCYup4yt7dg9d9P4fJ6T/dpPsnFR+7N9pfDPTOQXhcxsV3XBtiXpOQdSHzcX7RGvuP8A3dngof8AX9hvBNSVF9lKLHKbHC0qVfjcuXM/Uy5OQ3nX+mJK0tc1baAHcFuPbfN7vA+zbNTQeaisdDzb1n3WxDTkeePnJf2j+Ler5SSaPItJNlvNm9Cy2Om/hiOGIZY41CKOQHyGSeujD3IKp1yCOdt2Pz7/AMT/AAwJqWZJ4zxQ7uw8vkJqR5JIklFi0Rs1sGGlU9Y3eR9Wbv8AAHyjOAQGtqB+APJ6T/dpPsnEoZgC9OwUE7zdT93gkqKiQRQxi7McF0iCTVThVQXIUWtr6hqcU9HGbiJbFvpHifb4KKaokEUQLAu24XUj7/AaqYF9ciIvzm5dm7CJPI8rytnnl5Lx4achw3YjhiXLHGoRV5AYk2bZaqfzcVjqObeocedseMPrDR2kP7XzfqJ9X6JNV1F9lGLnKLnHiyKYqZyoWlj6xY9p468O7CSVDJQxng+r7vo/E4lqpulh1R1VMFs7cF9LFXCD5t4M5HaGFvtH5DPV1EcAsSM7am3IccZaSCSsN/SPm1t9fuxIsKQUwJ6rBczKPXp7sfn3/iT4Yipuk1SRXYL4yLIV7W4cuXr8hkdQ6MLMrC4Iw2Qgz0r6Fhoy/wBQffjbQ1UfobRkdgGjHHMOGDSUpKdHqe4ynmezs/AatqozHWTCwUn0E03jn/TtwvjNTDT5vR2rhb+3wNWdGLtTM2aWEtazcwT+PuEQlqaPq9WGZdLX3hWwqsZq6fgoGi7he25Ru1wy5tpVTWMz8NNwHt8BSC8sMR2MCprmPEjvPutiGBhad/OTftHhv4aDTl+iNSRufFaU5cvAvxP3ertxL0vUx7auWBpQjabMZb5e/n+LlFWmgY/4kaG49pIx/jV0i+pE09i7sXNpa2Qecl/lHZ9flXq5wr2uIl1dvV6t+7GXo9RRxA+kwDO33Dh8ceMLHJMHIBqZ205bzvtbhi9bWSO5A0gAUA8d97+7DwJRQVlRc+ay7Qg9rNe27d7sOI5fEoTuSHfv+lv9lsKiKXdjYKouScHJD3Fjj0Yx22wCzZ14rYYfqZCvbfC5m2FUnozhb6cjzGAqQx1ItfNFKLf91sJU9K2kkBDLTqbqP2ufdu04+CkmK+aeDIG7Qxv9oY6Oky5LRbO1/o9X7vJk2b5aqfzcVjqOber67YavlW8FL6NxoZP6b/4fLFHUVkcNRa9m3DvO4fK1NTlz7GNpMt7XsL4h295kUmeYsb3tzvvube3wf8No/wDoL8PK2lZUJAp3X3nuG878NFQr4lDuz75CNfZ6vbjbuDDDIcxqKje1+IG87+7tx14RWzEdZ6gXHqXcPr7fBJJS5xI7CMyp/hg8fu9fgDUtHJIhFxIeqp9Z0xHVVjCorEN1Vfyacj2n8duGMKXjOo1GmLFAvaWx5yUDXcowdRHEgzMzH3nCywyLLG250NwfI6MTMM4EhK31A6vwOGR2usU7InYLA/WT5L7N81LB5uKx0PNvXz5WxDTkWmbzk2vzz8NB6vLqa/o9BNDIdo0Knrg727/67seLWz06Mc1NMLZTfW3EHf7d2Fz01Ur26wUKQD7cB2qXgY/MeJrj2XGI6inkEsMgurD5GoVzZpWRE7Tmv9QOOkZiPOIqID2G9/sjy5KWjUVFYhszN+TXmO0/jswxBeqkHpyyN1YwT7uOgxDPUf3quWzZj6CN/wAo+PLh5BimjSWNt6OLg4DRdH0ysrZw2yFwb33+VSUi3AncsxDcF4e/3YqP3lvsr5FVVQ1MJilbMNu7Zu70fZ2WxT0SOZREPTbiSbn6/IqJI3yTyWijOu893G1zgVDjzFHaQ/tfN+Pq8mSTI0mRS2SMXY9gwBW08cdExteO5ZOV+f47sCalmSeM8UO7sPLwNPMskdQ5BaWJ9TYW43GGamtXwgXumj7vo/C+GR1KOpysrbwcSS02RtouVkkvl79D+L4An6OeOLi0cuY+ywx4xSTCaK9rjn5cH7yv2WxUvlGc1BBa2pGVfifKNF0XMUjU9epjNix5KeXbx7t61FZnpqJlzKRbO/dy9eI6enjEUMYsqj5SQx+lTtt7cwAb+439WJ6GrmSGnk84jyGwDcfaPq7cGWeVIYhveRsoGGWjgimoQwALqQ7jjx+7y5Mj5qaDzUVjoebb7anjytiGnYefbzkv7R+Gg9XkrUQxpJI8gjGfcN5+7Fky9D9KNZViA8xMbcPom/4Jwk7JNRNwkQ3U9lxpw3YC1cMdYLekPNtf6vdhEjqNlM26KYZTvt3H1Hw7Osp0nUbr7x3HeN2HNDmpJ8vUGYshPbfXBWVXiAYJUwHXMvx5H7sRTxHNFIodTzB8qD95X7LYn/eW+yvky9FUhyxr1Z5QfS/5R2c/Z3p0h0gn933xQN/idp7Pr7t/yc1VOcsUS5j8MTNSiRTEQGWVbHXd9/gaTooCaBjpCWsyes7x78bAdHT57kXZbL/Fuwtb0jkkqRYxxDURnme3yvE6dmWtnW+ddNmt9/1jH9pSqDT05sl+Mnd2Xv328pujhFC1NBYssiBtoSL+rfbTXtxT9J9HiSooJUzNxaEj0gdN2h1/B8W6UeSeK/Un9Jl7+Y9+I5EiEFwCJKOyhh9XHfh5KdkroxwTR930fgcbMTOyp1DT1NyBbS1t49WCKyNqFuB/KKfYPuxSHxqdhNUXaBWJQKT1uryAv3eCiqlFpJ1ZX7ctrHv19wwkeTLsJXjvf0vnX/7vKg/eV+y2J/3lvsr5D7NstTP5uOx1HNvV9dsF6hXFDF6RGmdvo3/Hv8moq3taJC1i1rngPXishrHQyDzsYVLWXiPVp26+XVUZteVLLc6Bvmn22xnmQ9QtDPGtr24+8e7EdRTyCWGQXVh8ht5+tI2kcIOrn4duMjylppTnlnYXyjn93sxDSwi0cS5R29vf5ElRUSCKGMXZjjxikmE0V7XHPFSziyyqjp2jLb6wcVBdQKd5rxmxuTax9Wg9+GqeiVtIWu1LcBf9PLu/2wUS+yz+epZOJ3HuPwxlRthU8YJDqdPm8+ONnWU6TqN1947jvGEqoKgy07vs8kg6ym3PjuOC9PIYZCLZ00Ydx4Yan6VkmnXMNnNbMV55je9sdGSZTkBkUtwBOW31HCM971EhmylbWG4fVf1+VB+8r9lsT/vLfZXyEhWQsIIQCnBWOv1ZcUKxD8pGJWPEswv+O7wolZUiF3FwuUsberAihrqaaRtyJKpJwtKLZ6p7WI+aup9+X24hRQLTo8bX5Wzfy/IGqpQE6QUdwlHI9vb+A0Uqu1KWtNTnQqeY5HEdRTyCWGQXVh5WzoSlfUc1Pm103349w7d2N8lVKTZpX9GMG51+iN+nsxsYOtI2ssxGrn4dnk0fR4vYDbvpoeC/ze3HSv8A9X8+EqaWPPV0/wA0b3Tl2/74ljmR5aSXUhN6tzGI6inkEsMgurDFMlCEmI1eo2Vm7Fudban3YZ0U1lIBfbRru53HDvwsNaDWwX9Mt50evjx+ONmJXKKysQvVkjb8XHLFbUxbaWdcmVpX9HrAHdbn4IvH5xPSMcr3QDKPpdUcMQ1UJmO1rMw+kI827uy+7yoP3lfstipTMM4qCSt9QMq2+o+Q8mbNt41ktbd83+XHRv7tH9keFJo5dhVxrlBOqsOAPLXjgR1lO8DHdfce47jv8Ec0RyyRsHU8iMBH6QkABv5q0Z9q2xT1Jed6FjZzVyEqU33AOvcR8fLNVSgJ0go7hKOR7e38CSNVzRlvO00txY8e48MZlro4TYErOdmR7d/qx/xOj/66/HAd+kYGBNvNNtD7Fvi1Ksla9riwyL7Tr7sLT59mr9TY0qkZzu7ze+7Cy9JXpKcj0AfOnl3evXTdgQ0sKQRDgg39p5nyq6fMjLtSqtHuKjQe4Y6QjnZIomjWQzSPlC2Nv58ZVkkqzcg7BNB6za/qxJUU1C9LIzattNH7cttD68ZqWS8et4X1Q+rniWqpqzxvpFbyTUZU57XOoPzv64pkeSd3ZFcrUG5jOUdUcgOWNqn9zqeLxqLNrvYcTvxtyHRV3VVMxsPXvG+2uHpOmYR0jRynrW8241HLu/rhf7E6Q/vK3Zoqv0vm6abgOYB1O/Ecy/3yLL51ol/Jn4duIK7ZbbZX6mbLe4I+/CJHUbGZt0U4ynfa3I+o+S1LMSmuZHHzW59uGhV9tI7Znly2vyH47fID1FJBUOBbNLGGNvJmimQSvUgxRxm2/wCl6vrt4LOpU2BsRwOoxMklLCa+Fg4kcXJXgRfdY8rfN+RzOuwqv8+MC50t1ufDBENXTSRcGkzKfZY4/OKP+Nv/AFx+cUf8bf8Arg+M9ISSpbdFGEN/fhvE6VIWbe+9u651tp5f97q44Ta+QnrW3ejvx/daTxmJlkVJnjbZN1SLgnjyvxtgpEBoLkswVVHaToMZUnE5v6SKQttN17Hnw4YMVHFny2zuTZU78NTTs087NmNSvVYcgN/49VhNrJCpzLUw36uul/ondgR9JxeML/nRAB/ZuPDljaUdQk6jfbeO8bxikd1BSaoWAsz5QgN+t7sVC0pjgqkNmNKw6h3WKbvq3Y26AzQxnMKin3r2kbxu7u3ASoPj8PKU9f8Ai+N8XZhT1jEC5tHNfh2Nu7cFqa1dCBe6aN/D8L4anz7RU6uxqgTk+8bt2AtTegmJtZ9V/i+NsK6MHRhcMpuCPk2d2CIouzMbADHmQWhQ7KBFuc2u+3M/DC1PSQSonK/m5F0Tv+kcCKrizZb5HBsyd2Fq4q3NSBW6hFnN/mnhyN+zd+hMtTUjbAfkU6z/ANPXgDo+kQRD51TqT6gdOPPDRioqqk5LNFTg2K9qr34BaBKVCuYNO/usLm+DSrMZyzmRnItr3erFNLRx7eBYsuTaAZWvv156ezH5j/5k+OJI4ap4GNs3i0+h9amx34i6Ork2s75stQthwvYj2693fhkdQ6MLFWGhGHkpy9DKeCapv+j8DhWfPFb0KqAnLqNwbnvxSeORzVFPO2eELCLsQD9EX3XxNAkRir4ntVZlsb3OX/bv8DySU+ynbfLCcp337vdgvAPH4ecQ6/8AD8L4SPP41SjTYy8B2Hhu7uzCCREepG6OXqS8dLjeN50OC/R1QJUtfZT6N7dx92DF5ymN7mCUdR9d9vVvHtwIvzas18wxvfuPH/f5KenckJMhjJXfYi2JDTB3kffLKbtbl+hLtOkYWzf5R2n2b4Io6F36ujztls3cL6evBgMxVJTlFPTLa/C3M35YVplSijOU3lPWsewcew2wrTK9bKMpvKercdg4dhvgRQxpDEu5I1sB4akVvjGw2zmn2pJTLf5vqt4FNPXUKzMSPF5ZCsnstgVb1D1VSt8htkVbjlz3+3wVNTlz7GNpMt7XsL4hnWpd5B1pmdtNlfrC3r3eGeocEpChkIXfYC+I6/bJHTOobPIwAXsP1Yaajl20atkJyka+vvxF0iHSOTNs2TJrJ6wOzj/uXp6SeoQG2aKMsL42T/32m4JI2q6fNPDhgxkBnseo9hNFu1Hu7MSTIDPQX0m4rfg3x3Y8X6VZ5o9Ak4F2X9rn37+/AlhkSaNtzo1wf0WpqainNU9yzSrIpzdw3+q2Figi6VhiXckayADC+M1PSVPm9HavIt/bgRrXuVH+YqufaRfEr09T0jUG922DPZb9g3Y87S1k0qr6dVcacrv34RqmWGkQ3zC+d19Q09+Eaplmq3F8wvkRvUNffjJSU0dOLAHIuptzPHykoFoqepWRBIxqesp1OmX1DBlWmhpb70guF95PgV69pXZpC0TTG7FNPvvv8FeUYqbKLg8C4BxXVOb8nGI8tt+Y3/l9/gkmlOWONS7HkBiWlooJlaZSjSS2GUd2t7i+KXooVTpCzEhXYlE0JJC+3EdPTxiKGMWVR4JqdJjMjzGQFhYgWAt7sTV1VG8csf8AjU6jOb2UX58MU9bGZEglGaCqjOUkEX4eibHd378LT9K2jkJCrUKLL/q5fVrwx450QUjYrnESnzcvHqnh9W7dgxecpje5glHUfXfb1bx7cLFXL4lNuz74zu9nr9uFdGDowuGU3BH6UYp66mhlG9JJVBGAj9IRsSL+avIPat8fn3/hf4Y/N6z+Bf8A2w2Towsl+qWnsSP4cTRRJDSo+geO+0Ud99/bbGaqnmraV9JEkfMR2rfjha6lr4YK6JcuyqJBHmG/Lrx137vujeWCGopiw2scsKyXXsvuws9BBDDFOoe8UeTMOH1+Gu/0fbXHSn/1fz4iSWeOJ5TZFdgCx7OfgpxSQRtVbZQZWsCqWbjyud2IpejKp2q4+s20OXMf+Xlr/viM9IwmGsUlHvbrWO/T8erwyU/Rcc09TJ5yeILojHjf19wv7HPS0tUUhfYrQzG8RAG+x0I5W+jja9EnOnGnkbUa/NPLv5ccCLM5hRrPSS7uNx/y4SnrI1SfhFPob6eg3f6zywZOjJfGF/yZSA/t3HjywUjaSnIN3p5l0O7ep9WuGiET09Si5zGdRa/A+z2/oslRUSCKGMXZjgijoXfq6PO2WzdwvcevBCzJSoVylYEt67m5viE5KyuTXZySklO3rHThjK8MdMLXzyyi3d1b4PjNZBElt8QLm/ux/wAV/wD1/wD+sHxisnle++IBBbu1w+0Saqvu2slsv8NsFE6OgIvfzq7Q+1r4zZPFagLlWWHTcNLjcf6YbJU0rJfqliwJHswueppVS/WKliQPZino49ViW1+Z4n2+Gu/0fbXHSjdGJBJUAxErUXsRZ93buxSz1MccctE91hyEAMDrfjwHswxk6KkoqfJmWZ30blbQX7/I6SpY06tEyptL+kdb6dhHgklWNFkktncLq1t1/BDNTRU0lI/UJkvmD69vLH9lnog9Jz3Z45Vm2bRLbnbde510ucZuj2FZET6LEK6/cfxpgwJUT0pQW2Eg0W+votuwY+leiaasHzMjGMrz11PLdikKqSEDliBuGQj7/wBFnotrsdrbr5b2sQfuwxrJHr24D8mo9hvf14HitJDAwXJnRBmI7TvPy1d/o+2uOmKejXNXSQB4bmw0uP5hiZulXNVWyknaI56naOZ464io4Wdo472Mm/U3+/wW6MqkpakMDmdbgjluP4GJqTpCiCVeTLHL6JB1GfiG15aY6RarkKlwmRVUktbN8eOFjgkelUt1IaW4Y8td5Ov9Mf8Ax+l3j/ZO/wBxK+88MMr1lVHNaxSpu1vU2KWOtKTLAxa1su0/at7NLYNUnQkKVbrkd4JTGhHYliBuGFz01Ur21ChSAfbingpkkjp0OdtqBdm4ff7cTTVlPtpVnKA52GmVeR7cO9HTCJ3Fi2Ysbev9S9Jfu0n2Thlc2aWBkTtNwfqB8CtWSFS4JRVUktb/AH446W6RkyjxmfcvBtSfV1/BV1td0lDDTalG1bIg3AjQDTCoil3Y2VVFyTgVVUA/SDDvEQ5Dt7fwfB/w2j/6C/DDKKQxki2ZZXuPfj8vWfxr/wCuBKwmqrbknYZfcB+qKLbI6sJDEVtqGIK/WfBBeohp6GJbD5z3O82t3DfwxX7dqmNYGAjt1Y6gdbU8+Hd4PEaObNRJ6ZXdI3fxH45Yi6VqvyzreCMH0VI3ntt+OX6uqpeo0sdT4wvLU5x9fkCmp3K1lQNGQ6xrz9e728seO1Co1FA1sja7RuXduP6wzqTeeFZGvz1X+XECQ9In+0nhjHUY7TMCua54cd+/E8ksEMNJGoW6hrmT4b/aMSTSnLHGpdjyAwN5lnfLGp1Ea8u4fE4hpKcHZRCwzG5/WFDU5vychjy235hf+X34nWRnjrIpwS6f5fLl9LE+wrZJKKQaU8i7m061/wCn1Yi6MjPWm85J+yDp7/s9uJek5B1pvNx/sg6+/wCz2/rGspVuXdLoAbXYaj3jCws1o6ldnq9hm3jv5f6sQUEEb1CvJsjOp6tzutz1/BxK0R2u0k2UA3dW9l37r/fiGlhFo4lyjt7f1hI1K7ps2zShCbsliCNO/wB2Pzo1KXuUqevf17/fiSXIkedi2SMWVewYV0Yo6m4ZTYg4M5D5aeMsGG7MdLH1E+z9YyQyjNHIpRhzBxPEkomSNyokXc4B3+Gop6rzEkzC05PVPJTy3nX9Yno2mZDPMpEx3mNeXefxvHhKU8ElQ4F8sSljbAM5hpFzWIZszW5i2nvxT0m2efZLlzvvP6warleeOVgA2yYWPtBx+Y/+V/jhoeu1O2sMrfOHxGK6sIcZVESn5pvqfqX2/rSr21/Mo06EHcyg4hdSbzu8jX53y/y/pf8A/8QAKhABAAICAgEEAQUAAwEBAAAAAREhMUEAUWEQIHGBkTBAUKHB0eHw8bH/2gAIAQEAAT8h/jppHaloRJJti4PLOpXTTVrSZbz2/wAjD9K0nBDWTB66AOnNTZpzTDwb8NhyCnx/t/CIxhSFQYiAVE9iOBCcAsd1RNGT7cu0KQMnNbDMvJv9Ima04KkTZzuYwI1Z/QwwBbk9hIaw7MqLNeOLHrJYmTIgMUCQzN8TNacFart4E3OtgAp0/wBD1/CFImesQdGXauoeKSlU1rCpKECROb5UvMTCcQ3RjtVME2XpqIUWckT0LifZLGpCMok3Z+fe3bgIWAGmh0BE+oZrRgrEdPKLpfVgcErLSDVfvYP12CQR08+oVx4ywQ6KoPQWyAqyuUi612SiOOZ0Epu6irQsL9SZ2UVOZYCsvF4zaso4kz8C2uQ1ajLMU5o/L7GtKTom31DQbMnlP3zSFIw+HlDt1EaYMsGY9qwIpYDJMQQm1vgP/wCJZDkiMpIxZJe+ODuCx/nTenkRSKUqZIE+QWfPf7Ol26iNMGWDMfqHqfvmMqRhrJyVrAwQV2ljqQYPMbIGWv8Ab/2Ohtu8awIC/B7Gyc8Ikh8qgfPJIm0GKk6YU+32UQV0VC7tQn1+gEyNImbKjYRMpCryX+Wdk8wBhGzVH2M26D9yVJlln59FXxmcRUWWfN2YeG10ynMPRek/CfpxbFjJTBK1lD74czQFAmWJiVBGGSYS2YqOYYSsjxoB1EJAD2BDIS8PQ2NMRQTsROqqT00VuWH+ugLXndT0VNtlo9C0lEVISBPaKPhvlGUD15okthv5OHJ9NCsRMnqPLtGAGXoej4V7ITpQZi1MJRt6N8ymhDc2fBSsbxx2GjQIZJjTO00dehnOUgSwUZdHltC+Z7EbNxtr/wBEAHt+l/74h+wnY4t14IQ7IyCmKwevTWwhYf66DfGZ2kmBZmGnbs0b5EzCajEEr4DksakIwmDVn55o4MeCgeM5phQk8jyVBgBqkkDRv59+ONTKGU5UY6P0kn+ihWq6ORFY48uHoKymYIybrOWxnyfL0MVmoYVvEwAC4ld8Snn18ERref6W9BIicKAHaodXaHDCfVScxbxGqYzh5CJ6IImQwuUYoM2PSCrVBuj2wD4D1HIOWKRUww3k1J9csIXECYEMJEvp9esNmSJD5BUWRg5uYZnfv0+IzzGIoqDynwJvGPWFCORYUWPkWysmUPWeVnLrlq4DLjRIpyjgX1i/Aqj7YRcDNaMFQBo9O/fzIAtDkOS/Sp5Qk4xp9LrR+CEOyMgoisHr9T86fhnOJi9ouI9+qJ+OrjTmmH4cLvQqii4kpE+G/SdWGJKFgzANFscQ7zopqCYyzFTm4ykiCmxGaPy/PCQCpSATdiGaeGN/SD5pV/V5Tk9IVhQBi29atG0YiS+X5sTpg4in3QIETT6XpLOYwjShJOvzAVwhAgRjS3uV+PdJOKCVhxQC9SXpsV5LcmB6JXln4BYJyp0wRq4UbkLjZ3wGmJWIgEnmBk5PK1luhu4BJnZKj6OEc4BUADtUOruONQQeASguDKaMq4oG7ewDDSgDxtlfZE7KKjMkjWTgZrRgqANH6W2703mJXk/PJS8IZ2gcWcg73O6YGemXBBr3weI9C0NSgRjGnIu0AGssWEp+e3oWWhkVtW1q6gR28MsCZAiD9r+uL2hOwLFNMwGfJ40KsJQBPmuBcmgV+UIcefRIXaQACfF+qCDvoDIKWBm4zvgT7RQqANHu/Gzsk91rj1idlFRmSRrJwJPmEhiUWDry9vEzWjBWq6OCBvWCqnkip4OBZSkr77RGrJLGXuD+xlOlRkLEwkWnRL49NLimUaRAsPSCzFEZRRaCGbRzDlzslCZrJQrvBHsrm3YyiMNZOfNV5Yr43+vVM/p4WVdvJOGkxKQEvl4kLsJQBPmuTRSywyCZcKPr0kidQLo9Eg+U9HEpHArCMbGxiB+HgZeqIIgkYJWTEnPJfi2Q1+DFw4xfqERQD0ptuaLjDIuDTfJCRF7CRLpYy/ZJmtGCtV0ceDK4DtFCKTcDLHB22AwA0shlNsQzxJckEgZVEmErTSjM/ojpbQqwlIUeK5sMbs3NOorIfYPB08QBERHsQequeajpyJeSpcMU5TwNk6AsdfkNfL0UvfihaoUdCZHIw9e/FJlKApm2/wBQgpeoYBJMTrFWWXFt3jWBAX4P0Mm1IzQBZfMb69P+wqUlAwZUJTZP6G6FTclJhpw+F41izWL1ICjQV9qvFTVyWEKD0/0HuHQVeoNUO2E/A+mityw/10Ba8iJXiABUTBJBpYOYJnhNsasSlian0zTOIBpdEi8G/RbM7s4pLSzL1tgZDAjkcqRSKLIo5dIwpgQH45sSaiQqhE6GFyUaABOWsczn8Dv9o4RzgFYAO1Q++GsvgDqiWrwqkk8tBHo7DCssQCXXYDoxOIKLW/FW4OQ1ajLMU5o/L+hKjwCNrYvBd8RiwslhKijIYIRu+7adqlJiVJil/p6f57BUFaMi7FAX2En2mBUibOEyBiJRSksHek8jgE8GNblJQ7x5hOWy5OC0aGvs3BxnpW5ajQklmYDBhzZ96bzEryfn0nvEHSWsBC5JpxTHFs8yEXKCgsbDvzzwod24CNhQ28ZFYC3lB6JXln4DjX98rQAUXAmRy3ZysETkKGg23+0x/YmnM5BZe+ScQgUU32QplKKsN8RmghK7XxFm+f8AACT9JvjKbeY2QMoP6f7ZOg9qhl8mkR0ZEolvmDHhEvMiCYAonFCdC0bZjQ2iPjkpGHWEkQdMYeaMXEiLkiszFlJM5cwx3RJYndgaGMW8DtacFQBt4kEhn/kT/vjkEpgpP74f8mwp+w4jSSSaH/48d8Vnl7KTZF0/KKgvRAesk/UXxKcWGUw4tOnyBj0hiWUW4IzR+XmoKzM2f3f79uxBqJiqESGzC5+zEZYyQw7CKve0WSugwvaQYUbOz9X8z+TQnWOGQMALN4il2y+f0Cx3AC5vBhwQmCuW8M2YZnExOwlcDr+MwlyFJDlfA6qAuqdBJJmyegxKTBnTjWL1CYZj0XBQCEMMRJnQzT1xgWLBNEghmWicSHETbaA+GuBArrQfieKh3iVJ8v8AxywEAwAFp/vMJHr1qT2FFTJYQFOn+h5FYuQGN+bV9+3fQ3ETVCjLTI5Ny7coJMpQoptv3th0qnigOKk0zMSnMArTrArpUslMniS9HaoBjJ5g+OTmgly7cfwPNHLlp/jpMj+jDNkQtKPi1fXJas02nYjGfw+8g+LBdkIhEQUTmR4hulRf0RKjzQVyOWSGx8mS5WE9iI20Xemx88cQ2W54Qkhx1URHuq0qoIFzYpfk/Xsl2Dh+MhNEQMg0OFy0kIWEjRKg6259mFPGltRwBJ7D45ah8i15HI7GLn2+O0AJMG10cdpqHNtrVNwDs7djUpKBgyoSmz0W7wkGARQgMGvnkiUMGALKbcwJMeY45wYECYRNI8L3KyVyQCy4fLvltkQJ1UKS/PAauoyIMiMI/Okd/oWwqJ0oEBej+x92ORbS0MeD4cMsWqJZiDMNyLIiZk1sIWH+u13+o1tN1rF1IhPwb5TYXp5AWNBljHPEWxYyUwStZQ++VK9lYo4TcT0k17/mdwixWVlIkdeVMbPUSZShRTbftuDtu4RQiajkz4jhVLzKb0LoES6gZ4qEwC2Zmit59WcBnVuZZUGAVANX35M3boByNQI38+oPPCYvFhyQmG+OEaCzJ5rYprMOzYDBkaJhgzJvMmaoNmkKRhsp/Ut28F+DVxcWBjZZoDx43jPof9T9QTYQWSXoTtaDanDweMAlKlIY2mvj0piveXMkgB3ORmF43LOQczfGqZvXIESL6SZMJaIoiRWE9uC/i6Jl2wEWWyQS6VIkNZGDQsrRh91UkpkZSkxEFll0iHXsSNNESHBCvh4XqslmW2TvIzYxdQOfQ1hSoWNqTRccvxHAx2WVZIhFquneoruRJh6RgnHHDNLKircg29KzrhUNaT1CUEfATqfSkvrAPBNIy66OfODBZ+rX68/p27e9rrJMVjTZhc57DrqI3FMqYOoPtkknY0jL2oD55E+PLhE2ITGdluvcPiSjALSLgC+uL0xnJSgDiQGm6TDzRW5af46Rsf0HKXIV/nG9PKgkU6oRyg+gUYVyNc84DDKgCTKu1fZorcsP9dAWvAauoyIMiNj86R3yWLIkaU/FL641bDgjF6aBGznNY4gmM2CDLJi6iBVPcwQBRZPmNiZiOGMJBxrPhI8qsK4dyAmLxYckEw3wKExYjSlDwkVnPLjsdAkezCJIqTC8nXGIZW5IKSBc+Dg0SOFAhPb/AEPHPApTYPmkIeh8v6Vu25inTIKgxKrHjrkBiXghCWC8weAa9Vdw+chAwfPT1yb0VWxLAM4OKJGpkyDDgaL7fTVKRkDB5kf3+hfPkYL8HT6NQ8RuRrC0g4YhNxDpNFblp/jpGx9ztBk8lCNXSzUlcp8virljGUD4GuGSqQD/ADi408qr7EKlgIkZicyR40zrgXfZIky0VKGwncZU5niapMUChdDjBdQ62ELT/HSa5lUPIgmODIpEab5KkAuVW2oIemL1zdRzQqTbE0bvAON7looqkRMOxXGJMA6ebxAAybT4j0j/AOCuzhZOUXJJmEUeleuYULSqYA699sZuZbAEnT/Q9ez4cuKlPO337D1siFtmi8mHbTUZYIDNYsOCEw16U/fMZUjDWTljLkSPgUXjHB/iVCBKSYiBE7j3r5+jBfg6fRqHil6oYBiG4xl+Rg5FlZHB1iI3JP69E1UOm0+xFZxxLD5EnmxGDF0jH0AtPbNMG3AWhquNE8Q1oYQhu/AhM86x4SUBJmALbfc08SZNzjLMgzueDEWxS4zVo3rzyF83dB5AHSnkg5cVK75WlhGcrPLRS3VoWyNCyGuq5TfTZ2mt8KHhmBkvXpCQoKFaXxi88y52TR8iRu5g4Z7I4gGkfDJuJ4I6GsCVUYkxSrpxqzRbjkgZNAwaQ/6kkaiWTNZ1xcVP/egwsP8ATmPI7ToDkZIEb9qgbt5CENLRPOmEuN/y9YLoNK2rfsgD0hmUSmLfz7XYeBo1dVZD4Jn0Ks4YTAJ8IifPD8NAnMsyoUZZP6LzMkOajwkeVQJfCNCLjq5AFzt9cWJUwcRdpWuaj75NpgyMMZTgVMe9elxcNMEEono09cc+A92iIlnCLo48C6PolQEoW5QynIyupUYIUZJQiMLZrMAEEFiV+aJWGBjhypl1Vg0CG5mVnXDHKDgZH/4iWljnwA6WW66FUDfHuAFxehLkhiS+RV/UFSmMHw+eSsdJEKBMCSyEyvPF1/GUymQhJcK4o09kGcXcphQgjgM2FVJgbVCdR0vCccVzAFlNuYEmPMc2g96SCCyCloLrjOEL7qoQCjEoBPieBmtOCsR2fppPtECtV0c3hw8e14aiYhccx/oG1zsQR4Gc08y0pJBIlfimRgkY5CC9IsghsH2MGf2UgWGGyBBDJkihfMkJhaZpEHtKscQV2QYLBBMFTYcrIHA6tAWwhhnkfDoNQFWggbbn45fArmYsGFGx88HoXa/Rmrk1jyq+W2TvMUsDBR9M8J2hOEZE2csBjs7LLxUAFVyVIhB0poRCEU51zvcIqRGwROSUmJ453YMvYaATCEpcvPLyrdimSyyqb+OU5e+DGbuV8BLHDsIyZigjPSAWJ47FwMZiTFQgTLDwxWUjGFBJFnMCvngsT2wFHQjXNBXG7meIG4wprNUgn9Jq1aQEKPN8wOOEHXAATfneCP2O2gibHuEM74yeNRyM3MNJvGeHU5ALTXJSyUVxjm5UhrNEzA/pdxqQAYzRFZ/pdzO7CInMsBWX1aJK+UkXrKpih9CooZcCd8kXJP8ATEQgdEQbKwhbHSSfT8z+TQnWOeKC5y6ULQECjUSehVq0lIUea5mcekKkmBGxOa5I3WiyDEA0ORPowIiU3FDfEFQ1xKHpDOkhmz88DHDEuGQC2HSEqol4EiiBbdqDuSqHZyMlWLoBp1SUmFjkHSYmCfDF5ZbVM7KKjMMJ5P2r9Kgkm56QVZEEa5LQi85lgKyvNj/puJiV5PzyU2sNpZty+3ggDDYkvCKYKKrllQGhLgxvCe+bixLIjYmvi/XPu6QhGhFfM+uT4VIgcNittlt9wpF4TR6blpn63zOIIpkrIQM4IKIPSWuqARDLaYekRUeiVOuEwJ8Io/PMa/y0E+P/AFHpTt81hSsHg5b+WWYkBoGo5HNIYEgImPGVyTPNbCFh/rtd+kmbHJAXb2r4OR5FWRaCwEXJGCk5EL6gBhlNAy1BBPLbsCNZvadnwAnjtb5JIwZImLEf2HgsS2wFHQjXqK4v6Ysl0zmavQC+A32nBWI7P3UJkbliSRZxywHSInwCaxn1LxDJVeEhSUPiX55M8ES704Uqjsjh2AiY3YqHWHDpKzylG5GAFQcpJ3xPaTksneoVISwmSuFhHHlktAHG3ftuDHGplRAcrMdnoFSZGSTlAGE/HK95gwSRYXhQuykXF4lQGhGTT9TaFPWyAwFkbcDFZ8jXCTor8uppSaBlLprOpSxmglMHae0uADV6oYEkqVa3kccYZYFRpqMZgQJBpz4AdLBVdi6BvgnSLLpdBQohjfJ3xVnqnbmUh0mJ/a6K3LD/AF0Ba8vfKjkZ6hpN4zyekZDtSAtkTByDGKjPFYa5yRniE/chTJ2JvqKzwEw8zdISubn69DDMdeJki+bniWB4PROO2dzjiQ/bafYisY5gegEIQOEK6YBIcCTq6VCk4fEvzwJejtUIRl8SfPE3u1zfErEpYmp9lzmICGGwSzSUM8jt/wAELCtpJJ+l8CnC0VYgGIksJK8nsh2JJL6ZBEk3OfShm6iNEmWDE+mF0zdlEApFfDOuJvBoVAkAUClhVtNR4RLzKAEXIq8SvqZXhZEzMxvzwj7iV1zDYjIYueSXOYFyeiQfKftf+j1pZP7cnmElgK1oHpeN8nlYQUa8IZb/AF7qhIKFeRm5YEhi0zyTPQOJmxiaZJSdZWx06C5FMBvp6VONHfipgzMw6b5Ut0G4AgoThR1hNejq3BBA0yC+CwOMQlAnQEEDBx/wD/8AXNR0KHRJMgYd8OLj/wDeIUayBrowokMkLdzxHLZGEY/0AJ3xAejvUAxk8wfHFZ6jODQmAHa9KOYqoTQkQG1wN/HbRApg+MwdH8N+jqwIXC/Fr69Ek+tWMhBA0yC+avqZqMJ8EW4fQ1Qq0ZLEQkrOZzK8DNqYFQBt5VbkZIyb7fQqX2CyS7CaeyUSeRPSlS4KilIykjGGS2v4jV4LjSIxEWffo6IfLkqDKRay3HJzWxyoILIigmzfpOUHrt/BqNLLduSy6MPupKjQ3f8AHoopZa4BcOBP37HCahQKT2mw+wzxeNi0mSHSUmmiGWP49IySKgYPED7njbBDDUktdoZRPMMM1mJhWNkZNueU/fNYUrB4ONxf7ZGqGGVY/s4KTmAKyq9qr9/yGJf5aCfH/qOOrdRmQRB2i+RDVNxOhF6jukoUkzDmgjpxU2bM0ycNhHfi5p2Ysk/kQ3GoJIjlqEB8dcsg9UeTG0if/hMsM4CIXXKFrEnFLLnUNAT8ji1jkSp5wmG0AStr2v8AIR3RJPJEKLNQuuG4nVRSRlU3AP8A955GOAMwaDRwO1owWEdPKATj1KRtIf4f5G27xpAhK8PM9bvgAbac59ZUz1ZErQippbKQT/IbJ0EJpiFfIasPSWNSEYTBqz88+KJFiQmvhGNc/wDYnZ0GA0AfyEfrwShAxYgCuvT+OFHJ5SSpGHGmCTmKIeyIcWndXyP5Q1E4ylUfJkfnuHjdJQqBg8QP7/d//9oADAMBAAIAAwAAABDzzzzzzzzzzzzzzzzzzzzznz7zzzzzzzzzzzzzzzzzyD3z/wA89ft/88888888/wDPPPPPt3vvvPPPPNPfPPPPNT64XxuPfJuP/PdPrbfPPP8Az3zzxj7gPfXzfzyU17P7z73zzzzT/wC8+04t+RXw88888t/88882888vf886XCUbWdf88888888841+884tnykIF+h883888Vw388849/wDOa/vvTfL3tnvPdfLHvL10/PPPPOZ//PPPO073/Pv/ADmfn/51Xzzz/wAs8188888bZ6xQ+8pY9cjZQ/8APHvPvDPf/PPPPTffO+uT8/dtl4GvfPPPPbfvPPPO8t2Vn/E+P/vdu3OvPfPPPPPPPPPOdamhH31fcr7tfAurfPPPPPPPPPLTvPPPPPPj97jbuLjfPPPPPPPPPPPPPPPPPPF3fXfPPPPPPPPPPPPPPPPPPPPPPK9fPPPPPPPPPPPPPPPPPPPPPPPOvdvfPPPPPPPPPPPPPPPPPPPPPPOXbvfPPPPPPPPPPPPPPPPPPPPPPPefPPPPPPPPPPPPPP/EACkRAQABAwIEBwEBAQEAAAAAAAERACExQVEQYXGBIJGhscHR8OFAMPH/2gAIAQMBAT8Q/wAsNAJQUJmJTEwxOzwAvLZvhO5PMPb/ABzmQyEpgOIQtyEWCCUvlZTe4zeFkjURqBQ4Sxl0Jx58WduxJeDc0eZDzhaDrYJar4MKaxMkimGgzgEYktpEn8ouEJFzaYvYh1i5z/4qGfEEWYDshK7N6KlF8WxyWUtixzlAOdbACAsTGsZ9LVApWdOmPd86UKWktMMDbuDybhOakpOJzkOkEhyMIjTjssgKqWOC7OcNgvxAeIukpyoYu3lnMQQUUJA7G0kCHQAC2kSB4V9x0N2hpssBp9e+9OZ4ST5/dFRZZZu6ftcF6Gky029/IooTOD970OJBNM5EmP3XTgeTJnqYwaQXuahRKiSIhsw2g1IxpTRHDabPVjHLOq8JuzgP2AyuAu2o/UBTrPmgjUIJArUVxgyXibTGJtNPpBgheCWCQZgDaJq/XBMJskZC6TGSZqWUtpMG0t3q8Z4w8uF/cAkiSOEnI6OPBDoldkGnlR2ErYKLi4OMbc9/OiPkduph9KsZKTNJilY+abKnHY/sUZuhHGBWAJmK6uxPc3giOspEKwdBZkjcIsQUwArISTpZcNi6xEq1GkKH0MMJm2XvAjLI8Uq6rwBJCYaVUufCISu3CHMNxgSgiX0vRxn87SDCZOdvKh/EwCQTe4PBIxs0BI1Wiw9L/FM2bDoffwUJXl+B90KPIDvU2c3fXpVj3/R5tzkeIQysZWQkHNuoF1CzTNsmBsYOurGVXwH2DCMPRCkqAyJCdTwSxCy4HKd9ILxJIJIZrwYBJESIItiO8SPCTGpVhJH92c0yrUX95Heo2zpttUqzKRO3b1tUkhPU/r6FFjbB+5fVACyu/Xb7/wCRz4E+kxC9huKXsF6Zj2AVjmkTiWzatveTEWIn6m8BjHGJi9N2MK9jTEFHOIdkkI5yGkwUSWm0UtgkeZzXowGZcQDzBbbQQcqCRmTDJ8lntUXLnEAazgf17UjiXU19yoTx7tELhPWB++Kky97+t+h4IOUrJt89vKsCfOPehAI28BMFnOigttiWmOIAaSBPkhpg2moaWGJV9+BdGpsG66H7NECNM3ByWcGLzteiT7Fbc3L6HKoUnLuBx/GzhtT2XN3b0KNkaGwbBocMuEf++lQS7HlTrZ1ev1RUfieBPGwffaocLrv124S7CATE99KHtAGxzXMNFQRYUQsZyFtbWOaFKy/Jn6fNaIwbEfDPvRuCMx+/HH0X2aas/kHCwTqmNR58jV2JR5BG2VnWHXc4wE0fpKurdVyr4WEzhjrUpYOT5+OtT0eQiO21RZt79eCd4CulrP3m9jgoJcU4SAkOYSbsLG2S5NIeMKe5qGgELBhsvN9D7y2Bzv3Be6wlJgYEEhIArEdIKTStUuD99UkPaL9w+ePovs1+XkcLNbzuBLJ0h1JOCoxikkgyjEIa7UUyUqyQLC4lAxmInbwBxDk3/u+/XMYZqj8frUglPmVZVPY+8UCba/Pd7FBzgOFuZUdi/wBedORxWzrKyJlEdPiooTBmFhFpLxYwbBZqJZJrIwzYgLwQm0wVBhmPC6WJcxqaJjCGXMpZ0IzrEx0pRFu3AIUqAjdBgOK5XICdycnMtwsOGWL23jlUv6UKRNjRvnE3jPCcYGFJDACA5fngBZiBiTMbC7DtZFE26Ojqa2qDeyKAgZEN5FepNvDkG3LVynzp+THfL5vCEQOdJwr5fvipM3ntN+m1Q0jWG/SP72p3hvQcQbThnNzGZiowGOduo3DHRT5D+grLCsQXdEMlqWtLEZpBdSoi9ySLUQ3SEG4tm5yMjUdD5o5pYeYjnT8JnLjle2N2DvFLlHo+44TZJHR8IKUMDADkunk2trQBhkDHRDB6uqs1BhRKAsS3gJb58MOIfuXzRAA8v5Qj+zsNaLHPPrL6etXQP0+/KKEgzgUzpvJ860cMrzNrctP2KwAeSB80/GjCMljvAkRY6zSNwBgAcy9JbS7EpSpOqBBIQuFTBCFmdJdS4XtlLJYtslZnlFLyyARskmDAoRUELohUmFjEcX3FCUtewLtGluN0OYZ6y7Uiof8AmoJcVlT0H3Fa6Te/180XQxz/APFEL2H01OLh0/W9KIGGuML+zvRQEB1DE61+Vto6C8uQfKFgnETqU1ALkSC5hkO09aulAGQRmZEENrRocow8CjkA8xGMmGGoAQCd5heGLGhGATKas59CgDIgG2JS5ZWoFGc2QO29rEu6hduhgiTFkXEMSTc3/wCanGkPesBnv558dpJkx0DVtxGWgIBsaDLMrihwgERTQ4NsXVBuRatIxFUpFiZG4zsWUoIEe0bhYixYUlNMwzoKwF2hZCYhbetUJhYgpCWcwxcsWsUDlihlJSNWNbS2zmnS6m6BtFh8j/IDDPZik/3AAAjIzlbhneSGidzwtdxIbrbayoAaRa7YNNAcie7K3f8AQ5ECMDQCRlEomzk3B9NGEgCbDC+4sA3ow8ODCCMd/Vzf6BAYkmCaG2EyLIhUYuQjS2dpAkO6biy7ShgRpTxSiq6rl/0W1GzOwyto94C6URSUM786K7FxCxEgupMWhO8NRnhQYiSbMOJzGn+hW1SJNrMehTlJGJBWCRsjnRJ1hMDCh0AeqL3/AO3/xAApEQEAAQMBBwQDAQEAAAAAAAABEQAhMUFRYXGBkaHwELHB0SBA4fEw/9oACAECAQE/EP1UwQL6ZSS9YZ72/TkXu7n85Uu+287+kgx6wQxh1KI2zi/xUUyaQBEbf+ILj8tFbxo2d3X62cXtSDIdVqCZoGWPoxEL7+FXK1nYbePrG4bQtW3J6bTpj8ho4nd97KgkhZUXjlfjjZVt2bxEZjdLKJGkztL1ELwMT0HlltQJIxe93s8Vy9aSkgLx50rIQgoy8BMaufbX0aQmxh6x8lOXs3zSo7sPx6IVgKNiAt43e/CrzXpOJrBBLAWt06xQHTXxp1M0BHLtxPqgkPpISmPwBCBWply774pprAJXzWrvIGc6OkbOlOVAK8Eudw5UHkgwGNjpQHwE8sf2gwI93+TT60L6ufKbbjy1N9Gx5vKASEZipLRkdvmnWjBwHogkNY/HODbn32Xq4tKbd6YAh2ehE3ae9NC6D5aZBr/VFhGfF+j5pwmAe79U+XUpbTj97qTQe9r3tyoRPe8Ne1uL+U8NVRbZQrJrvfwQgDxKAJUn4JCxm7Xd58VfLHc9Zx0RoZpIueamKM3pjZ5l5GaYZks7XWoMwkK5vHXZDV2M7PyQBnexyg3bH46vaaY2bByym5ZTdH/KOpke/wDDXPBGafEz8PxWUj2n1tF7YZ841lb4f7boNLX3Yl9h81I3F8+atV+lZTbrUqxGiY9mrlWw0OBTPZXaU+vW5bez/Du/gKxSRg/zgyb9pEvIJ9qRIQn4bghaHbE9rVhUekz4Kmwd958c2rwN4/VOXCmav+39v1UT4PTDJP8AnRoNckm+dntR7MaP7szTDbC2wj0g3Qq0zcg5+w1MqyxwNebfhHogEBRPMTT64MbrUCY6Tw2PC9RVvbMdbPxUmCWkbdRPZpxXoH4dPv17auy9GnJNrWkzQ3fXvrvAEGPxLBZE8Nag7uB+J7m6o/kNWZ57fffUxpexw9D4SsFTHnK3eWObr6IgJWlGSMw9886y4Jf6GprHSvpC/j71L0gQ6T5o8qhTwIgy/wCUD3sOSOxt05evbV23oleCeb88PSRRSaTelQthN74fijIVhdP5sdOGIsT4S/Xy9CoReDVte86fcUtv2dHm45tP0lZfS68GebYpzghg3HS5rQgueBDHRtwb7KEeICJg5X5kyUJexJk+udtJpcwk32QB7cL1C6EY4llCxe+duTUICbvTKKlFJcPR5i2jv6TwkZ2PLznRJcoVR4Lptn8eRCbxW/OlFyE2YOhB6TjLsL0KYO/P13oEnAnZadnHbUoV2FuM37HOg5GRNULyuoZg3N4s7sNj/bPPrTLYGwtznHeiOBPE79u+ofHWj5O4VOlsV9nXvSuFtPI9qmdJ+M4o1/FYqHM68aZHBM/iBQdf9n4KZxusfNIWAz4pRIle6Y5RBuutWZ97pjrNLSvF9Y9InhFkSHlpTZBTERCjv14Rtlac/QXoVeRVL8U9YmLG2sWAs7V3bCihMDnj/aVcJ2+OdqZM41WemvkLTPXsecaz/wA0QErV6Ab36mlKITJYe7PKM02KHd/r3q4Y6na53pg2czafvvTyPRHIbMdLSUhCJNlymk8q7t7UQIyzwym9Bhq5e5T6cExpxqYMl+zFM2yttPO26giPDY0iIFv03VLqNiz9PlqSXFEI283PH/mjMIycSlGR3adC3577h9qDHIBoa423qCn8AtRodlbNo23xTStUm9m/wVN3Qidu/JffVjTtJLbGRqewDmjqCNYJ/Usv4W+KUN109o7UMyOfPfpR3/d2/sQEydbZvbfn3oQuCebpjSle97P9z0/YiQmCYo00yDs8mlu2dGh/L5oADB+xhY991KLJQHPG61JQiH9i8kUpArrEcmiGa36/9v/EACoQAQEBAAIBAwMDBQEBAQAAAAERIQAxQRAgUTBQYUBxgZGhscHw4dHx/9oACAEBAAE/EPtzA59k7bODmtM65RS+cmphQZIspt9wmmCgTAbUhQWgqIem5H9tkVLaTAl4T8QRKKJBEEVpQQv2OrtE3b5Urd0t1x9o2dRAo7oIISQ9GIVJaLDvCULJP0h70uCFNRBRHEeGaET27s4qMFF9tBRY6CCiBXRPSIQsJJh3uK0IYUuvS4gU1UVVdV4fZTD3MwQBcWXb7IQUS7YuAQKpqlSPNrwUxqFbRBSOJXHL8pwCZK+oeY/PoCNlBJfgNlV6pCCdQCaoASdwfPvRgu+oiLYPaaRH0VSlxABoggiaJx55TkmRmSQQIj9afOjwBG0EETRONhL6Q/WOmoDGHoyhiQAgCAEZQBQlzNhBCsIPymAoerv3kbCDBUWGqvnmHUhoWQD5qnQUcrlWu6jkQIMe4wnqdZ7zpkFkPfVRPrupNGYEoMQS6c+eORMDdQYMNfakEolUCYTAQKpxGnTHOQzuEjwQo3LwkBNO9QMqk3gFEYGVeAQSkIIGw/ovmjkTA3UGDDX6nTr/AKTRmFQGIjNOB6aRg4fNaNAgDx44noiHUiiKogQFeeGFZgW5KwCqua+w9KMmDK+CHagNThCXQINd+TjYXS+xnBTUijfM8oBUhxiDaEIalPICBBrSVlkya+hVxR7SmVqvgFOJRa1e30WqE8FBRKvQIUaCoVIlo2CMiRaRQfSXXYSoLMFArqDt4vjBaxlVSTqrRw7DqE4o1QDHETxwcTDz5o4DRIp4fPjR5G8nIq0kF45kK7qgAVRQBQAKhwjzTYQAsJQs+FyRRVlsBDabQA4gnA7L/bCXdGQSog4HC+QgTEQRExH1SWSgTB0EFF8D7KJ8jCC4AVHwElFo1XfR8UFAtSwxY0mdzGKhofJicBWjLigdRQMIaQ74hYTDLMRWAAAA9n/U9/uX5uuT8qNlAFoEzWwUcee1d1QAKooAVIAqHCHgCusJyIjeAwDrzgbCjVgCuAHjiQgnWImKgQvVHzy+EUEVDAlwEdAKfRF5FpFG9QYiAvdtmT8kMWj2G/AfSOMFIQpoACq4By7Dw4ACUJjoLCfhQbUEJJkUcCBowekosziASUE6gBinirJNut142LpcQ1w76LBiN4OVBoBQqj0a7dSlI+I2AsAiTDXWYtIMJwfAOgjOfJTsXQer8XhD/RExAuRDbRGCSEbVOwVYoC77AUj7WAxXdiqg0E/miMcmcYZ8i3J81ISiY554WvUShKy9ehIKIXrDLKIWBn2D0xtloFNO4gyouMatDhOdW4AygQ4a9LgADQAAAwD0/gQIHAEVidfHIBjh10WVCAA1Iq+iUTqjZQBaFM1sFH0/65E6cY/dX889870Iw1prZbSYUvIFwtRJEDkETQR6T5olg80BQbRB5DIqnlUgOe1CAnDAboxoeihUB6qkJ7pAiIO9K6auJx/yeID+IehT1fBpSZgkYghTgWpMPR+Q8JJ47485GyiNEERHRPTuzIRSGxSiXoAbwfAa+KfHC1RVOntXqYsyKprUyAA47uv+RrJO9AqWHHVu6pK0D+UoUEbYS+m0wiIOMi8HDoMNQII8RC1F6DBWK/DDxFQNICmh1/aoMAOzVQrFup2wIVpAz0CvtB17yNhBqgJTEHxw96XBAGgAABgH0v6IDMdJjSdR88QkcapBakoL1U+eIk/IuyxxEAoRRPex8GtCkhYQAWAhfsMv0EwlVo5NJxr81gYMmgslIa+g2e2GE9Tu/KR8qDWq2TVhh7BUghw9SYlIMAUCoBfBxwjUGiiJhnhBRKPF9BoUgwUCLBZ4fUZfHI6pMCFLSVRxkpCANAAADAPd/wAELafuSXZfXr3kbCDVASmIPjgwUkmQIIWVBQBoe9LghTQAFVwDhTVW4DsQQIIRyDoQVs5DSojkEClX6CJ0vGSYBVjwCqsiDzskGUM3NlN9KwVh9FSW5riJUETTaUVBdsEKAQAevhrh9DhUBiIzR4/8RcRrTHbHxm09H1BaUKsqiqrqvDOm6k24AoFWcH0HJSDAFAsAvg5q8ZsxThhCuoXt9Fxo6SI58lORdp6LF1CAL7paRCaeg1xYhaU6jBBAJxZCT8ylBuljr7CfTKkeWHABTg/i5RGwkEoauioo6h+iT3pcEKaAAquAcvzRrngdjRBkQUebGISAESiDACoiV1Q4lcEKGSFR9D+LfUmJSDQVKigzwc2w/sgipJCaBvsHCHSr/DDxFFMAKOXkaXdJYCwhB4SYRNR+EInWEdGPRgSUab5A4x6RThOioG0+kjZqAxV9TKKdlGWjSArfZzDCswLclYBVXNfoCAq92A9SWSykWpweZ5iWjxQVJkFD6AQ0ncAKCIkiIPPJCwgIsiIqABVnoA063FHuwREMW3R7v4+E9BON+Cjauh9HMhXdQABVFAFAAqHIGkgcYJDEiRGAjCCARbD7zjAwPRyI7IyfmCKgNQCnBn+GRDUpBWQKCcHXxyMnfBAQIcAMoeCBlyqwBVXOf54lLA4ZXSI8DJBuswAAUYlKP0gYKxQI+PEVAaUBTMEMK9ISBAFQVuhowpsZwUFBqwHC9AnJpaJUlFAhRK5Vruo5ECDHuMJ7xxUcKBNohMVIAqcFse3NYpwBOiA4otSPaSXMoI1QUnB7OHHA3RKg6OV7DjJWEKaiCiOI8JDM7VwewQdFSg9hvNMZpOL+KJCNJrFSCMioreQPgLYG3EoykeCQHI/pT8x0mNJ1Hz6GnRn9V0OVaTtwNX2Vy0jtljZu1kfWso4FR61ovbx1JLOhayd7BUsPRlmGijxSdcKQlOKjYk14YdtKqyC/0ljC3CcYSqFQCIO5rRYSEWMpqhS+Or0gRCxqAquKRif9kY/bf/ejxxKTUOoIoKoIAAe79cFOKOPRAeEI40oGaxJLIJyBQziRiZKzu7CFgksPWG2jnbLIUIlOnCqCkSumAQE+ml+BIjxOh10IHSAu7S4oA1UUANV4+mAuEvYhqeDHy8XYAK5Oysp+ROYMQEnyQI/D8+EzhODCwQYjDbhPjdzPxDB6eTBCFqM4pIwp3lLFMA4oG2gexPsFsQmCvyKjk/gXXccihUB6rGfzBgeww/FmYrK+z/LApIFLLaRHnXWO6L1rpKBxj7walTMLcuoLCJ3fU/1IOCTrFjL08COX+FDrJOIdGj30K7wqQwK0aqGsN58Vf5rT5qwN5TlgDWlTtqLYoIeOl216kmGhRALA+iORBahcWoFJIAHg8w57jcKijwNRs7YtQlNQiCZ4Rk4yqnuVgGyEk23lXjVHC7S/GD3zMLFbq2Sj35THeudh12IH4ECqgB4DDBco6iqpGImPY+wd6/FG2UQBSLLp53YAWCaAezU9LAD2EnRr8JEpZDSKc7OBnQRySlZ0FV91ys1BLaIiW4QKhizokWSVBUUkLx82cOhTcmgkTY651ywPKCtQBK4ljQB6Fd1REYgigAIgifR65pTJdAnZoOlqDn9S7uWVUVKdErfcHJkNZtNQUk3ylSRKSKCuCIbC+gJ7QDTdomMDUEsPXoUMfAIkYBKdg8O4wkBFSoYECAAB7RhX/SKMmz2CosfaVgzx2IogTIsDxwOYYJBgzsgpAoK9QOPwU20FhaB2sT+Woe9eskuh0TPt/Mr020zFKFUOHiZHLRKAYXrCnCT5iWjxQVJlCh6CszIRaY+BWndo+ct+xBkZhBgjh+VJWZdoCI6Ijyf5g55VZkEWA8I43lmsPqBpgrsi14OpTH6wwCgUC+9iatDh2MxRAcGXb7hxpGlluUU1KiLTyafpCjGDoneExzHPaO6qqtUVQqlVVX6nRFV4JQPkZchRLoJAdEsFATCQOi6/CVBZgoFdQdvBEQIpCCL7ogqt95WboXCugpg7Tl47or2p1j3kvwQr2lyKm6olbAGga6ZjNJABAJVATYXNvh2EzdNs7QqE4CdWtaxUoaJ00PHSh1ZQWq2pQOlAelTCVqalKMVBGmciyuctURx3njEnBriCwBWADGhiIXHVzL6MwUBiCXQ+oxYtv2KEFydG8AhwH/3ebLftV/N6H1JrY+ETsgZrCAa8PTkYhCvoXgXQFPB/eYeEClAC6pKLXtrNWiKVQnFNKHVzkUU1YhT1L3VmJ04nQWPxw7E4lcvo6SVgw0E04e3Y4VC8sGwqILxrVO00Nh9DsiByAmGmuNqnoDAgKOVnJZ+YMkQsdJOZejKZW8UQUaNJxDyJEWFSt+RIJHGf4kASCRqsYF7Ti4bEOJalCjhIgctKSWvRAtLRwE7AGBNdDUQ7KK+d6H0mLGp+FEboF1ltIjz8oosSGimioBWL2AYKyaIMZ0UtICocq9Yy1kGiNVaR96DKFWTHcQbSERRKaNvlnqG4Cwo8MyFd1REYgiAARBE+hhZhAwUHdAsIEzgjpXRk5AgAQ9UUyXn5IU7MRIFCb7HMhXdUACqKAKABUONeTqUx5MjAKBgXrqlMlwK9GY9pEU0KECXb2UNTPRLj1k2JewRxEPmHEpwNVFDyMGEnlH76BQ6vzgAFjarvSpTQpRmoI0zguiqDBIjllrjNhmlDQSJ7sJUpGJ/eWh8xCMBidAj72tqAegDB5PieCQp+MmrM6gxCH0xi5V+Ahf5OxWlALYM6f0omBiGajS69QzCNAMs6NBALC1O7eR0KOoiwwF8cU6RFiDQB7lIDF4II0CqSykj7SOLE955mioMgzAQBsgPwGsgbIMt7YQOBOczIV3VERiCIABEET3UC+LakQf4SgQG+LA+iYoGQDFLOMbuCGGKV0AlAu+0SqtOvHeGWMhNHorEQM40TBiNhAZSbHWsC+zCqReIc9q7qiIxBFAKEQROL0ThuAaBT6LVBwbk5ICR0WaZpCwe5pVFJDoCiaQQA3vqbvpC1p3JeBN/ziCyxBmIaW8fyG+iFxYeEYAuC3KsUPAyoikOvvYl2ew9zEEAXFl29n9VhxN6/lZ4zK+vTI2TnTBsKyUL8y1ZZUtSsGqpjHPTr/pNGYVAYiM04HJITAJFLC00xlBFZdNAA6ErhA4veeYoqDIMwEAbID8AixkGOppWCE6zwkDOcau6QJFM2KeHhxRAhQtYMLQBhags5VxlIwAOkopXgBIEIfFrl7y1dbiCrpSVpMYKUAhlRlk80AS7Aq+570jh+cKEGJIAgRchgrSnpSiII9O9S0K5akWBQtiKnBT63ljeyoUIfczR1olgM2wiqVotYBQWK5i+a+VFKDRmsjkKVOubJYEN/OBi5tECi8KbUR0PQuiG866POMdz5EkKBOOqvHSjHcNmuBwaunZ03lEIQNTfBUvG/0FJNdrJlvClq0CypVtFA6UQ9gXUbKkKwgx7FHgLT+guIAoumEkyeoBwneKCpApOqvn23cpBCb2Gs0JGwHFIpSYIj5IdIExOFhZyBJHLIMFHT6MPdjx4lPjCgJNrd3BYbo8gSwDlhwewzj3SRlohEcxGDUYkc40W3jCKWdCyq+8/wfx+zWJUCjBJgREqUScD0JCEEgjLUSvAGkI1AuhrVaSBghlAeHWMTpub/ACkJ6HA/2ANYbTgaxSHA+Y3kBUFus1gHqZwd/wAhqSFWFvkV4VKYFKM1BSm8DuGxgGQilFErE4WYx8TqKlBLC02BNadm2Klqop4ERV84/UlU9Rh3hQ2hqSrEAaKHIwltF2FlFPghEcf4dJq+bOVu3eGTLVdgYAH5IUXBr0uAANRBETEfpnGSsIU0ABVcA4eFeagYAhu4GDitA7AfXGiAh1ZFh9/jGd0X+UlPAlayeCRlGUxTmQ/Qli3eYWrC0NtWChUsbFsZhI8RXVGzSZ7FXEbBRDIcUCRTtxhsnGhDBL2G8iwmA8UaUEFRxfgc1Yaq4grxuI+BsOlMC6VIxpzOIE6YBCHI9BbwAyYwoimIKI4jzEGoQ2lQSpEGjCL1QgTigHAzwFFCh6HoIQCClgBiEHIhXMNP7iRdHpcFFaxSFXwWDtRHVuvgn7koi9MOcmCEk/tUvzAGni0kSZUAtCiSAhD/AJYVYutCxYFdC9wrTABbu+JCXPy3R475kKyG1U/RAsihCDRAClEvh5hCZjkazlQFUpj+iftcuGL+Mpi7LGIpHfngji3j0OgkSHkwtdE348ANOOO0hzGpR1KNcafo9jCmPaQwKeHd/AmEGCosNVfPqBcGjNYkgJ+EBIcpbRSUoGl6qWohxh4Vwg2ahbF4T0f6kHBJ1ixl6eXdn1xIocTOMPWgLIoAgVBQpUL5ODfQCZNO10IXtzsaINHpnTQm948DsEIn9Muid/yVBhO4BWoAJO4PnmboFqHxU6iJDRDUnaNJrZq/YfE2ZMlGkDiCEKNYfXmcyKgIjTGpQjz3yyNhRoxBjiJ4/SmbDRKwBSAMCAwHaUQHEGCoYaq+ef6/qwJjSdR8nNf4E1F2EWUhAgAVxBXwgrcwABAhDqPhPWIA6Ia1B28+PrCaOhYLMU7OPGGJXR2JFZqnTwsFGlxHKY0aItW+1wocU0RCDWzs2HGvZIy2eLouIAbeIyMa/myxJK36GpFKTIlHoh0iHF5++jpr7WdMjfI7c796XRmFYjAVmHD9oVcBoaHVdiV6/MlO714QDREEDntHdVVWqKoVSqqr6ENI2pTFCDSBXhc6s5gCWxgilyJ0ZFoNFSyIUUCCYdUPKLqIp/IqhJiWaXSagp4AUXulZMAW7vioQ58lH5vT5owN1Xh5kuCANRBETEf1XQc0qEHFQlNEfPDLBAylIhaWmCMiL6NjyzpBVMMWIoAcrvlDO0kYPkCCag4kXwJ2WOti6B6HjmUwm0utAipMyclJ5BTQUIYjuWVsJxpzJNEBSsK+1JNsSfipg0OU35D0yssZEr2wVXRPCWOQCYee8I7tPBB0I8SYog2BFI30jaKCfRO28IJLyGT4J0EW6oCBhxngqJGCANhAq0iqPXDANq2SZEpTQ43oQmmwijvRw6Bz/IakpEgLeDOmpYYQSzwOAeMPIIKqIMYm2AEP0lzIXdQABVFAFAAqHHwlt+oI4N49DoI+ug8usNoI1CNUZgSt4p5d0zR1BwrZEjAAK4SUQqGCWWNt0IjCK0aBVpyGT5pFgiUarAiEJVqAMkdjOsL8WTa2baGKAhBwww1lVXNC9A6i5Jwi5A4eTOjRhmxFADSu082VGhTMiogXI75LiEtJbQ7TDAwPYkCI4kFXJKMS0mPBnKB3MnaAYYEpRkzpWiyN0AXzeumOJoTQtjoiGS8+OeRMDNRKMGHpnKGfcofZQO3Fsl6KGQnTZ2LI4lAG6xZZIFOAALxcSUYEUQISKjbFfF6Z6qMrtpEQYL5S0gBz5Kci7T9Kav8A6JpLjp3fE51xibyxuQdWhOw4ixiCzwGiysoFV366Vq2pKxRE1IkBePnIPbAIroAFzh8mYq44R0CDA/f0KEMxSTeBENRgSWgKlNKkIAhkio4KE34Rr6Imh2AUTgRo+6OEopAPP/zK/wCQNzKStcogqZYogxxisJLlUEu0AAdsuFcEWoZGDdGqUUcfclJgU3JoJE2OuCWOIAZL3MNnRZ6Ycr0YTsVLvfXFbLUAdTsxSBUuH2U51zTmSaBOzMOlqDwdCZjEvSCaHYBR2DZ/6PpMOjs6F5cF/ILSTwH2NcCryqKANqKAGq8VPiEhDsFFCyK77AoMgCYcIG6mmgaJnoRzZARknxZLqArJ9n2OFbCDqDDoGFI8yii85072Agh1c3TMkZx8ERs+qenjHbRK4Fx5AHCE1r2CQFx0p0tac+3sbAVmIdMIE2B6fZEIVvg5oCkDQEBIOglRKSn9GIv3BEUQUqCSEr7WtRAkHIFLnswbI0i1P+csZS4wCTE7td/1JozCsRgKzDj0Cn7m2zaHsEFUtgStCM9pEAFABA+3/vI6flrOmRvkdo9wrffo7kBcN40UUhog1AoUUhSTmYH95kUbaDAk4zA/vEijLAaBn3EIVQ05eKlJp3sZtt5BrWKJxqjvBTyyrDtW5SQUaVAs/Lo0O4Ji9KZkVyGTaxcQKjz9wVFNPtIdOXFR4PSBk4CWAQAlLGh/N9P703VKMAOI5S4oEZEEETROeNYcEXSOxArWgfcaYVmhbkJUKI7jyPZHsAEIgIiPb36r9RvVEB2RxH3HkYyG7AwyjTLHQPqSQgnWImKgQvVHzxSaMt7C+4giqsIvWQKvIsCurW7plft8eyywmxu1AhZVXg/WM8WQro8hHpC545SV3+ezIYbGh9zfaUmwl8IJkccgAkiClSSQlfa11ED9V//Z";
    const IMG_SUCCESS = "data:image/jpeg;base64,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";
    const IMG_ERROR   = "data:image/jpeg;base64,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S+gfzQIHt4suxkMeAqvozxVVhafBQV9oc03AY4RkJTISmzjT2czE6DMafIRObwApAJV/OCBWBf3lLMvpmlGPka4/XzD0tQUmvideP2P/2gAMAwEAAgADAAAAEPPPPPPPPPPPP6vPPPPPPPPPPP2fPLY/PPPPPPPPPPD/ALzzyXzzzzzzztbzxTzzzzzzzzvzzzzzzx3zitv/AM89i888J188807P2d88/v8APPPb/PPK64b/AHjwB/zznzrjn7n5dlffTzzzy3EWxI+/srY3wPz7gTknyf8Ac2d8899P1/I7Sx+e8c8u50j8+98c1Mt9vpMz08/9ue88q8883uCdhe9K88888r18818teM8OM88888888898+ep8AK988888888J8sY+ss/288888888g8888a87088888888/18+988b+888888888Y1vfzUe8888888888/wC9o/vPPPPPPP/EACgRAQABBAECBQUBAQAAAAAAAAERACExQVFh8BBxgZGhIDCxwdHh8f/aAAgBAwEBPxD7aZMDNt1p2EzML+t453b7h0EbcUIA8iYffvp4sxagiYiQTMv3CgpYE5Z7PehQwH0Bphjynq99I+wMklYOaMdGTPjl7r17/f2DQmKmDzB6+30IzNmC+Wx3ugyU9cxrxWKHVF/fwfl6QW0LSzabGpm1QCCg1TCMIgCSsWm4XSqzgXgLdeqzpUUdnUvCaScj/mTxZQlDHPSmCqxZ0dx9AVoU815rK4sd7Z+gcwr5vMT26qCvRSLAQvgcYAENLdUk9iEoQEhZ0IQhOLRUCLnKwAniYozDIBCASzHBRIvA1lQdvEMJ7aoTWbHt9C7/AJNKJH1b89viPcuDzWCotUWAuMELMsXtFBCoeMACgAGA2KrNW5lvNkJkBSRDTWBOSUqMcAyKUgks0gGJEhVi6XWRhKRZygWMgJbFNL+9jzxNSz0ETiaZYxDOOevx5VjAPFCcWYSjozGtOdZofzoJehF9w+lI8LYSAsLhbEiseuqLZXmGYOFgKbR02jECX05+bZtRdEfksoBZwuRDZOSp2RZMXaoiJLzHLVN28AL3Q9IJ5KQG5ELchsC1OUSbE4sfksibMkrl8ENl16yVzU4hudESXppauJBmsIgTC6EgFIb7fXyHUwmnoi/QRYGMYvZ63YeqXFEaKYs4um7icRDUGKL+GEFhARELIItpiid2sCww6TJIkwjDargh8GcBG+6CJ+VqGnkYyBIzI72pOrsnQkRlKImcQTYIRE2m58sWQlUgTCbDmSG9QsxsTZO/nvExZCZ7jxScQfN5TmOHytyNCQLBKV5EehMPUFXo1molQmJEOsTGUBR5sJ6hC8+Ync0FiYgAuVYy3hW9qTpAUOt6zL5d+tBTesxLwefnUrJiBvdM83/UUeItkxcTeAcTH4pvnW2PssHSHm7BUOq1I52k31eMQ1FMFIQxNyIjEVyGBIwIQoMKTDhGGLQ+tJvIADzNslcL3ZYakxJlWwCrEygmILTdKwmGcss2L/h9gKcsiCmeB0WcgLC0vbOHu46SyXLVBxJlMFnMjMYAQeatRNAHKQQjm+KccCynKasC2nOsRSR/Ccue/LmkWw2BVhgSeX/KOaYAjrcUDNgWaTn0z71yNb6a/np4NFj0IPbD7V1ZkPHpBi0o+tFwnNa1yAZ0p1muoKhKReA9QMcVx34h5tgvBdNkhKkSAeWJDhHJDAkk0RiOuJOpZn0MJIIyIf3MidEv/tTPZ4efPMShwL4AIEFapnsG++WogWCTZ1qGFZ5Yn01ep5Mrtv8A58Uvz0cVlxa0edJBFzw77aItNsA39dntQPDvI299a9/7zQBIyVbRnCw9pG72vuvOhGb4A+KSsLAsQQeME8NRTkn3v+/FU7IjpnwtqWqv/fdoQkFaUmfRnv2qOczC0/zv1QOBBx/2gmrMobevfNTnxScLj9e98fU57ElXnmAex4/I/bwSG3YPdv8AE0JxgKSmlJez4SAZh5i37ohDLD4/f2vkH78LD1D8/wArq2DSxvsW428evz4NM808Pfx5U3q3l5dfj80zkSv5Ef37L18g/Bf58ADgaeWNJp+Pz6VnliRgO+2nXLTv5d/NAwSULADoRT4Vn1ez7LmsEv1PevpkrKEzOp14OHgKCt7byO49Pshscnn3aod4Cdf9/vhih0x0J/76T4JngKV+D271u8UKCA+1GOTTnr3jPNAHke4phjajSO90JCnQM/53E0UsjH99+2KBAgPtpNml07s09/8AIqFZ+f7UMCRMDN5xnjkr32rP+fY//8QAKBEBAAEDAgUEAwEBAAAAAAAAAREAITFBUWFxgZHwEKGxwTDR8SDh/9oACAECAQE/EPxiWSkX04lSC3EJEh0zeJnTa/5FPh0M40xhvTDoakSdvOPrYadalWgqWIgiPv8AIFWGw7ZHeJ7bU4eV1f8AC6UN4RBqH84zY/AghMhdPPGm/E0R63Vi228+vwKXE60F/NcNur/gIwJkbMeLZ60hx5MTrGnUs59RUKWksdvR6GXAfdMGc3BbYv3tl40JIHK+vO9Bsb6oxAuduNG0GbNX3468P8S5AfbE/W5V8eW/Fj9f4igTgPunEq0nrs7fNIxMuGM/c0sOSXZkg254/dMcd82f15p6vFHbyJ0pAcSvf/AvkJubUCM3W8nQ/nqjwF6hrctsEx7fEb0eiRb3mSLZ05b0rCRgFydeNumjTIBATBGZnglLKPVjptrV9+/PQoOaUhzpi0eMYrOev4dD1nMb2ChSUjtCO8VKBx0zhcVunk/2h0c1PPijIWQea6Vo+xztiP8AlBKQPEn7KPNKQYtMWw6fPvT8eMRvudP5SQkyJc5YccbZikquFnc45zpzhvWGTqbP+SAjC4mSkmxIffhRQBvZyOlm6ae/IHbmmPfSpILEXcT87UDJI2tU4kTDHHHbVpAokymOu1YpItBZvuvFe6/dXWzcZq18GEt9p7POhMWOBED8+cawIaTdPONRrlS3oCZXztQJLgGxz/72ihQBnJBbOkHIzyqAnlECERDoQmNFmgwHpghmNHmc4vJhgdsvOtTJU4l/O9DAXYRxHG0usRQnidSloyiL/JvYvTpCXg58dKUwnfm0zRIMSXJ28mgjSNRHDhFvO9NAsU7UO2pxfiBrQT3TRt5rebFPoQSvn9oCFmGznpcjkMUtV1B4QHwaFKjLQpuOQ5ONdKspq9df319MAvSpUq+DWynB93pHEhOyAnXSO+lK1LRvE7yR90aFIstuZDrTy/0HTy3K1FMbQUbnojJLQeDGDP1obtS2DZB0eDRCYQTYsPF3qOIAvBa/z70Pw4NalIFrq8v5wJpkxY3cC/8A3hmn5oFUW6aPent6PF0POG9R+5BoeZachh41gr4+KUyI6taCPXlW+aXZmO1vr1Ku4k73Pie1ETeikAAEIvyJ5Ukme8CncJfw3bOxAisJrN13ZePAHMSS7/y1LnnApg4fviUhzGGxIz99rZ/0SZMpqPiJL3Z9fZfT0MiSEvYt7xTuTL/PaiAYIHjJ/wB9IrsoNJv9VMQ1b3Ludrfi9l9PS7MMnct7xTvaNWwtMO+hMYtidoPQXmjG/jy9IllkbGvnKimmKL6sx8fhOtyHKVt2PRCSEoQ1QcH5+OtQTDJXK/rH6IoBLo24PntSskNPy/Mr80ZayAci/wCz8Jj7rH7j37/5vwEOIvMa8dvQ0crip32EvFeT1/C+kYeTntnmUK53lpDpa1viG0+gTYpvAcUedY9BoSulCzAY4vOuC16cPK/icvD2nbjscFxBSUITyThWNl/Kn8XzBpVko1rBGae33L+vI55pwsrq/jGLlAofAmHX77zS2wOn6KnFJQKRaM432aH0uBj/AL5H4P/EACkQAQEBAAIBAwQCAgIDAAAAAAERIQAxQRAgUTBAYYFxkcHRUPChsfH/2gAIAQEAAT8Q+3f5C8LhQTcYlWiaHgpw3cFkNehMzuE9nmgVJFCh/wACYgi2jerh2mkHhUzGzV6kBmiDafbgB5QB1WXcJeqmUVcGfrJpAFlMl9/ZLL6MwqIwFZg8fN1DNhLFA6gIDmRkJIK4LMbUUH3C8VrawrZChe5ZOtCUEZEReHdMYnCQhSpBKAvElWkV7UX4sPDyAURIijxX/vgdc9un57yaFcDBErQCkVfty2ELaSWh2mGjg82ReD6WQKoLqBzBYtaOpWjVRBDPoMWzIP2snFoaIHO+PdFMa6JQw7tD7aoWvuVXjHmB4K8pGRJkPJDoRcH6MfCgz5WoLYokIkcxeQtVVVFUFCiq/Y/2S8x2msL3Xx7n2GjwRFWh0wx0+6OMlYQpoACq4Bwbb50gsGmtqJjnAmIY6kGSxiojZ09xIhoCFGFDDYNCAO03yDZyAR0k+wGDQYwMRIOwi3U4neGWFBJOBRGJh6kxCQIoAFBS+X39dIpITZGHrp5XImlwIKj1nZJx/AR8IBFpVWgvJBHPbaCjLR4oAhyjE4jG9TrqISs6Flr7D3pcEKaAAquAcTn0/wDTqkpUhJeaQbqcBhhIBIvcemgo2oFwxCwQuFsp+FhYUag2GA+jFeLalEZ5ltUQI8PiSEpxoKlRQZ4PoHSgGKAHK7QVUwI2aqGjFRsS7Cn0jMuCwAoO6hYQJnBL7kAbJsCHouhgU8zIMggDZws4jjZCcKU8pNpEeABy1wIi1EFEcR5/vcfDfuPzTscj1eygR1PT0AOsL3T0QxOWNwgxf5AWCkpxBdA0FSYfj/G+AJmMWFlh73zJCEFEHRPBRUKgVWv6imbcKqJE4SAiiVC5UAYYSNn0ehSKEBWUxMKACgCBkbStd39N1ToOe0d1VVaoqhVKqqvI5C8MUwmxCLphQ97/AK6NFtNrpOJz2ruqIjoIoBQiCJzEhjbNONnYOlCeWQk476ggnKQAhYDwqD8rR2PIkAPaoQPIGfEVgB5FNuw4cTakDBUGHF0AGFHGgqVKDPB7wMvMCIpiCiOI8Fh7VgkEtCnRhyNeebAS6kAcQsAnv64+Qo7MTMSB55NO2zz3CgUi1QcWzVgMQtKlJEkLKTr3YTGJmdVOF7ZZHwo0YEY4ieOAAJTeW/4aEtznU5/LAfBcwEScEkbXnEunxVBaOZI03ujk/OSI9OCFJOcRwZSUHBH3BSzfoBAlCEWKgwfoiB1POxEOcj9B+Ba2+wISDERIeGe0bs5WkQR6wJb/AIRABIxiws6Pf3VZlOLpKUhaYCUQ1A4DS6IEFU4R23eoaoWQkqYAtyRoQdr2AwPIii5UCHxKTIMgdKWw90XQDSICEY6XhSaYx+Ts0aJ5xzNIQwQCAiSEUr8azQixIBRBrxyhZXUT1aEZoEKi5VILZc4ViA6JSy9gOlIxCVPQvImCCBOcWloXcZZQ4a5CQvSLLCijRQSV+l0qOKhR1gBXADx7xWMl8IlAOF0oVE5q3jADdJqRSuB1MCP+7AwsKueDntXdURGIIoBQiCJzw8Y2LoWRwXBJC5KhWRL+IFRQMSkErqYCWsheCLD8g3gyrQQREhKCPj9iDA4AosHv44kRqbQC3ZRYxOz06KxYJaq+AztVX6tJUAUhAz2CnqCeBWb6gqt+n2zZv/oBygrgADl19B1wU4o49EB4Qjhyz4PLo8oJYByw4QpkxjwROFRoQGrDrFgYUFagAriWNDEeBxEqNgyFIlEed+8zYUaogxxE8ey6NMmEeGUMhNHG3NmOCqHEMdCwnlrC4p1MYsQ1aP0bBWJAQ4AAa6obz/eQRFc2s7O/D9mmrKyk3Rsgqi0QEUQbyNiRHAdxY0lESUgFZTcm/kWDAyXLahS+BqLC1gHm/I0CqMEcRXoEJSqqRtdoiQtqkih6VRcGvS4AA1EERMR+ky4BAVAIqRYjPJxadJwCIlQgKzXgk2GEXQZVAAGoCqiDCdwqFSBSdVfPKx1lkMTYI441aQxI/wBCEuCpVhAB9Z8Nl2eUkMENlTIpvQRgQPFKY0GIFwkLUHxxCSEnCYBjmcoRrFayCtzSuigTHeExGHDX0rPHeK5jNJV4mdDwUEURRV4PzSBZqV5CkQSODCdXOuWZsINGAlMQfHN/CWPCIGjDtRy94oMrajGIKohQMfJv5I/0T5NUwVlgGHikoqJIYDk0ygNEErNEQkSoOwJDzc39KKOhz4Eo+W9iMsBgS+bOC+qIThHQOTsEi5JHBMUtJQfTp18bQpLj9/sdOGAENAar+WNA5FOQFSmBCUEQFKKfPMBxDXWi00O1ScX7YWbNss4lYsXB+KtKBjLdcFEKHHe6Jz5bXE7ACcZgZ2gFRZCSYjgsjKmMrWEJwkUY8NhAXJZWukAZlwUcxIAiIC/BR3P0ReRaRA/QSIKB9sI60Q0ebYRQI4a545oenYSglkozJHYVgD1cHVo+uVRSNphLxaO/ABQZuC0ONQoiBCDxMCnOBhW4mRg0vCJKs3Axdm1KrGdCowVSSS0IhgcoG0CZ0KoIhruHM2xrnzCRMBtNBy8u4jyUuRcSh5x4M01CCuju2Eeie8bGrriBEqQslXxw5QITQCO20oUCHCkGUXKpUwB0clIjc8nemBJ9SmbKWYBVzLZCZ6MUzsFXsORY7fA9kFFCEONi1k7kUoQnlFuQRYD6DcpNJEM+CnY3Sc6k7huuUlt7xVOPvkUkCzaEhU+XMBRBUTJYBWjbQLikfEg9CKgaQFLnIjpGZiFjNC4hQgoa6nJ4BFgj2beud6Jk9TiO0dG8wbDQmGiTFpjChxJgrBGDRg+IKgCc4pQAMHCQAcCoPqmzTJAHlkHaDD7FXC0kXmwqDYMHhvg34KSYt2NHQhjLFIoqqLCMQI49P/KN/DciarcFiFTMYiRiOVIC4BIhukLZO3HxClBkVKw+3zhqj04D4dNAgIBSjETWfY2idIbpNgUzhyGEr4K9igFGtIhwGVwQdCO9IPynont7AgjnY5WUSjxA9SRVSkMeAhRPQKFKQAF2CHQVQp431pH7sEAcoR2SBSrCCGIoaBYKhwsJwiVoTcjh2vMILzSJQ9IRAgwrKYlmh+1UEqMHn9oDMdprC918eybcYPCBEiYzpnOQDE8T13ngsBce9SkKlzgAALjl/TapTET0ixkSp6SQgRAEQUA51Y1MQXRsJslZcvx1lX/AvyYCj0XG4bW3CEwjqjzJxo3dBemq8DL4anRPgBL+iR5Tmf8AGqRodAN1UXnLeMoauCCaIIKvPnRSvvISrB4UY8/ZEY5M4wz5FuT9YFGn5X2uId3ECv8AHnwkgBKQCkcAcRETlo0IFgIUogs+iRTr+dJA9phA/h2YbDRHR9JUojKQQQQQ46/o7fAAZ2LeAC7gzknpA76DE6EtbBye9LggDQAAAwDgKECHZJGdqFMIOUfdLd1caAILDC+3Yji0PNe5yBOd0TD0SUAIPjJSeoIwTqBBFAAXuL55j27Rq/AWzdEs5ooXTXoXYSFbmTwiSlJgQCZqawVwDOHjF+aBQV3EGVBw0d1PFEE+WDU5Pd/AmEGCosNVfPpfKbO90JNRNhOSmtFlxQIJWRlGAusHFqQ7MLYkVg+8byIg0qoEBjxS+PVUzVRkABOePH5AhBV5yGBQudc57ivWukoHGvsRaG78ss8hO4vRy3RLRuhu9EaBCjvmWxVHTOWHS9BMiIQTT0Yt1ViFEkCchlNDYXdiHATJkmBE4VGhAasOrogkekYYs3rDkx3JDoJataIDJ79ZGwo1RBjiJ45P9iIHpChrGVQEADiTbSm5UBNgucQkl6MElskyItILlzWGTHs21inLU9mkhsC9BKPUFfYQdIQp1sBCTviwVqvyw8RADAAA9KA+DZ/8IIfCOTMhHdVVWqKoKFVVeHmSoIU0ABVcA4SxDEtsBQKsZRz3ZeBVZgY3iAsxDL/pCloCMKYqjlwCFvlTFookc4nG8a44QoA7G0DB4AhmqoAS4ylq78eZyOWXyIVaS4WWMOEKq1MDtwD6HVvMmEGqIKYg+OLHWMGJICDqFKKaxftlIGg/hhdcyT3gP2c1+BLdk5UELFxWhsYwRXia0T+aC5AqZRkYxQvisVwk2EUOHcaEMSEo/mSBeMWINUgtSUF6qfPFUk/sAF+4/mdjmj+vAjKjKA2mse5f45+IZjC9V8+jTX/ygixtoMCJw4PaHYi3xAJaiMDpGToCgOXlImonA+KjLfqiiS1zE41OuBICGwGFcA3gs39ZI0A0jxJm9LhILqgDlKMTm/Rs9PAcBmAPktMowKQyRDGh7U8CyD5lkFQFtd9KuErSxK0aqGsN407edSIYIMCaugpUREQeruPKa5DFaIWlIOD/AHwNlv4q/Feh7P8AJg44v4ympslaKZFYMNDhVaEEjaJnRCHDQ81GBFaSPGP4CCobHQAGlv0BVWviQLgQxQrA+in/AIhhSKmiiBiDAOffOl0Yx3ELCE4+P+iB9d7p2/HOLxMLI6lXlFHAAILL16jcUR284LxQ9LWfCgKqoXnZ7JZfRmFRGArMHg4XJYWMWzVKCJPoTLi4YCApUNifDwySMsHzZ4O94dP1cLmxx7DG+LNxS080gI7/AD2JDTYUfaKJAhtAH3mCzGepe26+CoSMYCnAVnHvyt8KZBh1sBYJzNyHu0qq1RVBQqqrwvwdSmPxlgUQxI0vlOEUfwQDICmpCZy8DUgIgrkPJmWxREjIQCgp4NlDdwPcMEhZQlIGhv2p4DTtqU4FCXGwYuzIKUiEYVCsUIey7YbQHhtPiqAKeAk3THLqThnlbJhAMKIn6nVtB5BpF7ND0sUev6GyLCp0Yh2lVfoq3G1vIGFtF0EAgds8PozCojAVmHG2DaL3IUCgLVBwGBJe07kAAA4WVTEAUlpgAdEHZ4Bg8maCNrvxEh6wmwtQzDNi30JQ8VMunYftGdYCTBpogjWZIHuNaOOJNbHhiCCiBqAo+4KzlbjkpWFToVGhKGcEwogTKQEoGdARPNR3JdHuhWwPYOArAAsAQZMEFoI+gwFpCiCDRREVgoil0q/yw8RADAAByjjJCQlNTDF2BgsH6fVfQmEGqAlMQfHDjJWEAaAAAGAfTF/ziwfhrO21vgdl9JCU40FSooM8Ho9pGTBkTsh0glE9hnATYqZVOFQ7HojShYkHRBZ0KUOCl0P8Oj0OghmR+4VrqIosSJBEEVpQQoWg9DQKQALFL5fT8r8u+ebCFhlVI8XoRhr2rDQlhHok7xyeox4IIdSYlxC2kbQ6TDQxOdksvozCojAVmHO5L9LjFLW8SFujnWMLwBgDQlDuVPtlf4MgCLqezgWGsO5ssHiWpKBVDnoHGSsIU0ABVcA4GpWVeTBYtQaYT6LAQ5w+41kK1EAP1wtGwMaQzG8V5jJCYVQVHSiw04dK1RdBgsNKC0B9GAwDxUYAxmGLpT7HXSD3csihQF7rCPiSEJwqoFQVZ5fZh4CVR0YiYNDxyyQQ3wlIYAigkbh120DKCtQAVxLGhlALd1NqSVBGCDlas1dqFaZMgUppxguzzgTVACTuD59DI72sApXdiqg0IuS6dqh0QAIxK3mQS+ryhUYiGkI2cF7CO7lkaFSnVK37A2RrGgokgpUawNZ++MB2GH4szFZX2UWbSrqUoxUQ0zjIImXKjugwwDVd5QFbVFQfDOolfRNN0tF0gWELBZz+LeD8Nb22k8Hp7JGsYnAKsZICy9HC/FdURSwlfKR1aH2EIc8MB+1LBoUhP65QQhBI9UdMPoj05GISp6FtEwQR53CoqratBCQz2qiWdE6z5JYYQXgMcdXBtMKGCh2/Y9k8vo3AKhREuJxbFKZoYtMlDWUP0FatUtEFJEtUoioGZTOpEQsWAgWIfWsMpspaJU7FAmgaTm5Wf3ZAP3A/ZD5/B/JFW2B0LeFG59EE/XbBqwQ+pIoYwiKwmBhQAUAO9p8kCyCKcsdxdl+UKAEhgKwAAAH2VMITwktDpMMHF5FAgaGIUxgFQgfacZKwhTQAFVwDjJawwEWFVhQvY8PcKjFSAkQol9iR+lAZLERRACheAO0GRNoLB3ALoE4538MEMUrsBKBd+1Ea1pci2aIiq01ZdSdNoaZhb8a617GubXeIU8wYCVEmZZUCkFECBIJUhgcZKwgDQAAAwD18gCFkFHYesLrmq6C6EEtrLU6JYlZP2lHswFSgDAH21oRYhEMLBEhESs+ZqIDKkQEFLIJ+4Ep4LxCkNDeoGgNjD6QW39jSRKkBnKH6lVEBS+iVA+5AZeYERTEFEcR4RGZggEoiIOgUjglrA5SKSiGUBUmAmUmKrB+YZAPKccLsrgqSmMBTgKzhRAnDCjmoQ0HUN5jnJWlBZRzMRkeRjtpCCUsEAegO4YnTc39pKehfuw0p6h6zoQQCAlOCPpO+LDK8GK4kttlKsT7D/9k=";

    let viewer = null;

    function setStatus(state, message) {
      const img  = document.getElementById(\'status-img\');
      const txt  = document.getElementById(\'status-text\');
      const err  = document.getElementById(\'error-text\');
      img.classList.remove(\'spinning\');

      if (state === \'loading\') {
        img.src = IMG_LOADING; img.style.display = \'block\';
        img.classList.add(\'spinning\');
        txt.textContent = message || \'Looking up…\';
        txt.style.display = \'block\';
        err.style.display = \'none\';
      } else if (state === \'success\') {
        img.src = IMG_SUCCESS; img.style.display = \'block\';
        txt.textContent = message || \'\';
        txt.style.display = message ? \'block\' : \'none\';
        err.style.display = \'none\';
      } else if (state === \'error\') {
        img.src = IMG_ERROR; img.style.display = \'block\';
        txt.style.display = \'none\';
        err.textContent = message;
        err.style.display = \'block\';
      } else {
        img.style.display = \'none\';
        txt.textContent = message || \'\';
        txt.style.display = \'block\';
        err.style.display = \'none\';
      }
    }

    function lookup() {
      const q = document.getElementById(\'query\').value.trim();
      if (!q) return;

      document.getElementById(\'table\').innerHTML = \'\';
      document.getElementById(\'compound-header\').style.display = \'none\';
      document.getElementById(\'viewer\').style.display = \'none\';
      setStatus(\'loading\', \'Consulting the magician…\');

      fetch(\'/lookup?q=\' + encodeURIComponent(q))
        .then(r => r.json())
        .then(data => {
          if (data.error) {
            setStatus(\'error\', data.error);
            return;
          }
          render3D(data.sdf);
          renderTable(data.props);
          setStatus(\'success\', \'\');
        })
        .catch(e => setStatus(\'error\', \'Request failed: \' + e));
    }

    function render3D(sdf) {
      const el = document.getElementById(\'viewer\');
      el.style.display = \'block\';
      el.innerHTML = \'\';
      viewer = $3Dmol.createViewer(el, { backgroundColor: \'white\' });
      viewer.addModel(sdf, \'sdf\');
      viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.22 } });
      viewer.zoomTo();
      viewer.render();
    }

    function renderTable(props) {
      const skip = new Set([\'CID\', \'IUPAC_Name\', \'Common_Name\', \'SMILES\']);
      const hdr = document.getElementById(\'compound-header\');
      hdr.style.display = \'block\';
      hdr.innerHTML =
        \'<div class="cid">CID \' + (props.CID || \'N/A\') + \'</div>\' +
        \'<div class="cname">\' + (props.Common_Name || props.IUPAC_Name || \'Unknown\') + \'</div>\' +
        \'<div class="iupac">\' + (props.IUPAC_Name || \'\') + \'</div>\';

      let html = \'<table><tr><th>Descriptor</th><th>Value</th></tr>\';
      html += \'<tr><td>SMILES</td><td style="word-break:break-all;font-size:11px;color:#555">\' + (props.SMILES || \'\') + \'</td></tr>\';
      for (const [k, v] of Object.entries(props)) {
        if (skip.has(k)) continue;
        const val = typeof v === \'number\' && !Number.isInteger(v) ? v.toFixed(4) : v;
        html += \'<tr><td>\' + k + \'</td><td>\' + val + \'</td></tr>\';
      }
      html += \'</table>\';
      document.getElementById(\'table\').innerHTML = html;
    }
  </script>
</body>
</html>
"""


def pubchem_interactive_cli():
    """
    pubchem_interactive — start a local web viewer.

    Launches a small Flask server at http://localhost:5050.
    Open that URL in your browser (Windows browser for WSL users).
    Type a CID, compound name, or SMILES, press Enter or click Look up.
    Press Ctrl+C in the terminal to stop the server.
    """
    try:
        from flask import Flask, request, jsonify
    except ImportError:
        print("Error: Flask is not installed. Run: pip install flask")
        sys.exit(1)

    app = Flask(__name__)
    app.logger.disabled = True
    import logging
    log = logging.getLogger("werkzeug")
    log.setLevel(logging.ERROR)

    @app.route("/")
    def index():
        from flask import Response
        return Response(_VIEWER_HTML, mimetype="text/html")

    @app.route("/lookup")
    def lookup():
        q = request.args.get("q", "").strip()
        if not q:
            return jsonify({"error": "empty query"})

        try:
            resolved = resolve_compound_input(q)
        except Exception as e:
            return jsonify({"error": str(e)})

        if resolved["mol"] is None:
            return jsonify({"error": f"Could not resolve a valid molecule for '{q}'"})

        cid = resolved["cid"]

        # Fetch 3D SDF — only possible when we have a CID
        sdf = None
        if cid:
            url_3d = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
                "/SDF?record_type=3d"
            )
            resp = requests.get(url_3d, timeout=10)
            if resp.status_code == 200:
                sdf = resp.text
        if not sdf:
            return jsonify({"error": (
                f"No 3D structure available for this compound "
                f"(CID={cid or 'N/A'}). "
                "PubChem 3D conformers are only available for compounds in their database."
            )})

        props_raw = get_rdkit_dict(q)
        if props_raw is None:
            return jsonify({"error": f"Could not compute descriptors for '{q}'"})

        props_clean = {
            k: round(v, 4) if isinstance(v, float) else v
            for k, v in props_raw.items()
        }
        return jsonify({"sdf": sdf, "props": props_clean})

    port = 5050
    print(f"\n  PubChem viewer running at: http://localhost:{port}")
    print("  WSL users: open that URL in your Windows browser.")
    print("  Press Ctrl+C to stop.\n")
    app.run(host="0.0.0.0", port=port, debug=False)


# ---------------------------------------------------------------------------
# CLI: rdkit_check_descrpt
# ---------------------------------------------------------------------------

def rdkit_check_descrpt_cli():
    """rdkit_check_descrpt <cid|name|smiles> <func|key>

    Evaluate a single RDKit descriptor on a compound and print a labelled result.

    The compound can be given as a PubChem CID, a common/IUPAC name, or a
    SMILES string.  The descriptor argument can be either a DEFAULT_DESCRIPTORS
    key (e.g. TPSA, HBA, NumRings) or any dotted RDKit callable.

    Examples
    --------
      rdkit_check_descrpt 3033 TPSA
      rdkit_check_descrpt aspirin NumRings
      rdkit_check_descrpt "CC(=O)O" rdMolDescriptors.CalcTPSA
      rdkit_check_descrpt 2244 Fragments.fr_COO
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="rdkit_check_descrpt",
        description="Evaluate one RDKit descriptor on a single compound",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "DEFAULT_DESCRIPTORS keys:\n  "
            + ", ".join(DEFAULT_DESCRIPTORS) + "\n\n"
            "Supported RDKit namespaces:\n"
            "  Chem, rdMolDescriptors, Fragments, Descriptors, GraphDescriptors\n\n"
            "Examples:\n"
            "  rdkit_check_descrpt 3033 TPSA\n"
            "  rdkit_check_descrpt aspirin rdMolDescriptors.CalcNumAromaticRings\n"
            '  rdkit_check_descrpt "CC(=O)O" HBA\n'
            "  rdkit_check_descrpt cannabidiol Fragments.fr_Ar_OH"
        ),
    )
    parser.add_argument(
        "compound",
        help="PubChem CID (e.g. 3033), compound name (e.g. aspirin), or SMILES string",
    )
    parser.add_argument(
        "function",
        help=(
            "Descriptor to compute: a DEFAULT_DESCRIPTORS key (e.g. TPSA) "
            "or a dotted RDKit callable (e.g. rdMolDescriptors.CalcTPSA)"
        ),
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

    func_label = args.function
    if args.function in DEFAULT_DESCRIPTORS:
        func = DEFAULT_DESCRIPTORS[args.function]
    else:
        try:
            func = _resolve_rdkit_func(args.function)
        except (ValueError, TypeError) as e:
            print(f"Error: {e}")
            print(
                f"Hint: valid DEFAULT_DESCRIPTORS keys are: {', '.join(DEFAULT_DESCRIPTORS)}\n"
                "Or use a dotted RDKit path such as rdMolDescriptors.CalcTPSA."
            )
            sys.exit(1)

    try:
        result = func(mol)
    except TypeError as e:
        print(f"Error calling {func_label}(mol): {e}")
        print("Hint: this function may require extra arguments. Check the RDKit docs.")
        sys.exit(1)

    cid_display = str(resolved["cid"]) if resolved["cid"] else "N/A (not in PubChem)"
    val_display = f"{result:.6f}" if isinstance(result, float) else str(result)

    print(f"\nCID         : {cid_display}")
    print(f"IUPAC_Name  : {resolved['iupac_name'] or 'N/A'}")
    print(f"Common_Name : {resolved['common_name'] or 'N/A'}")
    print(f"SMILES      : {resolved['smiles'] or 'N/A'}")
    print(f"Function    : {func_label}")
    print(f"Value       : {val_display}\n")


# ---------------------------------------------------------------------------
# CLI: rdkit_default_descrpts
# ---------------------------------------------------------------------------

def rdkit_default_descrpts_cli():
    """rdkit_default_descrpts <cid|name|smiles>

    Print CID, IUPAC name, common name, SMILES, and all DEFAULT_DESCRIPTORS
    for a compound.  Accepts a PubChem CID, common/IUPAC name, or SMILES string.
    If the SMILES is valid but not in PubChem, identifiers are shown as N/A
    but descriptors are still computed.

    Examples
    --------
      rdkit_default_descrpts 2244
      rdkit_default_descrpts aspirin
      rdkit_default_descrpts "CC(=O)Oc1ccccc1C(=O)O"
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="rdkit_default_descrpts",
        description=(
            "Print CID, IUPAC name, common name, SMILES, and all default "
            "RDKit descriptors for a compound"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  rdkit_default_descrpts 2244\n"
            "  rdkit_default_descrpts aspirin\n"
            '  rdkit_default_descrpts "CC(=O)Oc1ccccc1C(=O)O"'
        ),
    )
    parser.add_argument(
        "compound",
        help="PubChem CID, compound name, or SMILES string",
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
        "CID":         cid_display,
        "IUPAC_Name":  resolved["iupac_name"] or "N/A",
        "Common_Name": resolved["common_name"] or "N/A",
        "SMILES":      resolved["smiles"] or "N/A",
    }
    for desc_name, func in DEFAULT_DESCRIPTORS.items():
        try:
            val = func(mol)
            row[desc_name] = f"{val:.4f}" if isinstance(val, float) else val
        except Exception as e:
            row[desc_name] = f"ERROR: {e}"

    display_table(row)


# ---------------------------------------------------------------------------
# CLI: rdkit_batch_fetcher
# ---------------------------------------------------------------------------

def rdkit_batch_fetcher_cli():
    """rdkit_batch_fetcher — batch-compute RDKit descriptors from a CSV file.

    Input CSV must contain at least one of: CID, Name, SMILES columns.
    All three may be present; see resolution rules below.

    Resolution rules
    ----------------
    Priority (highest first): SMILES > CID > Name
    Each field is validated independently.  Invalid fields are flagged and
    ignored.  If two or more valid fields point to different PubChem compounds,
    a MISMATCH flag is printed, but processing continues using the
    highest-priority valid input.

    If the resolved compound differs from an input field value (e.g. input CID
    was used but SMILES was also given and resolves to a different CID), the
    resolved values are written as CID_resolved, Name_IUPAC_resolved,
    Name_common_resolved, SMILES_resolved alongside the original columns.

    Output CSV columns
    ------------------
    All original columns (passed through unchanged) + CID_resolved +
    Name_IUPAC_resolved + Name_common_resolved + SMILES_resolved +
    all DEFAULT_DESCRIPTORS + any extra DescName::RDKitFunction columns.

    Extra custom descriptors
    -----------------------
    Add a column with header  DescName::RDKitFunction  to the input CSV.
    Example column header:   Chi0v::Chem.rdMolDescriptors.CalcChi0v

    Examples
    --------
      rdkit_batch_fetcher molecules.csv results.csv
      rdkit_batch_fetcher molecules.csv results.csv --formats xyz sdf
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="rdkit_batch_fetcher",
        description="Batch-compute RDKit molecular descriptors from a CSV",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Input CSV must have at least one of: CID, Name, SMILES\n\n"
            "Extra descriptor columns use '::' in the header:\n"
            "  Chi0v::Chem.rdMolDescriptors.CalcChi0v\n"
            "  NumHeteroatoms::rdMolDescriptors.CalcNumHeteroatoms\n\n"
            "Examples:\n"
            "  rdkit_batch_fetcher molecules.csv results.csv\n"
            "  rdkit_batch_fetcher my_list.csv out.csv"
        ),
    )
    parser.add_argument("input_csv",  help="Input CSV file path")
    parser.add_argument("output_csv", help="Output CSV file path")
    args = parser.parse_args()

    if not os.path.isfile(args.input_csv):
        print(f"Error: input file not found: {args.input_csv}")
        sys.exit(1)

    with open(args.input_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        headers = list(reader.fieldnames or [])
        rows    = list(reader)

    id_cols = {"CID", "Name", "SMILES"}
    if not id_cols.intersection(headers):
        print("Error: input CSV must contain at least one of: CID, Name, SMILES")
        sys.exit(1)

    # Separate extra computed columns (DescName::func) from passthrough columns
    extra_descs:      dict[str, str] = {}
    passthrough_cols: list[str]      = []
    for h in headers:
        if "::" in h:
            desc_name, func_str = h.split("::", 1)
            extra_descs[desc_name.strip()] = func_str.strip()
        else:
            passthrough_cols.append(h)

    resolved_extra: dict[str, object] = {}
    for desc_name, func_str in extra_descs.items():
        try:
            resolved_extra[desc_name] = _resolve_rdkit_func(func_str)
        except (ValueError, TypeError) as e:
            print(f"Warning: skipping extra descriptor '{desc_name}': {e}")

    output_rows: list[dict] = []
    total = len(rows)

    for i, row in enumerate(rows, 1):
        res = _resolve_batch_row(row, i, total)
        if res is None:
            continue

        mol         = res["mol"]
        smiles      = res["smiles"]
        cid         = res["cid"]
        iupac_name  = res["iupac_name"]
        common_name = res["common_name"]

        # Build output row: start with all original passthrough columns
        out_row: dict = {}
        for col in passthrough_cols:
            out_row[col] = row.get(col, "")

        # Append resolved identifier columns
        out_row["CID_resolved"]          = cid if cid else ""
        out_row["Name_IUPAC_resolved"]   = iupac_name
        out_row["Name_common_resolved"]  = common_name
        out_row["SMILES_resolved"]       = smiles

        # Default descriptors
        for desc_name, func in DEFAULT_DESCRIPTORS.items():
            try:
                out_row[desc_name] = func(mol)
            except Exception as e:
                out_row[desc_name] = f"ERROR: {e}"

        # Extra user-defined descriptors
        for desc_name, func in resolved_extra.items():
            try:
                out_row[desc_name] = func(mol)
            except TypeError as e:
                out_row[desc_name] = f"ERROR (extra args needed): {e}"
            except Exception as e:
                out_row[desc_name] = f"ERROR: {e}"

        output_rows.append(out_row)
        print(f"    → done")

    if not output_rows:
        print("No rows processed — output file not written.")
        return

    # Build fieldnames from union of all row keys (preserves order, handles gaps)
    seen:       dict[str, None] = {}
    for out_row in output_rows:
        seen.update(dict.fromkeys(out_row.keys()))
    fieldnames = list(seen.keys())

    with open(args.output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"\nSaved {len(output_rows)}/{total} rows → {args.output_csv}")


# ---------------------------------------------------------------------------
# CLI: fetch_xyz_batch
# ---------------------------------------------------------------------------

def fetch_xyz_batch_cli():
    """fetch_xyz_batch — batch-generate 3D structure files from a CSV.

    Input CSV must contain at least one of: CID, Name, SMILES columns.
    The same resolution rules as rdkit_batch_fetcher apply (SMILES > CID > Name).

    File naming
    -----------
    Files are named using the first available of:
      1. Code column value (if --code-col is specified and non-empty)
      2. CID_CommonName  (e.g. 3033_Diclofenac)
      3. CID only        (if no common name)
      4. CommonName only (if no CID)
      5. missing_id<N>   (if no CID, no name, and SMILES not found in PubChem)

    Examples
    --------
      fetch_xyz_batch molecules.csv ./structures/
      fetch_xyz_batch molecules.csv ./structures/ --format xyz sdf
      fetch_xyz_batch molecules.csv ./structures/ --format xyz pdb --code-col Abbreviation
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="fetch_xyz_batch",
        description="Batch-generate 3D structure files from a CSV",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Input CSV must have at least one of: CID, Name, SMILES\n\n"
            "Examples:\n"
            "  fetch_xyz_batch molecules.csv ./structures/\n"
            "  fetch_xyz_batch molecules.csv ./structures/ --format xyz sdf\n"
            "  fetch_xyz_batch molecules.csv ./out/ --format sdf pdb --code-col Abbreviation"
        ),
    )
    parser.add_argument("input_csv",  help="CSV with CID, Name, and/or SMILES columns")
    parser.add_argument("output_dir", help="Directory to write structure files into")
    parser.add_argument(
        "--format", "-f",
        nargs="+",
        choices=["xyz", "sdf", "pdb", "mol"],
        default=["xyz"],
        metavar="FORMAT",
        help="Output format(s): xyz sdf pdb mol  (default: xyz)",
    )
    parser.add_argument(
        "--code-col",
        default=None,
        metavar="COLUMN",
        help=(
            "Name of a CSV column to use as the output filename stem "
            "(e.g. Abbreviation). Overrides auto-naming."
        ),
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input_csv):
        print(f"Error: input file not found: {args.input_csv}")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.input_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        headers = list(reader.fieldnames or [])
        rows    = list(reader)

    id_cols = {"CID", "Name", "SMILES"}
    if not id_cols.intersection(headers):
        print("Error: input CSV must contain at least one of: CID, Name, SMILES")
        sys.exit(1)

    # Import here to avoid circular dependency concerns
    from mof_toolkit.molecule_manager import get_3d_structure

    total           = len(rows)
    missing_counter = 0

    for i, row in enumerate(rows, 1):
        res = _resolve_batch_row(row, i, total)
        if res is None:
            continue

        # Determine file stem
        code_val = (row.get(args.code_col, "") or "").strip() if args.code_col else ""
        if code_val:
            stem = code_val
        else:
            cid         = res["cid"]
            common_name = res["common_name"] or res["iupac_name"]
            # Sanitise name for use in filenames
            safe_name = (
                common_name.replace(" ", "_").replace("/", "-").replace("\\", "-")[:40]
                if common_name else ""
            )
            if cid and safe_name:
                stem = f"{cid}_{safe_name}"
            elif cid:
                stem = str(cid)
            elif safe_name:
                stem = safe_name
            else:
                missing_counter += 1
                stem = f"missing_id{missing_counter:02d}"

        print(f"    Writing {', '.join(args.format)} → {stem}.*")
        get_3d_structure(
            query=res["smiles"] or (str(res["cid"]) if res["cid"] else ""),
            formats=args.format,
            output_stem=os.path.join(args.output_dir, stem),
            source="pubchem" if res["cid"] else "rdkit",
            cid=res["cid"],
        )
