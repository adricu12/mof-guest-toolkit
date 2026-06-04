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
    "Aliphatic_OH":     Fragments.fr_Al_OH,
    "Aromatic_OH":      Fragments.fr_Ar_OH,
    "COO":              Fragments.fr_COO,
    "Carbonyl_O":       Fragments.fr_C_O_noCOO,
    "Aromatic_N":       Fragments.fr_Ar_N,
    "Amine_NH2":        Fragments.fr_NH2,
    "Amine_NH":         Fragments.fr_NH1,
    "Amine_N":          Fragments.fr_NH0,
    "Ether":            Fragments.fr_ether,
    "Sulfonamide":      Fragments.fr_sulfonamd,
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
  <title>Molecule Explorer</title>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
  <style>
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: sans-serif; background: #f0f2f5; }

    header {
      background: #1e2a3a;
      padding: 12px 24px;
      display: flex; align-items: center; gap: 16px;
    }
    .header-mascot {
      height: 90px; width: auto;
      background: #f0f2f5;
      border-radius: 8px;
      padding: 4px;
    }
    .header-text { flex: 1; }
    .header-text h2 {
      font-size: 32px; font-weight: 700;
      color: #e8edf2; letter-spacing: 0.3px; margin: 0;
    }
    .header-links {
      font-size: 18px; color: #7a9ab5; margin-top: 4px;
    }
    .header-links a {
      color: #7ab5d8; text-decoration: none;
    }
    .header-links a:hover { color: #aad4f0; text-decoration: underline; }

    .main {
      padding: 20px 24px;
      display: flex; gap: 24px;
      align-items: flex-start; flex-wrap: wrap;
    }
    .left-panel { display: flex; flex-direction: column; gap: 12px; flex-shrink: 0; }

    .search-row { display: flex; gap: 0; }
    input {
      padding: 15px 15px; font-size: 24px;
      border: 1px solid #ced4da;
      border-right: none;
      border-radius: 6px 0 0 6px;
      width: 300px; background: #fff; color: #1e2a3a;
      outline: none;
    }
    input:focus { border-color: #4a6fa5; box-shadow: -2px 0 0 2px #4a6fa520; }
    button {
      padding: 15px 20px; font-size: 24px;
      background: #2c4a6e; color: #fff;
      border: none; border-radius: 0 6px 6px 0;
      cursor: pointer; transition: background 0.15s;
    }
    button:hover { background: #1e3550; }

    #status-area {
      min-height: 96px;
      display: flex; align-items: center; gap: 12px;
    }
    #status-img { height: 140px; width: auto; display: none; }
    #status-text { font-size: 24px; color: #4a6080; font-style: italic; }
    #error-text  { font-size: 24px; color: #c0392b; display: none; }

    #viewer {
      width: 700px; height: 450px;
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
      font-size: 20px; color: #7a9ab5;
      text-transform: uppercase; letter-spacing: 0.5px;
    }
    #compound-header .cname { font-size: 24px; font-weight: 600; margin: 2px 0; }
    #compound-header .iupac { font-size: 20px; color: #9ab4cc; }

    table {
      width: 100%; border-collapse: collapse;
      background: #fff; border-radius: 0 0 8px 8px;
      overflow: hidden; border: 1px solid #d0d7e0;
      border-top: none; font-size: 24px;
    }
    th {
      background: #2c4a6e; color: #c8d8ea;
      padding: 8px 10px; text-align: left;
      font-weight: 500; font-size: 18px;
      text-transform: uppercase; letter-spacing: 0.4px;
    }
    td { padding: 8px 14px; color: #1e2a3a; border-bottom: 1px solid #eef0f3; }
    td:first-child { color: #4a6080; font-weight: 500; width: 46%; }
    tr:last-child td { border-bottom: none; }
    tr:nth-child(even) td { background: #f7f9fc; }

    @keyframes wiggle {
      0%,100% { transform: rotate(-6deg); }
      50%      { transform: rotate(6deg); }
    }
    .wiggling { animation: wiggle 0.5s ease-in-out infinite; transform-origin: bottom center; }
  </style>
</head>
<body>
  <header>
    <img class="header-mascot" src="/static/figures/glasses.png" alt="mascot">
    <div class="header-text">
      <h2>Molecule Explorer</h2>
      <div class="header-links">
        by <a href="https://github.com/adricu12" target="_blank">adricu12</a>
        &nbsp;·&nbsp;
        <a href="https://github.com/adricu12/mof-guest-toolkit" target="_blank">mof-guest-toolkit</a>
      </div>
    </div>
  </header>

  <div class="main">
    <div class="left-panel">
      <div class="search-row">
        <input id="query" type="text"
               placeholder="CID, name, or SMILES"
               onkeydown="if(event.key==='Enter') lookup()">
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
    let viewer = null;

    function setStatus(state, message) {
      const img = document.getElementById('status-img');
      const txt = document.getElementById('status-text');
      const err = document.getElementById('error-text');
      img.classList.remove('wiggling');

      if (state === 'loading') {
        img.src = '/static/figures/success2.png';
        img.style.display = 'block';
        img.classList.add('wiggling');
        txt.textContent = 'Consulting the magician\u2026';
        txt.style.display = 'block'; err.style.display = 'none';
      } else if (state === 'success') {
        img.src = '/static/figures/success.png';
        img.style.display = 'block';
        txt.style.display = 'none'; err.style.display = 'none';
      } else if (state === 'error') {
        img.src = '/static/figures/error.png';
        img.style.display = 'block';
        txt.style.display = 'none';
        err.textContent = message; err.style.display = 'block';
      } else {
        img.style.display = 'none';
        txt.textContent = message || '';
        txt.style.display = 'block'; err.style.display = 'none';
      }
    }

    function lookup() {
      const q = document.getElementById('query').value.trim();
      if (!q) return;
      document.getElementById('table').innerHTML = '';
      document.getElementById('compound-header').style.display = 'none';
      document.getElementById('viewer').style.display = 'none';
      setStatus('loading');

      fetch('/lookup?q=' + encodeURIComponent(q))
        .then(r => r.json())
        .then(data => {
          if (data.error) { setStatus('error', data.error); return; }
          render3D(data.sdf);
          renderTable(data.props);
          setStatus('success');
        })
        .catch(e => setStatus('error', 'Request failed: ' + e));
    }

    function render3D(sdf) {
      const el = document.getElementById('viewer');
      el.style.display = 'block'; el.innerHTML = '';
      viewer = $3Dmol.createViewer(el, { backgroundColor: 'white' });
      viewer.addModel(sdf, 'sdf');
      viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.22 } });
      viewer.zoomTo(); viewer.render();
    }

    function renderTable(props) {
      const skip = new Set(['CID','IUPAC_Name','Common_Name','SMILES']);
      const hdr = document.getElementById('compound-header');
      hdr.style.display = 'block';
      hdr.innerHTML =
        '<div class="cid">CID ' + (props.CID || 'N/A') + '</div>' +
        '<div class="cname">' + (props.Common_Name || props.IUPAC_Name || 'Unknown') + '</div>' +
        '<div class="iupac">' + (props.IUPAC_Name || '') + '</div>';

      let html = '<table><tr><th>Descriptor</th><th>Value</th></tr>';
      html += '<tr><td>SMILES</td><td style="word-break:break-all;font-size:18px;color:#555">'
            + (props.SMILES || '') + '</td></tr>';
      for (const [k, v] of Object.entries(props)) {
        if (skip.has(k)) continue;
        const val = typeof v === 'number' && !Number.isInteger(v) ? v.toFixed(4) : v;
        html += '<tr><td>' + k + '</td><td>' + val + '</td></tr>';
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

    app = Flask(__name__,
            static_folder=os.path.join(os.path.dirname(__file__), "static"))
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
