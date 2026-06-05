"""
rdkit_descriptors.py
--------------------
Fetch chemical information and compute molecular descriptors
from PubChem CIDs, compound names, or SMILES strings.

CLI
---
  interactive_explorer
      Launch a local web viewer (3D structure + descriptor table in browser).

  rdkit_descriptor -i <cid|name|smiles>
      Compute and print the default descriptor set for one compound.

  rdkit_descriptor -i <cid|name|smiles> -d TPSA MolLogP NumHDonors
      Compute specific descriptors by RDKit name.

  rdkit_descriptor -i <cid|name|smiles> --full
      Compute all ~210 RDKit descriptors.

  rdkit_descriptor -i <cid|name|smiles> --full --filter TPSA MolLogP
      Compute all ~210, return only the filtered subset.

  rdkit_descriptor -i <cid|name|smiles> [-d ...] [--full] [-o output.csv]
      Save any of the above to a CSV file.

  rdkit_descriptor --batch molecules.csv -o results.csv
      Batch mode: process a CSV with CID, Name, and/or SMILES columns.
      Supports --full, --filter, and -d flags.

  fetch_xyz_batch -i molecules.csv -o ./structures/
      Batch-generate 3D structure files from a CSV file.

Python helpers (scripts and notebooks)
---------------------------------------
  resolve_compound_input(query)
      Resolves CID / name / SMILES → {cid, iupac_name, common_name, smiles, mol}.

  get_all_rdkit_descriptors(mol)
      All ~210 RDKit descriptors for an RDKit Mol. NaN-safe, returns flat dict.

  get_rdkit_dict(query, full=False, descriptors=None, filter_keys=None)
      Descriptor dict for one compound (CID, name, or SMILES).
      - default (no flags): DEFAULT_DESCRIPTORS set
      - full=True: all ~210 RDKit descriptors
      - descriptors=["TPSA","MolLogP"]: specific descriptors by RDKit name
      - filter_keys=["TPSA"]: subset when full=True

  display_table(descriptor_dict)
      Print a descriptor dict as an aligned terminal table.
"""

import csv
import math
import os
import re
import sys

import pubchempy as pcp
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Fragments, rdMolDescriptors

# Suppress RDKit's stderr SMILES parse warnings during name/SMILES probing
_rdlogger = RDLogger.logger()
_rdlogger.setLevel(RDLogger.CRITICAL)


# ---------------------------------------------------------------------------
# Descriptor registries
# ---------------------------------------------------------------------------

# Default set — identity + drug-likeness descriptors shown in the viewer
# and returned when no flags are passed to the CLI or get_rdkit_dict.
# Uses RDKit's CalcMolDescriptors naming where possible.
DEFAULT_DESCRIPTORS: dict = {
    "MolecularFormula": rdMolDescriptors.CalcMolFormula,
    "ExactMolWt":       rdMolDescriptors.CalcExactMolWt,
    "NumHeavyAtoms":    lambda mol: mol.GetNumHeavyAtoms(),
    "RingCount":        rdMolDescriptors.CalcNumRings,
    "NumAromaticRings": rdMolDescriptors.CalcNumAromaticRings,
    "NumHBA":           rdMolDescriptors.CalcNumHBA,
    "NumHBD":           rdMolDescriptors.CalcNumHBD,
    "NumRotatableBonds":rdMolDescriptors.CalcNumRotatableBonds,
    "TPSA":             rdMolDescriptors.CalcTPSA,
}

# Subset used in the interactive viewer (same as DEFAULT_DESCRIPTORS here —
# functional group counts are excluded; add them via -d or --full if needed)
VIEWER_DESCRIPTORS: dict = DEFAULT_DESCRIPTORS


# ---------------------------------------------------------------------------
# PubChem fetchers
# ---------------------------------------------------------------------------

def fetch_cid_from_name(name: str) -> int:
    """Resolve a compound name to a PubChem CID."""
    compounds = pcp.get_compounds(name, "name")
    if not compounds:
        raise ValueError(f"No compound found for name: '{name}'")
    return compounds[0].cid


def fetch_smiles_from_cid(cid: int) -> str | None:
    """Fetch canonical SMILES for a CID. Returns None on any failure."""
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
    Fetch IUPAC name, common name, SMILES, formula, and weight for a CID.

    Returns dict with keys:
      CID, IUPAC_Name, Common_Name, SMILES, MolecularFormula,
      MolecularWeight_PubChem
    """
    prop_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
        "/property/IUPACName,MolecularWeight,MolecularFormula,CanonicalSMILES/JSON"
    )
    response = requests.get(prop_url, timeout=10)
    response.raise_for_status()
    p = response.json()["PropertyTable"]["Properties"][0]

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
            for s in synonyms:
                if not s.isdigit() and "-" not in s[:5]:
                    common_name = s
                    break
            if not common_name and synonyms:
                common_name = synonyms[0]
    except Exception:
        pass

    return {
        "CID":                     cid,
        "IUPAC_Name":              p.get("IUPACName", ""),
        "Common_Name":             common_name,
        "SMILES":                  p.get("CanonicalSMILES", ""),
        "MolecularFormula":        p.get("MolecularFormula", ""),
        "MolecularWeight_PubChem": p.get("MolecularWeight", ""),
    }


def _looks_like_formula(q: str) -> bool:
    """Return True if q looks like a molecular formula (e.g. C6H6, H2O, C9H8O4).

    Requires at least one digit so bare element symbols like 'CO' or 'NO'
    fall through to the name-lookup path instead.
    """
    return (
        any(c.isdigit() for c in q)
        and bool(re.match(r'^([A-Z][a-z]?\d*)+$', q))
    )


def fetch_cid_from_formula(formula: str) -> int:
    """Return the first PubChem CID matching a molecular formula.

    Takes the lowest CID (most well-known compound) when multiple structures
    share the same formula. Handles PubChem's async ListKey response.
    """
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/formula"
        f"/{formula}/cids/JSON?MaxRecords=5"
    )
    resp = requests.get(url, timeout=15)
    if resp.status_code != 200:
        raise ValueError(f"No compound found for formula: '{formula}'")
    data = resp.json()

    # PubChem may return an async ListKey for large result sets
    if "Waiting" in data:
        import time
        list_key = data["Waiting"]["ListKey"]
        time.sleep(2)
        poll = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey"
            f"/{list_key}/cids/JSON",
            timeout=10,
        )
        poll.raise_for_status()
        data = poll.json()

    cids = data.get("IdentifierList", {}).get("CID", [])
    if not cids:
        raise ValueError(f"No compound found for formula: '{formula}'")
    return int(cids[0])


def _lookup_by_smiles(smiles: str) -> tuple[int | None, str | None, str | None]:
    """Query PubChem for a CID matching the given SMILES."""
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


# ---------------------------------------------------------------------------
# Compound resolver
# ---------------------------------------------------------------------------

def resolve_compound_input(query: str) -> dict:
    """
    Resolve a CID, compound name, or SMILES string into a unified dict.

    Returns
    -------
    dict with keys:
      cid         (int | None)
      iupac_name  (str)
      common_name (str)
      smiles      (str | None)
      mol         (Chem.Mol | None)

    If the input is a valid SMILES not found in PubChem, cid/names are
    empty but mol and smiles are populated for local descriptor computation.
    """
    q = query.strip()

    if q.isdigit():
        cid = int(q)
        smiles = fetch_smiles_from_cid(cid)
        try:
            meta        = fetch_pubchem_metadata(cid)
            iupac_name  = meta.get("IUPAC_Name", "")
            common_name = meta.get("Common_Name", "")
        except Exception:
            iupac_name = common_name = ""
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        return {"cid": cid, "iupac_name": iupac_name,
                "common_name": common_name, "smiles": smiles, "mol": mol}

    mol_test = Chem.MolFromSmiles(q)
    if mol_test is not None:
        canonical = Chem.MolToSmiles(mol_test)
        cid, iupac_name, common_name = _lookup_by_smiles(q)
        if cid is None:
            print(
                f"  Warning: SMILES '{q}' is valid but not found in PubChem "
                "— CID, IUPAC_Name, and Common_Name will be empty."
            )
        return {"cid": cid, "iupac_name": iupac_name or "",
                "common_name": common_name or "", "smiles": canonical,
                "mol": mol_test}

    if _looks_like_formula(q):
        cid = fetch_cid_from_formula(q)
        smiles = fetch_smiles_from_cid(cid)
        try:
            meta        = fetch_pubchem_metadata(cid)
            iupac_name  = meta.get("IUPAC_Name", "")
            common_name = meta.get("Common_Name", "")
        except Exception:
            iupac_name = common_name = ""
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        return {"cid": cid, "iupac_name": iupac_name,
                "common_name": common_name, "smiles": smiles, "mol": mol}

    cid = fetch_cid_from_name(q)
    smiles = fetch_smiles_from_cid(cid)
    try:
        meta        = fetch_pubchem_metadata(cid)
        iupac_name  = meta.get("IUPAC_Name", "")
        common_name = meta.get("Common_Name", "")
    except Exception:
        iupac_name  = ""
        common_name = q
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    return {"cid": cid, "iupac_name": iupac_name,
            "common_name": common_name, "smiles": smiles, "mol": mol}


# ---------------------------------------------------------------------------
# Core descriptor helpers
# ---------------------------------------------------------------------------

def get_all_rdkit_descriptors(mol: Chem.Mol) -> dict:
    """
    Compute all ~210 RDKit descriptors for an RDKit Mol object.

    Uses RDKit's own naming convention (e.g. 'TPSA', 'MolLogP', 'NumHDonors').
    NaN values are replaced with empty strings so the result is always CSV-safe.

    Parameters
    ----------
    mol : Chem.Mol

    Returns
    -------
    dict mapping RDKit descriptor name → value (float, int, or '').

    Examples
    --------
    >>> from rdkit import Chem
    >>> from mof_toolkit.rdkit_descriptors import get_all_rdkit_descriptors
    >>> mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    >>> d = get_all_rdkit_descriptors(mol)
    >>> print(d["TPSA"], d["MolLogP"])
    >>> print(sorted(d.keys()))   # see all ~210 available names
    """
    raw = Descriptors.CalcMolDescriptors(mol)
    return {
        k: ("" if isinstance(v, float) and math.isnan(v) else v)
        for k, v in raw.items()
    }


def _find_in_full(all_desc: dict, name: str) -> tuple[str | None, object]:
    """
    Look up a descriptor name in CalcMolDescriptors output.
    Tries exact match, then case-insensitive. Returns (canonical_name, value)
    or (None, None) if not found.
    """
    if name in all_desc:
        return name, all_desc[name]
    lower_map = {k.lower(): k for k in all_desc}
    if name.lower() in lower_map:
        canonical = lower_map[name.lower()]
        return canonical, all_desc[canonical]
    return None, None


def _build_id_block(resolved: dict) -> dict:
    """Build the CID/name/SMILES identifier prefix shared by all outputs."""
    return {
        "CID":         resolved["cid"] if resolved["cid"] else "",
        "IUPAC_Name":  resolved["iupac_name"],
        "Common_Name": resolved["common_name"],
        "SMILES":      resolved["smiles"] or "",
    }


# ---------------------------------------------------------------------------
# get_rdkit_dict — main Python helper
# ---------------------------------------------------------------------------

def get_rdkit_dict(
    query,
    full: bool = False,
    descriptors: list[str] | None = None,
    filter_keys: list[str] | None = None,
) -> dict | None:
    """
    Return a descriptor dict for a single compound.

    Accepts a PubChem CID (int or numeric string), a compound name, or a
    SMILES string. When a valid SMILES is not in PubChem, CID/names are
    empty but all descriptors are still computed from the SMILES.

    Parameters
    ----------
    query : int or str
        PubChem CID, compound name, or SMILES string.
    full : bool, optional
        If True, compute all ~210 RDKit descriptors using RDKit's own naming.
        Default: False (uses DEFAULT_DESCRIPTORS).
    descriptors : list of str, optional
        Specific RDKit descriptor names to compute (e.g. ['TPSA', 'MolLogP']).
        Uses RDKit's CalcMolDescriptors naming — case-insensitive.
        When provided, overrides DEFAULT_DESCRIPTORS (full=False is ignored).
    filter_keys : list of str, optional
        When full=True, return only descriptors matching these names.
        Case-insensitive. Ignored when full=False.

    Returns
    -------
    dict with keys CID, IUPAC_Name, Common_Name, SMILES, and descriptor
    values; or None if no valid molecule could be obtained.

    Examples
    --------
    >>> from mof_toolkit.rdkit_descriptors import get_rdkit_dict, display_table
    >>> # Default set
    >>> display_table(get_rdkit_dict(3033))
    >>> display_table(get_rdkit_dict("aspirin"))
    >>> display_table(get_rdkit_dict("CC(=O)Oc1ccccc1C(=O)O"))  # SMILES

    >>> # SMILES not in PubChem — CID/names empty, descriptors still computed
    >>> display_table(get_rdkit_dict("C1CC1"))

    >>> # Specific descriptors by RDKit name
    >>> display_table(get_rdkit_dict("aspirin", descriptors=["TPSA", "MolLogP"]))

    >>> # All ~210 descriptors
    >>> display_table(get_rdkit_dict("aspirin", full=True))

    >>> # All ~210, filtered subset
    >>> display_table(get_rdkit_dict("aspirin", full=True,
    ...                              filter_keys=["TPSA", "MolLogP", "MolWt"]))

    >>> # Loop over a list — single compounds or batch
    >>> queries = [3033, "ibuprofen", "CC(=O)Oc1ccccc1C(=O)O"]
    >>> results = []
    >>> for q in queries:
    ...     r = get_rdkit_dict(q)
    ...     if r is not None:
    ...         results.append(r)
    >>> # Save to CSV
    >>> import csv
    >>> fieldnames = list(results[0].keys())
    >>> with open("my_results.csv", "w", newline="") as f:
    ...     writer = csv.DictWriter(f, fieldnames=fieldnames)
    ...     writer.writeheader()
    ...     writer.writerows(results)
    """
    try:
        resolved = resolve_compound_input(str(query))
    except Exception as e:
        print(f"  Error resolving '{query}': {e}")
        return None

    mol = resolved["mol"]
    if mol is None:
        print(f"  Error: could not obtain a valid molecule for '{query}'.")
        return None

    id_block = _build_id_block(resolved)

    # Case 1: specific descriptor names requested
    if descriptors is not None:
        all_desc = get_all_rdkit_descriptors(mol)
        result = dict(id_block)
        for name in descriptors:
            canonical, value = _find_in_full(all_desc, name)
            if canonical is not None:
                result[canonical] = value
            else:
                print(f"  Warning: descriptor '{name}' not found — skipping.")
                print(f"  Tip: run get_all_rdkit_descriptors(mol).keys() "
                      f"to see all ~210 available names.")
        return result

    # Case 2: full ~210 descriptors
    if full:
        all_desc = get_all_rdkit_descriptors(mol)
        if filter_keys:
            lower_map = {k.lower(): k for k in all_desc}
            all_desc = {
                lower_map[f.lower()]: all_desc[lower_map[f.lower()]]
                for f in filter_keys
                if f.lower() in lower_map
            }
        return {**id_block, **all_desc}

    # Case 3: default set
    result = dict(id_block)
    for name, func in DEFAULT_DESCRIPTORS.items():
        try:
            result[name] = func(mol)
        except Exception as e:
            result[name] = f"ERROR: {e}"
    return result


# ---------------------------------------------------------------------------
# Terminal display
# ---------------------------------------------------------------------------

def display_table(descriptor_dict: dict):
    """Print a descriptor dict as an aligned terminal table."""
    if descriptor_dict is None:
        print("  (no data)")
        return
    print(f"\n{'Descriptor':<30} {'Value'}")
    print("-" * 60)
    for key, value in descriptor_dict.items():
        if isinstance(value, float):
            value = f"{value:.4f}"
        print(f"{key:<30} {value}")


# ---------------------------------------------------------------------------
# CSV writer helper
# ---------------------------------------------------------------------------

def _write_csv(rows: list[dict], output_path: str, total: int):
    """Write a list of dicts to CSV, building fieldnames from all row keys."""
    seen: dict[str, None] = {}
    for row in rows:
        seen.update(dict.fromkeys(row.keys()))
    fieldnames = list(seen.keys())
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nSaved {len(rows)}/{total} rows → {output_path}")


# ---------------------------------------------------------------------------
# Batch row resolver
# ---------------------------------------------------------------------------

def _resolve_batch_row(row: dict, row_index: int, total: int) -> dict | None:
    """
    Resolve one CSV row with CID, Name, and/or SMILES columns.

    Resolution priority: SMILES > CID > Name.
    Invalid fields are flagged and ignored. Mismatches are flagged.

    Returns dict with mol, smiles, cid, iupac_name, common_name — or None.
    """
    label = f"[{row_index}/{total}]"

    raw_cid    = (row.get("CID",    "") or "").strip()
    raw_name   = (row.get("Name",   "") or "").strip()
    raw_smiles = (row.get("SMILES", "") or "").strip()

    if not raw_cid and not raw_name and not raw_smiles:
        print(f"  {label} SKIP — no CID, Name, or SMILES found in row")
        return None

    valid: dict[str, dict] = {}

    if raw_cid:
        try:
            r = resolve_compound_input(raw_cid)
            if r["mol"] is not None:
                valid["CID"] = r
            else:
                print(f"  {label} FLAG — CID '{raw_cid}' gave no valid molecule; ignoring")
        except Exception as e:
            print(f"  {label} FLAG — CID '{raw_cid}' is invalid ({e}); ignoring")

    if raw_name:
        try:
            r = resolve_compound_input(raw_name)
            if r["mol"] is not None:
                valid["Name"] = r
            else:
                print(f"  {label} FLAG — Name '{raw_name}' gave no valid molecule; ignoring")
        except Exception as e:
            print(f"  {label} FLAG — Name '{raw_name}' is invalid ({e}); ignoring")

    if raw_smiles:
        if Chem.MolFromSmiles(raw_smiles) is None:
            print(f"  {label} FLAG — SMILES '{raw_smiles}' is invalid; ignoring")
        else:
            try:
                r = resolve_compound_input(raw_smiles)
                valid["SMILES"] = r
            except Exception as e:
                print(f"  {label} FLAG — SMILES '{raw_smiles}' could not be resolved ({e}); ignoring")

    if not valid:
        print(f"  {label} SKIP — no valid input remained after validation")
        return None

    cids_found = {k: v["cid"] for k, v in valid.items() if v["cid"] is not None}
    if len(set(cids_found.values())) > 1:
        details = ", ".join(f"{k}→CID {c}" for k, c in cids_found.items())
        print(f"  {label} MISMATCH — inputs point to different compounds: {details}")

    chosen_key = next(k for k in ("SMILES", "CID", "Name") if k in valid)
    chosen = valid[chosen_key]
    print(f"  {label} Using {chosen_key} as primary input (CID={chosen['cid'] or 'N/A'})")

    return {
        "mol":          chosen["mol"],
        "smiles":       chosen["smiles"] or "",
        "cid":          chosen["cid"],
        "iupac_name":   chosen["iupac_name"],
        "common_name":  chosen["common_name"],
        "input_cid":    raw_cid,
        "input_name":   raw_name,
        "input_smiles": raw_smiles,
        "chosen_key":   chosen_key,
    }


# ---------------------------------------------------------------------------
# Interactive viewer — Flask + 3Dmol.js
# ---------------------------------------------------------------------------

def interactive_explorer_cli():
    """
    interactive_explorer — launch the local Molecule Explorer web viewer.

    Opens a Flask server at http://localhost:5050. Type a CID, compound name,
    or SMILES into the search box to display the 3D structure and descriptors.
    Press Ctrl+C to stop.
    """
    try:
        from flask import Flask, jsonify, request, render_template
    except ImportError:
        print("Error: Flask is not installed. Run: pip install flask")
        sys.exit(1)

    _pkg_dir = os.path.dirname(__file__)
    app = Flask(
        __name__,
        static_folder=os.path.join(_pkg_dir, "static"),
        template_folder=os.path.join(_pkg_dir, "templates"),
    )
    app.logger.disabled = True
    import logging
    logging.getLogger("werkzeug").setLevel(logging.ERROR)

    @app.route("/")
    def index():
        return render_template("viewer.html")

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
        sdf = None
        if cid:
            resp = requests.get(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
                "/SDF?record_type=3d",
                timeout=10,
            )
            if resp.status_code == 200:
                sdf = resp.text
        if not sdf:
            return jsonify({"error": (
                f"No 3D conformer available for this compound (CID={cid or 'N/A'}). "
                "PubChem 3D conformers are only available for compounds in their database."
            )})

        mol = resolved["mol"]
        viewer_props = _build_id_block(resolved)
        for name, func in VIEWER_DESCRIPTORS.items():
            try:
                v = func(mol)
                viewer_props[name] = round(v, 4) if isinstance(v, float) else v
            except Exception:
                pass

        return jsonify({"sdf": sdf, "props": viewer_props})

    port = 5050
    print(f"\n  Molecule Explorer running at: http://localhost:{port}")
    print("  WSL users: open that URL in your Windows browser.")
    print("  Press Ctrl+C to stop.\n")
    app.run(host="0.0.0.0", port=port, debug=False)


# ---------------------------------------------------------------------------
# CLI: rdkit_descriptor
# ---------------------------------------------------------------------------

def rdkit_descriptor_cli():
    """
    rdkit_descriptor — compute RDKit descriptors for one compound or a CSV batch.

    Single compound
    ---------------
      rdkit_descriptor -i <cid|name|smiles>
          Compute DEFAULT_DESCRIPTORS (identity + drug-likeness set).

      rdkit_descriptor -i aspirin -d TPSA MolLogP NumHDonors
          Compute specific descriptors by RDKit name.

      rdkit_descriptor -i aspirin --full
          Compute all ~210 RDKit descriptors.

      rdkit_descriptor -i aspirin --full --filter TPSA MolLogP MolWt
          All ~210, filtered to a subset.

      rdkit_descriptor -i aspirin -o results.csv
          Save any of the above to CSV.

    Batch
    -----
      rdkit_descriptor --batch molecules.csv -o results.csv
          Process a CSV with CID, Name, and/or SMILES columns.

      rdkit_descriptor --batch molecules.csv -o results.csv --full
      rdkit_descriptor --batch molecules.csv -o results.csv --filter TPSA MolLogP
      rdkit_descriptor --batch molecules.csv -o results.csv -d MolLogP NumHeteroatoms

    Descriptor names use RDKit's CalcMolDescriptors naming convention
    (e.g. TPSA, MolLogP, MolWt, NumHDonors, NumHAcceptors, RingCount, ...).
    Run: python -c "from rdkit.Chem import Descriptors; print([d[0] for d in Descriptors.descList])"
    to list all available names.
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="rdkit_descriptor",
        description="Compute RDKit descriptors for one compound or a CSV batch",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Descriptor names follow RDKit's CalcMolDescriptors convention:\n"
            "  TPSA, MolLogP, MolWt, NumHDonors, NumHAcceptors, RingCount, ...\n"
            "  (case-insensitive; ~210 available)\n\n"
            "Single compound examples:\n"
            "  rdkit_descriptor -i aspirin\n"
            "  rdkit_descriptor -i 3033 -d TPSA MolLogP\n"
            "  rdkit_descriptor -i aspirin --full\n"
            "  rdkit_descriptor -i aspirin --full --filter TPSA MolLogP MolWt\n"
            "  rdkit_descriptor -i aspirin -o results.csv\n\n"
            "Batch examples:\n"
            "  rdkit_descriptor --batch molecules.csv -o results.csv\n"
            "  rdkit_descriptor --batch molecules.csv -o results.csv --full\n"
            "  rdkit_descriptor --batch molecules.csv -o results.csv -d MolLogP"
        ),
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-i", "--input", metavar="COMPOUND",
        help="PubChem CID, compound name, or SMILES string",
    )
    input_group.add_argument(
        "--batch", metavar="CSV",
        help="Path to input CSV with CID, Name, and/or SMILES columns",
    )
    parser.add_argument(
        "-o", "--output", default=None, metavar="FILE",
        help="Output CSV file path (required in batch mode)",
    )
    parser.add_argument(
        "--full", action="store_true",
        help="Compute all ~210 RDKit descriptors instead of the default set",
    )
    parser.add_argument(
        "--filter", nargs="+", default=None, metavar="NAME",
        help="With --full: return only these descriptor names (case-insensitive)",
    )
    parser.add_argument(
        "-d", "--descriptors", nargs="+", default=None, metavar="NAME",
        help=(
            "Specific descriptor names to compute (RDKit naming, case-insensitive). "
            "Overrides default set. Ignored when --full is used."
        ),
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Single compound
    # ------------------------------------------------------------------
    if args.input:
        result = get_rdkit_dict(
            args.input,
            full=args.full,
            descriptors=args.descriptors,
            filter_keys=args.filter,
        )
        if result is None:
            sys.exit(1)
        display_table(result)
        if args.output:
            _write_csv([result], args.output, total=1)
        return

    # ------------------------------------------------------------------
    # Batch
    # ------------------------------------------------------------------
    if not os.path.isfile(args.batch):
        print(f"Error: input file not found: {args.batch}")
        sys.exit(1)
    if not args.output:
        print("Error: -o / --output is required in batch mode.")
        sys.exit(1)

    with open(args.batch, newline="") as fh:
        reader = csv.DictReader(fh)
        headers = list(reader.fieldnames or [])
        rows    = list(reader)

    if not {"CID", "Name", "SMILES"}.intersection(headers):
        print("Error: input CSV must contain at least one of: CID, Name, SMILES")
        sys.exit(1)

    output_rows: list[dict] = []
    total = len(rows)

    for i, row in enumerate(rows, 1):
        res = _resolve_batch_row(row, i, total)
        if res is None:
            continue

        mol         = res["mol"]
        cid         = res["cid"]
        iupac_name  = res["iupac_name"]
        common_name = res["common_name"]
        smiles      = res["smiles"]

        # Pass through all original columns
        out_row: dict = {col: row.get(col, "") for col in headers}

        # Append resolved identifiers
        out_row["CID_resolved"]         = cid if cid else ""
        out_row["Name_IUPAC_resolved"]  = iupac_name
        out_row["Name_common_resolved"] = common_name
        out_row["SMILES_resolved"]      = smiles

        # Compute descriptors
        if args.descriptors:
            all_desc = get_all_rdkit_descriptors(mol)
            for name in args.descriptors:
                canonical, value = _find_in_full(all_desc, name)
                if canonical:
                    out_row[canonical] = value
                else:
                    print(f"  Warning: descriptor '{name}' not found — skipping.")
        elif args.full:
            all_desc = get_all_rdkit_descriptors(mol)
            if args.filter:
                lower_map = {k.lower(): k for k in all_desc}
                all_desc  = {
                    lower_map[f.lower()]: all_desc[lower_map[f.lower()]]
                    for f in args.filter if f.lower() in lower_map
                }
            out_row.update(all_desc)
        else:
            for name, func in DEFAULT_DESCRIPTORS.items():
                try:
                    out_row[name] = func(mol)
                except Exception as e:
                    out_row[name] = f"ERROR: {e}"

        output_rows.append(out_row)
        print(f"    → done")

    if not output_rows:
        print("No rows processed — output file not written.")
        return

    _write_csv(output_rows, args.output, total)


# ---------------------------------------------------------------------------
# CLI: fetch_xyz_batch
# ---------------------------------------------------------------------------

def fetch_xyz_batch_cli():
    """
    fetch_xyz_batch — batch-generate 3D structure files from a CSV.

    Usage
    -----
      fetch_xyz_batch -i molecules.csv -o ./structures/
      fetch_xyz_batch -i molecules.csv -o ./structures/ --format xyz sdf
      fetch_xyz_batch -i molecules.csv -o ./out/ --format sdf pdb --code-col Abbreviation

    Input CSV must contain at least one of: CID, Name, SMILES.
    Resolution priority: SMILES > CID > Name.

    File naming priority
    --------------------
    1. --code-col value (e.g. Abbreviation column)
    2. CID_CommonName  (e.g. 3033_Diclofenac)
    3. CID only
    4. Common name only
    5. missing_id<N>  (valid SMILES with no PubChem match and no name)
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="fetch_xyz_batch",
        description="Batch-generate 3D structure files from a CSV",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  fetch_xyz_batch -i molecules.csv -o ./structures/\n"
            "  fetch_xyz_batch -i molecules.csv -o ./structures/ --format xyz sdf\n"
            "  fetch_xyz_batch -i molecules.csv -o ./out/ --code-col Abbreviation"
        ),
    )
    parser.add_argument(
        "-i", "--input", required=True, metavar="CSV",
        help="Input CSV with CID, Name, and/or SMILES columns",
    )
    parser.add_argument(
        "-o", "--output", required=True, metavar="DIR",
        help="Output directory for structure files",
    )
    parser.add_argument(
        "--format", "-f",
        nargs="+",
        choices=["xyz", "sdf", "pdb", "mol"],
        default=["xyz"],
        metavar="FORMAT",
        help="Output format(s): xyz sdf pdb mol  (default: xyz)",
    )
    parser.add_argument(
        "--code-col", default=None, metavar="COLUMN",
        help="CSV column to use as filename stem (e.g. Abbreviation)",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: input file not found: {args.input}")
        sys.exit(1)

    os.makedirs(args.output, exist_ok=True)

    with open(args.input, newline="") as fh:
        reader = csv.DictReader(fh)
        headers = list(reader.fieldnames or [])
        rows    = list(reader)

    if not {"CID", "Name", "SMILES"}.intersection(headers):
        print("Error: input CSV must contain at least one of: CID, Name, SMILES")
        sys.exit(1)

    from mof_toolkit.molecule_manager import get_3d_structure

    total           = len(rows)
    missing_counter = 0

    for i, row in enumerate(rows, 1):
        res = _resolve_batch_row(row, i, total)
        if res is None:
            continue

        code_val = (row.get(args.code_col, "") or "").strip() if args.code_col else ""
        if code_val:
            stem = code_val
        else:
            cid         = res["cid"]
            common_name = res["common_name"] or res["iupac_name"]
            safe_name   = (
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
            output_stem=os.path.join(args.output, stem),
            source="pubchem" if res["cid"] else "rdkit",
            cid=res["cid"],
        )
