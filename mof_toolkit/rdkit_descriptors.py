"""
rdkit_descriptors.py
--------------------
Compound resolution, RDKit descriptor computation, and CLI tools for the
MOF-Guest Toolkit.

Supported input types for all public functions and CLI commands:
    CID       — PubChem Compound ID (integer or numeric string)
    name      — IUPAC or common compound name (e.g. 'aspirin', 'cannabidiol')
    formula   — Molecular formula (e.g. 'CH4', 'c6h6', 'NaHCO3')
    smiles    — SMILES string (e.g. 'CC(=O)Oc1ccccc1C(=O)O')
    inchikey  — Standard InChIKey (27-char, e.g. 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N')

Python helpers
--------------
    resolve_compound_input(query, input_type)
        Resolve any supported identifier → {cid, iupac_name, common_name,
        smiles, formula, inchikey, mol}.  Always returns canonical SMILES.

    get_all_rdkit_descriptors(mol)
        All ~210 RDKit descriptors for an RDKit Mol.  NaN-safe, returns flat dict.

    get_descriptor(smiles, descriptor_names)
        Compute one or more RDKit descriptors from a SMILES string.
        Accepts RDKit CalcMolDescriptors names (case-insensitive).
        Returns a dict {descriptor_name: value, ...}.

    get_smiles(query, input_type)
        Convenience wrapper — resolve any identifier to canonical SMILES.

    get_rdkit_dict(query, full, descriptors, filter_keys)
        Full descriptor dict for one compound; supports default / full / custom sets.

    display_table(descriptor_dict)
        Print a descriptor dict as an aligned terminal table.

CLI commands
------------
    mol_explorer -input <query> [-format <type>] [-descriptor <...>] [-output <file>]
        Resolve a compound and print / save identifiers or descriptors.

    mol_explorer -batch <csv> [-descriptor <...>] [-output <file>]
        Batch-compute descriptors for every row of a CSV file.

    interactive_explorer
        Launch the local Molecule Explorer web viewer (Flask + 3Dmol.js).
"""

from __future__ import annotations

import csv
import math
import os
import re
import sys
import time

import pubchempy as pcp
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Fragments, rdMolDescriptors

# Silence RDKit's internal stderr SMILES-parse warnings — we handle errors ourselves
RDLogger.logger().setLevel(RDLogger.CRITICAL)


# ---------------------------------------------------------------------------
# Descriptor registries
# ---------------------------------------------------------------------------

# DEFAULT_DESCRIPTORS — the concise drug-likeness set shown when no flags are given.
# All use RDKit's CalcMolDescriptors naming where a direct equivalent exists.
DEFAULT_DESCRIPTORS: dict[str, object] = {
    "MolecularFormula":    rdMolDescriptors.CalcMolFormula,
    "ExactMolWt":          rdMolDescriptors.CalcExactMolWt,
    "NumHeavyAtoms":       lambda mol: mol.GetNumHeavyAtoms(),
    "RingCount":           rdMolDescriptors.CalcNumRings,
    "NumAromaticRings":    rdMolDescriptors.CalcNumAromaticRings,
    "NumHBA":              rdMolDescriptors.CalcNumHBA,
    "NumHBD":              rdMolDescriptors.CalcNumHBD,
    "NumRotatableBonds":   rdMolDescriptors.CalcNumRotatableBonds,
    "TPSA":                rdMolDescriptors.CalcTPSA,
}

# VIEWER_DESCRIPTORS — subset shown in the interactive web viewer.
# Identical to DEFAULT_DESCRIPTORS here; customise independently if needed.
VIEWER_DESCRIPTORS: dict[str, object] = DEFAULT_DESCRIPTORS


# ---------------------------------------------------------------------------
# EBI ChEMBL cross-reference lookup
# ---------------------------------------------------------------------------

def fetch_ebi_sources(inchikey: str) -> list[dict]:
    """
    Query the EBI Unichem API for all database cross-references for an InChIKey.

    Parameters
    ----------
    inchikey : str
        Standard 27-character InChIKey.

    Returns
    -------
    list of dicts with keys 'src_id' and 'src_compound_id'.
    Returns an empty list on any network or parse error.
    """
    try:
        url = f"https://www.ebi.ac.uk/unichem/rest/inchikey/{inchikey}"
        resp = requests.get(url, timeout=10)
        if resp.status_code != 200:
            return []
        data = resp.json()
        return data if isinstance(data, list) else []
    except Exception:
        return []


# ---------------------------------------------------------------------------
# PubChem helpers
# ---------------------------------------------------------------------------

def _pubchem_get(url: str, timeout: int = 10) -> requests.Response | None:
    """GET a PubChem REST URL, returning None on network / status errors."""
    try:
        resp = requests.get(url, timeout=timeout)
        return resp if resp.status_code == 200 else None
    except Exception:
        return None


def fetch_smiles_from_cid(cid: int) -> str | None:
    """
    Fetch the canonical SMILES for a PubChem CID.
    Returns None on failure.
    """
    resp = _pubchem_get(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
        "/property/CanonicalSMILES/TXT"
    )
    return resp.text.strip() if resp else None


def fetch_pubchem_metadata(cid: int) -> dict:
    """
    Fetch rich metadata for a PubChem CID: IUPAC name, common name, canonical
    SMILES, molecular formula, molecular weight, and InChIKey.

    Returns
    -------
    dict with keys: CID, IUPAC_Name, Common_Name, SMILES, MolecularFormula,
    MolecularWeight_PubChem, InChIKey.

    Raises requests.HTTPError on non-200 responses.
    """
    prop_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
        "/property/IUPACName,MolecularWeight,MolecularFormula,"
        "CanonicalSMILES,InChIKey/JSON"
    )
    resp = requests.get(prop_url, timeout=10)
    resp.raise_for_status()
    p = resp.json()["PropertyTable"]["Properties"][0]

    # Best-effort common name from synonyms list
    common_name = ""
    try:
        syn_resp = _pubchem_get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
        )
        if syn_resp:
            synonyms = (
                syn_resp.json()["InformationList"]["Information"][0].get("Synonym", [])
            )
            for s in synonyms:
                # Prefer human-readable names: skip pure digits and registry numbers
                if not s.isdigit() and not re.match(r"^\d{2,}-\d{2}-\d$", s):
                    common_name = s
                    break
            if not common_name and synonyms:
                common_name = synonyms[0]
    except Exception:
        pass

    return {
        "CID":                     int(cid),
        "IUPAC_Name":              p.get("IUPACName", ""),
        "Common_Name":             common_name,
        "SMILES":                  p.get("CanonicalSMILES", ""),
        "MolecularFormula":        p.get("MolecularFormula", ""),
        "MolecularWeight_PubChem": p.get("MolecularWeight", ""),
        "InChIKey":                p.get("InChIKey", ""),
    }


def _cid_from_smiles(smiles: str) -> int | None:
    """Query PubChem for the CID matching a SMILES string. Returns None on failure."""
    try:
        resp = requests.get(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/TXT",
            params={"smiles": smiles},
            timeout=10,
        )
        if resp.status_code == 200 and resp.text.strip():
            return int(resp.text.strip().splitlines()[0])
    except Exception:
        pass
    return None


def _cid_from_inchikey(inchikey: str) -> int | None:
    """Query PubChem for the CID matching an InChIKey. Returns None on failure."""
    try:
        resp = _pubchem_get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
            f"{inchikey}/cids/TXT"
        )
        if resp and resp.text.strip():
            return int(resp.text.strip().splitlines()[0])
    except Exception:
        pass
    return None


def _all_cids_from_name(name: str) -> list[int]:
    """
    Return all CIDs PubChem associates with a name (up to 100).
    Returns an empty list if nothing is found.
    """
    try:
        compounds = pcp.get_compounds(name, "name")
        return [int(c.cid) for c in compounds] if compounds else []
    except Exception:
        return []


def _all_cids_from_formula(formula: str) -> list[int]:
    """
    Return all CIDs PubChem associates with a molecular formula.
    Returns an empty list on failure.
    """
    try:
        compounds = pcp.get_compounds(formula, "formula")
        return sorted(int(c.cid) for c in compounds) if compounds else []
    except Exception:
        return []


# ---------------------------------------------------------------------------
# Formula normalisation helpers
# ---------------------------------------------------------------------------

# Known 2-letter element symbols (Hill notation: first letter upper, second lower)
_ELEMENTS_2: set[str] = {
    "He","Li","Be","Ne","Na","Mg","Al","Si","Cl","Ar","Ca","Sc","Ti","Cr","Mn",
    "Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Zr","Nb",
    "Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","Xe","Cs","Ba","La",
    "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",
    "Ta","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra",
    "Ac","Th","Pa","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
}


def _looks_like_formula(q: str) -> bool:
    """
    Heuristic: return True if *q* looks like a molecular formula.
    Requires at least one digit (to avoid bare element symbols like 'CO' or 'NO'
    falling through to this branch instead of the name-lookup path).
    """
    return (
        any(c.isdigit() for c in q)
        and bool(re.match(r"^([A-Z][a-z]?\d*)+$", q, re.IGNORECASE))
    )


def _looks_like_inchikey(q: str) -> bool:
    """Return True if *q* matches the standard InChIKey pattern (27 chars)."""
    return bool(re.match(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$", q))


def normalise_formula(raw: str) -> str:
    """
    Normalise a user-supplied molecular formula to proper Hill notation.

    PubChem's formula endpoint is case-sensitive; this handles:
        'ch4'     → 'CH4'
        'c6h6'    → 'C6H6'
        'nahco3'  → 'NaHCO3'
        'MGSO4'   → 'MgSO4'
        'CaCl2'   → 'CaCl2'   (already correct)
    """
    has_upper = any(c.isupper() for c in raw if c.isalpha())
    has_lower = any(c.islower() for c in raw if c.isalpha())
    mixed = has_upper and has_lower

    result: list[str] = []
    i = 0
    while i < len(raw):
        c = raw[i]
        if c.isalpha():
            if mixed:
                # Mixed case: 2-letter element only when first=upper AND second=lower
                if (
                    c.isupper()
                    and i + 1 < len(raw)
                    and raw[i + 1].isalpha()
                    and raw[i + 1].islower()
                ):
                    two = c + raw[i + 1]
                    if two in _ELEMENTS_2:
                        result.append(two)
                        i += 2
                        continue
            else:
                # Uniform case: try greedy 2-letter match
                if i + 1 < len(raw) and raw[i + 1].isalpha():
                    two = c.upper() + raw[i + 1].lower()
                    if two in _ELEMENTS_2:
                        result.append(two)
                        i += 2
                        continue
            result.append(c.upper())
            i += 1
        else:
            result.append(c)
            i += 1
    return "".join(result)


# ---------------------------------------------------------------------------
# Random compound picker
# ---------------------------------------------------------------------------

#: Upper bound for random CID sampling (covers ~119 M compounds as of 2025).
_PUBCHEM_MAX_CID: int = 150_000_000


def get_random_cid(max_attempts: int = 20) -> int:
    """
    Return a random PubChem CID that resolves to a real compound.

    Samples random integers in [1, _PUBCHEM_MAX_CID] and validates each
    candidate with a lightweight single-property fetch.  The sampling is
    intentionally biased toward lower CIDs (more common compounds) via
    exponential weighting on ~40 % of attempts.

    Raises RuntimeError if no valid CID is found within *max_attempts*.
    """
    import random

    for attempt in range(max_attempts):
        # Bias the first half of attempts toward smaller (more common) CIDs
        if attempt < max_attempts // 2:
            candidate = random.randint(1, 5_000_000)
        else:
            candidate = random.randint(1, _PUBCHEM_MAX_CID)

        resp = _pubchem_get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{candidate}"
            "/property/CanonicalSMILES/TXT",
            timeout=8,
        )
        if resp and resp.text.strip():
            return candidate

    raise RuntimeError(
        f"Could not find a valid random compound after {max_attempts} attempts."
    )


# ---------------------------------------------------------------------------
# Core resolver
# ---------------------------------------------------------------------------

def resolve_compound_input(query: str, input_type: str = "auto") -> dict:
    """
    Resolve any supported compound identifier to a unified metadata dict.

    Parameters
    ----------
    query : str
        The identifier to resolve.  May be a PubChem CID, compound name,
        molecular formula, SMILES string, or InChIKey.
    input_type : str, optional
        How to interpret *query*.  One of:
          'auto'     — detect automatically (default)
          'cid'      — integer PubChem CID
          'name'     — IUPAC or common name
          'formula'  — molecular formula (e.g. 'CH4', 'c6h6')
          'smiles'   — SMILES string
          'inchikey' — 27-character standard InChIKey

    Returns
    -------
    dict with keys:
        cid         (int | None)
        iupac_name  (str)
        common_name (str)
        smiles      (str | None)   — canonical SMILES
        formula     (str)
        inchikey    (str)
        mol         (Chem.Mol | None)
        all_cids    (list[int])    — all CIDs found (only populated for name/formula)

    Notes
    -----
    * Auto-detection order: CID → InChIKey → formula → SMILES → name.
    * When the SMILES is valid but not in PubChem, cid/names/inchikey are
      empty but mol and smiles are populated for local descriptor computation.
    * Always returns canonical SMILES (RDKit Chem.MolToSmiles).
    * When a name or formula maps to multiple CIDs, the first (lowest) CID
      is used and all_cids is populated.  The caller should warn the user.
    """
    q = query.strip()

    # Shared empty-result skeleton
    def _empty() -> dict:
        return {
            "cid": None, "iupac_name": "", "common_name": "",
            "smiles": None, "formula": "", "inchikey": "", "mol": None,
            "all_cids": [],
        }

    def _from_cid(cid: int) -> dict:
        """Build a full result dict from a known CID."""
        smiles_raw = fetch_smiles_from_cid(cid)
        mol = Chem.MolFromSmiles(smiles_raw) if smiles_raw else None
        canonical = Chem.MolToSmiles(mol) if mol else smiles_raw
        try:
            meta = fetch_pubchem_metadata(cid)
            iupac   = meta.get("IUPAC_Name", "")
            common  = meta.get("Common_Name", "")
            formula = meta.get("MolecularFormula", "")
            inchkey = meta.get("InChIKey", "")
        except Exception:
            iupac = common = formula = inchkey = ""
        return {
            "cid": cid, "iupac_name": iupac, "common_name": common,
            "smiles": canonical, "formula": formula, "inchikey": inchkey,
            "mol": mol, "all_cids": [cid],
        }

    # ------------------------------------------------------------------ CID
    if input_type == "cid" or (input_type == "auto" and q.isdigit()):
        return _from_cid(int(q))

    # ------------------------------------------------------------------ InChIKey
    if input_type == "inchikey" or (input_type == "auto" and _looks_like_inchikey(q)):
        cid = _cid_from_inchikey(q)
        if cid is None:
            if input_type == "inchikey":
                raise ValueError(f"No PubChem compound found for InChIKey: '{q}'")
            # Fall through in auto mode
        else:
            return _from_cid(cid)

    # ------------------------------------------------------------------ Formula
    if input_type == "formula" or (input_type == "auto" and _looks_like_formula(q)):
        norm = normalise_formula(q)
        all_cids = _all_cids_from_formula(norm)
        if not all_cids:
            if input_type == "formula":
                raise ValueError(f"No PubChem compound found for formula: '{q}'")
            # Fall through in auto mode
        else:
            result = _from_cid(all_cids[0])
            result["all_cids"] = all_cids
            return result

    # ------------------------------------------------------------------ SMILES
    if input_type == "smiles" or input_type == "auto":
        mol_test = Chem.MolFromSmiles(q)
        if mol_test is not None:
            canonical = Chem.MolToSmiles(mol_test)
            cid = _cid_from_smiles(canonical)
            if cid is None:
                print(
                    f"  Note: SMILES '{q}' is valid but not found in PubChem — "
                    "CID/IUPAC/formula/InChIKey will be empty."
                )
                formula_local = rdMolDescriptors.CalcMolFormula(mol_test)
                return {
                    "cid": None, "iupac_name": "", "common_name": "",
                    "smiles": canonical, "formula": formula_local,
                    "inchikey": "", "mol": mol_test, "all_cids": [],
                }
            result = _from_cid(cid)
            # Preserve the user-supplied canonical SMILES to avoid confusion
            result["smiles"] = canonical
            return result
        if input_type == "smiles":
            raise ValueError(f"Invalid SMILES string: '{q}'")

    # ------------------------------------------------------------------ Name
    all_cids = _all_cids_from_name(q)
    if not all_cids:
        raise ValueError(f"No PubChem compound found for name: '{q}'")
    result = _from_cid(all_cids[0])
    result["all_cids"] = all_cids
    return result


# ---------------------------------------------------------------------------
# Descriptor helpers
# ---------------------------------------------------------------------------

def get_all_rdkit_descriptors(mol: Chem.Mol) -> dict:
    """
    Compute all ~210 RDKit descriptors for an RDKit Mol object.

    Uses RDKit's CalcMolDescriptors naming (e.g. 'TPSA', 'MolLogP',
    'NumHDonors').  NaN values are replaced with empty strings so the
    result is always CSV-safe.

    Parameters
    ----------
    mol : Chem.Mol

    Returns
    -------
    dict mapping RDKit descriptor name → value (float, int, or '').

    Examples
    --------
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


def get_descriptor(smiles: str, descriptor_names: list[str] | str) -> dict:
    """
    Compute one or more RDKit descriptors from a SMILES string.

    Accepts any name from RDKit's CalcMolDescriptors registry (case-insensitive).
    This is the low-level workhorse used by get_rdkit_dict and the CLI; it is
    also useful for quick one-off calculations in scripts and notebooks.

    Parameters
    ----------
    smiles : str
        Valid canonical SMILES string.
    descriptor_names : str or list of str
        One or more RDKit descriptor names (e.g. 'TPSA', 'MolLogP',
        ['TPSA', 'MolLogP', 'NumHDonors']).  Case-insensitive.

    Returns
    -------
    dict mapping canonical descriptor name → value.
    Unrecognised names are warned and excluded from the result.

    Raises
    ------
    ValueError
        If *smiles* is not a valid SMILES string.

    Examples
    --------
    >>> get_descriptor("CC(=O)O", "TPSA")
    {'TPSA': 37.3}
    >>> get_descriptor("CC(=O)O", ["TPSA", "MolLogP", "NumHDonors"])
    {'TPSA': 37.3, 'MolLogP': -0.17, 'NumHDonors': 1}
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: '{smiles}'")

    if isinstance(descriptor_names, str):
        descriptor_names = [descriptor_names]

    all_desc = get_all_rdkit_descriptors(mol)
    lower_map = {k.lower(): k for k in all_desc}
    result: dict = {}

    for name in descriptor_names:
        canonical = lower_map.get(name.lower())
        if canonical is not None:
            result[canonical] = all_desc[canonical]
        else:
            print(
                f"  Warning: descriptor '{name}' not found in RDKit. "
                "Run get_all_rdkit_descriptors(mol).keys() to see all ~210 names."
            )
    return result


def get_smiles(query: str, input_type: str = "auto") -> str | None:
    """
    Resolve any supported identifier to its canonical SMILES string.

    Parameters
    ----------
    query : str
        CID, name, formula, SMILES, or InChIKey.
    input_type : str, optional
        Same values as resolve_compound_input.  Default: 'auto'.

    Returns
    -------
    Canonical SMILES string, or None if resolution fails.

    Examples
    --------
    >>> get_smiles(3033)             # CID
    'OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl'
    >>> get_smiles("aspirin")        # name
    'CC(=O)Oc1ccccc1C(=O)O'
    >>> get_smiles("C9H8O4", "formula")
    'CC(=O)Oc1ccccc1C(=O)O'
    """
    try:
        res = resolve_compound_input(str(query), input_type=input_type)
        return res["smiles"]
    except Exception as e:
        print(f"  Error resolving '{query}': {e}")
        return None


# ---------------------------------------------------------------------------
# Shared output helpers
# ---------------------------------------------------------------------------

def _build_id_block(resolved: dict) -> dict:
    """
    Build the standard CID/name/SMILES/formula/InChIKey identifier prefix.
    Used as the leading columns in all descriptor outputs.
    """
    return {
        "CID":         resolved["cid"] if resolved["cid"] else "",
        "IUPAC_Name":  resolved["iupac_name"],
        "Common_Name": resolved["common_name"],
        "SMILES":      resolved["smiles"] or "",
        "Formula":     resolved["formula"],
        "InChIKey":    resolved["inchikey"],
    }


def _find_descriptor_in_full(all_desc: dict, name: str) -> tuple[str | None, object]:
    """
    Look up a descriptor name in get_all_rdkit_descriptors output.
    Tries exact match first, then case-insensitive.
    Returns (canonical_name, value) or (None, None) if not found.
    """
    if name in all_desc:
        return name, all_desc[name]
    lower_map = {k.lower(): k for k in all_desc}
    canon = lower_map.get(name.lower())
    if canon:
        return canon, all_desc[canon]
    return None, None


# ---------------------------------------------------------------------------
# get_rdkit_dict — main Python helper for single-compound descriptor fetching
# ---------------------------------------------------------------------------

def get_rdkit_dict(
    query,
    full: bool = False,
    descriptors: list[str] | None = None,
    filter_keys: list[str] | None = None,
    input_type: str = "auto",
) -> dict | None:
    """
    Return a descriptor dict for a single compound.

    Accepts a PubChem CID (int or numeric string), a compound name, a
    molecular formula, a SMILES string, or an InChIKey.  When the SMILES
    is valid but not in PubChem, CID/names are empty but all descriptors
    are still computed from the SMILES.

    Parameters
    ----------
    query : int or str
        Compound identifier (CID, name, formula, SMILES, or InChIKey).
    full : bool, optional
        Compute all ~210 RDKit descriptors using RDKit's own naming.
        Default: False (uses DEFAULT_DESCRIPTORS).
    descriptors : list of str, optional
        Specific RDKit descriptor names to compute (case-insensitive).
        Overrides DEFAULT_DESCRIPTORS when provided.
    filter_keys : list of str, optional
        When full=True, keep only descriptors matching these names.
        Case-insensitive.  Ignored when full=False.
    input_type : str, optional
        Hint for how to interpret *query* (default: 'auto').

    Returns
    -------
    dict with keys CID, IUPAC_Name, Common_Name, SMILES, Formula,
    InChIKey, and descriptor values; or None on failure.

    Examples
    --------
    >>> from mof_toolkit.rdkit_descriptors import get_rdkit_dict, display_table

    # Default set
    >>> display_table(get_rdkit_dict(3033))
    >>> display_table(get_rdkit_dict("aspirin"))
    >>> display_table(get_rdkit_dict("CC(=O)Oc1ccccc1C(=O)O"))

    # Specific descriptors
    >>> display_table(get_rdkit_dict("aspirin", descriptors=["TPSA", "MolLogP"]))

    # Full ~210 descriptors
    >>> display_table(get_rdkit_dict("aspirin", full=True))

    # Full set, filtered subset
    >>> display_table(get_rdkit_dict("aspirin", full=True,
    ...                              filter_keys=["TPSA", "MolLogP", "MolWt"]))
    """
    try:
        resolved = resolve_compound_input(str(query), input_type=input_type)
    except Exception as e:
        print(f"  Error resolving '{query}': {e}")
        return None

    mol = resolved["mol"]
    if mol is None:
        print(f"  Error: could not obtain a valid molecule for '{query}'.")
        return None

    id_block = _build_id_block(resolved)

    # Warn when the identifier mapped to multiple PubChem entries
    all_cids = resolved.get("all_cids", [])
    if len(all_cids) > 1:
        shown = all_cids[1:11]
        suffix = " … and more!" if len(all_cids) > 11 else ""
        others = ", ".join(str(c) for c in shown) + suffix
        print(
            f"  Warning: '{query}' mapped to {len(all_cids)} PubChem entries. "
            f"Using CID {all_cids[0]}. Other CIDs: {others}"
        )

    # Case 1 — specific descriptor names
    if descriptors is not None:
        all_desc = get_all_rdkit_descriptors(mol)
        result = dict(id_block)
        for name in descriptors:
            canon, value = _find_descriptor_in_full(all_desc, name)
            if canon is not None:
                result[canon] = value
            else:
                print(
                    f"  Warning: descriptor '{name}' not found — skipping. "
                    "Run get_all_rdkit_descriptors(mol).keys() to see all ~210 names."
                )
        return result

    # Case 2 — full ~210 descriptors
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

    # Case 3 — default set
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

def display_table(descriptor_dict: dict | None) -> None:
    """Print a descriptor dict as an aligned, readable terminal table."""
    if descriptor_dict is None:
        print("  (no data)")
        return
    print(f"\n{'Descriptor':<32} Value")
    print("-" * 62)
    for key, value in descriptor_dict.items():
        formatted = f"{value:.4f}" if isinstance(value, float) else str(value)
        print(f"{key:<32} {formatted}")


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------

def _write_csv(rows: list[dict], output_path: str, total: int) -> None:
    """
    Write a list of dicts to a CSV file.
    Column order is determined by the union of all row keys (first-seen order).
    """
    seen: dict[str, None] = {}
    for row in rows:
        seen.update(dict.fromkeys(row.keys()))
    fieldnames = list(seen.keys())
    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nSaved {len(rows)}/{total} rows → {output_path}")


# ---------------------------------------------------------------------------
# Batch row resolver (shared by mol_explorer --batch and fetch_xyz_batch)
# ---------------------------------------------------------------------------

# Column names accepted as each input type
_BATCH_COL_ALIASES: dict[str, list[str]] = {
    "smiles":   ["smiles", "SMILES", "Smiles"],
    "cid":      ["cid", "CID", "Cid", "pubchem_cid", "PubChem_CID"],
    "name":     ["name", "Name", "NAME", "compound_name", "Compound_Name"],
    "inchikey": ["inchikey", "InChIKey", "InChI_Key", "INCHIKEY"],
}


def _get_batch_column(row: dict, field: str) -> str:
    """Return the first non-empty value from *row* for all known aliases of *field*."""
    for alias in _BATCH_COL_ALIASES.get(field, [field]):
        val = (row.get(alias) or "").strip()
        if val:
            return val
    return ""


def _resolve_batch_row(row: dict, row_index: int, total: int) -> dict | None:
    """
    Resolve one CSV row that may contain any combination of CID, Name,
    SMILES, and/or InChIKey columns.

    Resolution priority: SMILES > CID > InChIKey > Name.
    Invalid fields are flagged and skipped.  When multiple valid fields point
    to different PubChem compounds a MISMATCH warning is printed.

    If the row has a SMILES and other identifier columns, the SMILES is
    authoritative; resolved PubChem metadata (CID, name, etc.) is updated
    from the SMILES lookup and a warning is printed.

    Returns a dict with mol, smiles, cid, iupac_name, common_name, formula,
    inchikey, and bookkeeping fields — or None if no valid input remains.
    """
    label = f"[{row_index}/{total}]"

    raw = {k: _get_batch_column(row, k) for k in ("smiles", "cid", "inchikey", "name")}

    if not any(raw.values()):
        print(f"  {label} SKIP — no CID, Name, SMILES, or InChIKey found in row")
        return None

    valid: dict[str, dict] = {}

    # Validate each present field independently
    for field in ("smiles", "cid", "inchikey", "name"):
        val = raw[field]
        if not val:
            continue
        try:
            r = resolve_compound_input(val, input_type=field if field != "smiles" else "smiles")
            if r["mol"] is not None:
                valid[field] = r
            else:
                print(f"  {label} FLAG — {field.upper()} '{val}' gave no valid molecule; ignoring")
        except Exception as e:
            print(f"  {label} FLAG — {field.upper()} '{val}' is invalid ({e}); ignoring")

    if not valid:
        print(f"  {label} SKIP — no valid input remained after validation")
        return None

    # Check for CID mismatches between valid fields
    cids_found = {k: v["cid"] for k, v in valid.items() if v["cid"] is not None}
    unique_cids = set(cids_found.values())
    if len(unique_cids) > 1:
        details = ", ".join(f"{k.upper()}→CID {c}" for k, c in cids_found.items())
        print(f"  {label} MISMATCH — inputs point to different compounds: {details}")

    # Choose highest-priority valid field
    chosen_key = next(k for k in ("smiles", "cid", "inchikey", "name") if k in valid)
    chosen = valid[chosen_key]

    # Warn if SMILES caused a metadata update
    if chosen_key == "smiles" and len(valid) > 1:
        print(
            f"  {label} Note: SMILES is authoritative — "
            "other identifier fields have been updated from the SMILES lookup."
        )

    print(f"  {label} Using {chosen_key.upper()} as primary input (CID={chosen['cid'] or 'N/A'})")

    # Warn about multiple PubChem mappings
    all_cids = chosen.get("all_cids", [])
    if len(all_cids) > 1:
        shown = all_cids[1:11]
        suffix = " … and more!" if len(all_cids) > 11 else ""
        others = ", ".join(str(c) for c in shown) + suffix
        print(
            f"  {label} Warning: query mapped to {len(all_cids)} PubChem entries. "
            f"Using CID {all_cids[0]}. Other CIDs: {others}"
        )

    return {
        "mol":           chosen["mol"],
        "smiles":        chosen["smiles"] or "",
        "cid":           chosen["cid"],
        "iupac_name":    chosen["iupac_name"],
        "common_name":   chosen["common_name"],
        "formula":       chosen["formula"],
        "inchikey":      chosen["inchikey"],
        "input_smiles":  raw["smiles"],
        "input_cid":     raw["cid"],
        "input_inchikey":raw["inchikey"],
        "input_name":    raw["name"],
        "chosen_key":    chosen_key,
    }


# ---------------------------------------------------------------------------
# Interactive viewer — Flask + 3Dmol.js
# ---------------------------------------------------------------------------

def interactive_explorer_cli() -> None:
    """
    interactive_explorer — launch the local Molecule Explorer web viewer.

    Opens a Flask development server at http://localhost:5050.
    Users can search by CID, name, formula, SMILES, or InChIKey.
    The 3D structure is fetched from PubChem and rendered with 3Dmol.js.
    Press Ctrl+C to stop.
    """
    try:
        from flask import Flask, jsonify, request, render_template
    except ImportError:
        print("Error: Flask is not installed.  Run: pip install flask")
        sys.exit(1)

    _pkg_dir = os.path.dirname(__file__)
    app = Flask(
        __name__,
        static_folder=os.path.join(_pkg_dir, "static"),
        template_folder=os.path.join(_pkg_dir, "templates"),
    )
    # Silence Flask/Werkzeug request logs — only errors are shown
    app.logger.disabled = True
    import logging
    logging.getLogger("werkzeug").setLevel(logging.ERROR)

    @app.route("/")
    def index():
        return render_template("viewer.html")

    @app.route("/lucky")
    def lucky():
        """Pick a random valid PubChem compound and return the viewer payload."""
        try:
            cid = get_random_cid(max_attempts=20)
        except RuntimeError as e:
            return jsonify({"error": str(e)})

        try:
            resolved = resolve_compound_input(str(cid), input_type="cid")
        except Exception as e:
            return jsonify({"error": str(e)})

        if resolved["mol"] is None:
            return jsonify({"error": f"Could not resolve CID {cid}"})

        resp_3d = _pubchem_get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
            "/SDF?record_type=3d",
            timeout=12,
        )
        if resp_3d is None:
            return jsonify({"error": f"No 3D conformer available for CID {cid}."})

        mol = resolved["mol"]
        viewer_props = _build_id_block(resolved)
        for name, func in VIEWER_DESCRIPTORS.items():
            try:
                v = func(mol)
                viewer_props[name] = round(v, 4) if isinstance(v, float) else v
            except Exception:
                pass

        return jsonify({"sdf": resp_3d.text, "props": viewer_props})

    @app.route("/ebi_refs")
    def ebi_refs():
        """Return EBI Unichem cross-references for an InChIKey as JSON."""
        inchikey = request.args.get("inchikey", "").strip()
        if not inchikey:
            return jsonify([])
        sources = fetch_ebi_sources(inchikey)
        return jsonify(sources)

    @app.route("/lookup")
    def lookup():
        """Resolve a compound query and return SDF + descriptor payload for the viewer."""
        q          = request.args.get("q",    "").strip()
        input_type = request.args.get("type", "auto").strip().lower()

        if not q:
            return jsonify({"error": "empty query"})

        valid_types = {"auto", "cid", "name", "formula", "smiles", "inchikey"}
        if input_type not in valid_types:
            return jsonify({"error": f"Invalid input type '{input_type}'. "
                                     f"Choose from: {sorted(valid_types)}"})

        try:
            resolved = resolve_compound_input(q, input_type=input_type)
        except Exception as e:
            return jsonify({"error": str(e)})

        if resolved["mol"] is None:
            return jsonify({"error": f"Could not resolve a valid molecule for '{q}'"})

        # Warn about multiple mappings in the JSON payload
        warnings: list[str] = []
        all_cids = resolved.get("all_cids", [])
        if len(all_cids) > 1:
            shown = all_cids[1:11]
            suffix = " … and more!" if len(all_cids) > 11 else ""
            others = ", ".join(str(c) for c in shown) + suffix
            warnings.append(
                f"'{q}' matched {len(all_cids)} PubChem entries. "
                f"Showing CID {all_cids[0]}. Others: {others}"
            )

        cid = resolved["cid"]
        sdf = None
        if cid:
            resp = _pubchem_get(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
                "/SDF?record_type=3d",
                timeout=12,
            )
            if resp:
                sdf = resp.text

        if not sdf:
            return jsonify({"error": (
                f"No 3D conformer available for this compound "
                f"(CID={cid or 'N/A'}). "
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

        payload = {"sdf": sdf, "props": viewer_props}
        if warnings:
            payload["warnings"] = warnings
        return jsonify(payload)

    port = 5050
    print(f"\n  Molecule Explorer running at: http://localhost:{port}")
    print("  WSL users: open that URL in your Windows browser.")
    print("  Press Ctrl+C to stop.\n")
    app.run(host="0.0.0.0", port=port, debug=False)


# ---------------------------------------------------------------------------
# CLI: mol_explorer
# ---------------------------------------------------------------------------

def mol_explorer_cli() -> None:
    """
    mol_explorer — resolve compound identifiers and compute RDKit descriptors.

    Single-compound usage
    ---------------------
      mol_explorer -input <query> [-format <type>]
          Resolve a compound and print its CID, SMILES, name, formula.

      mol_explorer -input <query> -descriptor default
          Print the default descriptor set (default when -descriptor is omitted).

      mol_explorer -input <query> -descriptor full
          Print all ~210 RDKit descriptors.

      mol_explorer -input <query> -descriptor TPSA MolLogP NumHDonors
          Print specific descriptors by RDKit name.

      mol_explorer -input <query> -descriptor TPSA -output results.csv
          Save descriptors to a CSV or TXT file.

    Batch usage
    -----------
      mol_explorer -batch molecules.csv -output results.csv
          Process a CSV with CID/Name/SMILES/InChIKey columns.

      mol_explorer -batch molecules.csv -descriptor full -output results.csv

    Input types (for -format flag, optional — default: auto)
    ---------------------------------------------------------
      cid | name | formula | smiles | inchikey

    Examples
    --------
      mol_explorer -input 3033
      mol_explorer -input aspirin -format name
      mol_explorer -input "CC(=O)Oc1ccccc1C(=O)O" -format smiles
      mol_explorer -input BSYNRYMUTXBXSQ-UHFFFAOYSA-N -format inchikey
      mol_explorer -input aspirin -descriptor TPSA MolLogP
      mol_explorer -input aspirin -descriptor full -output aspirin_all.csv
      mol_explorer -batch mols.csv -descriptor default -output results.csv
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="mol_explorer",
        description="Resolve compound identifiers and compute RDKit descriptors",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Descriptor names follow RDKit's CalcMolDescriptors convention:\n"
            "  TPSA, MolLogP, MolWt, NumHDonors, NumHAcceptors, RingCount …\n"
            "  (case-insensitive; ~210 available)\n\n"
            "Single-compound examples:\n"
            "  mol_explorer -input 3033\n"
            "  mol_explorer -input aspirin -format name\n"
            "  mol_explorer -input aspirin -descriptor TPSA MolLogP\n"
            "  mol_explorer -input aspirin -descriptor full -output all.csv\n\n"
            "Batch examples:\n"
            "  mol_explorer -batch molecules.csv -output results.csv\n"
            "  mol_explorer -batch molecules.csv -descriptor full -output results.csv\n"
        ),
    )

    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        "-input", dest="input", metavar="COMPOUND",
        help="CID, name, formula, SMILES, or InChIKey",
    )
    mode_group.add_argument(
        "-batch", dest="batch", metavar="CSV",
        help="Input CSV with CID/Name/SMILES/InChIKey columns",
    )
    parser.add_argument(
        "-format", dest="fmt", metavar="TYPE",
        choices=["cid", "name", "formula", "smiles", "inchikey", "auto"],
        default="auto",
        help="Explicit input type (default: auto-detect)",
    )
    parser.add_argument(
        "-descriptor", dest="descriptor", nargs="+", metavar="NAME",
        default=None,
        help=(
            "Descriptors to compute: 'default', 'full', or specific RDKit names "
            "(e.g. TPSA MolLogP).  Default: print identifiers only."
        ),
    )
    parser.add_argument(
        "-output", dest="output", metavar="FILE", default=None,
        help="Output file path (.csv or .txt)",
    )
    args = parser.parse_args()

    # Validate format flag
    valid_formats = {"cid", "name", "formula", "smiles", "inchikey", "auto"}
    if args.fmt not in valid_formats:
        print(f"Error: invalid -format '{args.fmt}'. "
              f"Choose from: {', '.join(sorted(valid_formats))}")
        sys.exit(1)

    # Parse -descriptor into (full, specific_list) flags
    use_full       = False
    use_default    = False
    specific_descs = None

    if args.descriptor is not None:
        lowered = [d.lower() for d in args.descriptor]
        if "full" in lowered:
            use_full = True
        elif "default" in lowered:
            use_default = True
        else:
            specific_descs = args.descriptor

    # ------------------------------------------------------------------
    # Single-compound mode
    # ------------------------------------------------------------------
    if args.input:
        try:
            resolved = resolve_compound_input(args.input, input_type=args.fmt)
        except Exception as e:
            print(f"Error: invalid {args.fmt} '{args.input}': {e}")
            sys.exit(1)

        # Warn about multiple PubChem mappings
        all_cids = resolved.get("all_cids", [])
        if len(all_cids) > 1:
            shown = all_cids[1:11]
            suffix = " … and more!" if len(all_cids) > 11 else ""
            others = ", ".join(str(c) for c in shown) + suffix
            print(
                f"Warning: '{args.input}' matched {len(all_cids)} PubChem entries. "
                f"Using CID {all_cids[0]}. Other CIDs: {others}"
            )

        mol = resolved["mol"]

        # No -descriptor flag → print identifiers only
        if args.descriptor is None:
            id_block = _build_id_block(resolved)
            display_table(id_block)
            if args.output:
                _write_csv([id_block], args.output, total=1)
            return

        # Build the descriptor dict
        if mol is None:
            print("Error: could not obtain a valid molecule.")
            sys.exit(1)

        if use_full:
            result = get_rdkit_dict(args.input, full=True, input_type=args.fmt)
        elif specific_descs:
            result = get_rdkit_dict(args.input, descriptors=specific_descs,
                                    input_type=args.fmt)
        else:  # default
            result = get_rdkit_dict(args.input, input_type=args.fmt)

        if result is None:
            sys.exit(1)

        display_table(result)
        if args.output:
            _write_csv([result], args.output, total=1)
        return

    # ------------------------------------------------------------------
    # Batch mode
    # ------------------------------------------------------------------
    if not os.path.isfile(args.batch):
        print(f"Error: input file not found: {args.batch}")
        sys.exit(1)
    if not args.output:
        print("Error: -output is required in batch mode.")
        sys.exit(1)

    with open(args.batch, newline="", encoding="utf-8") as fh:
        reader  = csv.DictReader(fh)
        headers = list(reader.fieldnames or [])
        rows    = list(reader)

    # Check that at least one recognised column is present
    accepted_cols = {alias for aliases in _BATCH_COL_ALIASES.values() for alias in aliases}
    if not accepted_cols.intersection(headers):
        print(
            "Error: input CSV must contain at least one of: "
            "CID, Name, SMILES, InChIKey"
        )
        sys.exit(1)

    output_rows: list[dict] = []
    total = len(rows)

    for i, row in enumerate(rows, 1):
        res = _resolve_batch_row(row, i, total)
        if res is None:
            continue

        mol = res["mol"]

        # Start with all original columns, then append resolved identifiers
        out_row: dict = {col: row.get(col, "") for col in headers}
        out_row["CID_resolved"]         = res["cid"] if res["cid"] else ""
        out_row["IUPAC_Name_resolved"]  = res["iupac_name"]
        out_row["Common_Name_resolved"] = res["common_name"]
        out_row["SMILES_resolved"]      = res["smiles"]
        out_row["Formula_resolved"]     = res["formula"]
        out_row["InChIKey_resolved"]    = res["inchikey"]

        # Append descriptors
        if use_full:
            out_row.update(get_all_rdkit_descriptors(mol))
        elif specific_descs:
            all_desc = get_all_rdkit_descriptors(mol)
            for name in specific_descs:
                canon, value = _find_descriptor_in_full(all_desc, name)
                if canon:
                    out_row[canon] = value
                else:
                    print(f"  Warning: descriptor '{name}' not found — skipping.")
        elif use_default or args.descriptor is not None:
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