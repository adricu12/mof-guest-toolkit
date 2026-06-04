"""
molecule_manager.py
-------------------
Local molecule manipulation — 3D conformer generation via RDKit.
No PubChem dependency for the core logic; PubChem is used optionally
to download pre-computed 3D conformers when a CID is available.

Public Python helpers
---------------------
  embed_3d(mol)
      Add Hs, generate a 3D conformer (ETKDGv3), run MMFF94/UFF optimisation.
      Returns the 3D mol with Hs on success, raises RuntimeError on failure.

  get_3d_structure(query, formats, output_stem, source, cid)
      Obtain a 3D structure for a CID, compound name, or SMILES and write it
      in one or more formats.  Single-compound helper — call in a loop for lists.

CLI commands
------------
  smiles_to_3d <SMILES> [--format xyz sdf pdb mol] [--output stem]
      Convert a SMILES string to 3D coordinate file(s).
"""

import os
import sys

import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


# ---------------------------------------------------------------------------
# Low-level format writers
# ---------------------------------------------------------------------------

def _write_xyz(mol: Chem.Mol, filepath: str, title: str = "") -> None:
    """Write a 3D-embedded RDKit Mol to XYZ format."""
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), title]
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"{atom.GetSymbol():<3}  {pos.x:12.6f}  {pos.y:12.6f}  {pos.z:12.6f}"
        )
    with open(filepath, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _write_mol_formats(mol: Chem.Mol, stem: str, formats: list[str], title: str = "") -> list[str]:
    """
    Write a 3D RDKit Mol to one or more format files.

    Parameters
    ----------
    mol     : 3D RDKit Mol (with Hs, must have a conformer)
    stem    : full path stem, e.g. '/output/aspirin' → writes aspirin.xyz etc.
    formats : list of format strings, subset of {'xyz', 'sdf', 'pdb', 'mol'}
    title   : comment line used in XYZ header

    Returns list of file paths successfully written.
    """
    written = []
    for fmt in formats:
        outfile = f"{stem}.{fmt}"
        try:
            if fmt == "xyz":
                _write_xyz(mol, outfile, title=title)
            elif fmt == "sdf":
                writer = Chem.SDWriter(outfile)
                writer.write(mol)
                writer.close()
            elif fmt == "mol":
                Chem.MolToMolFile(mol, outfile)
            elif fmt == "pdb":
                Chem.MolToPDBFile(mol, outfile)
            written.append(outfile)
            print(f"      Saved: {os.path.basename(outfile)}")
        except Exception as e:
            print(f"      Error writing {outfile}: {e}")
    return written


# ---------------------------------------------------------------------------
# 3D conformer generation (RDKit)
# ---------------------------------------------------------------------------

def embed_3d(mol: Chem.Mol) -> Chem.Mol:
    """
    Add hydrogens, embed a 3D conformer (ETKDGv3), and run MMFF94 optimisation.
    Falls back to UFF if MMFF is unavailable for the molecule.

    Returns the 3D molecule (with Hs) on success.
    Raises RuntimeError if ETKDGv3 fails to embed.

    Notes
    -----
    MMFF94 return codes: 0 = converged, 1 = not fully converged, -1 = unavailable.
    A warning is printed when the geometry did not converge (code 1).
    """
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol_h, params)
    if result == -1:
        raise RuntimeError(
            "ETKDGv3 failed to generate a 3D conformer. "
            "The molecule may be too complex or the SMILES may be invalid."
        )
    ff_result = AllChem.MMFFOptimizeMolecule(mol_h)
    if ff_result == -1:
        # MMFF not available for this molecule (e.g. some metals) — try UFF
        uff_result = AllChem.UFFOptimizeMolecule(mol_h)
        if uff_result == 1:
            print("  Warning: UFF geometry optimisation did not fully converge.")
    elif ff_result == 1:
        print("  Warning: MMFF94 geometry optimisation did not fully converge.")
    return mol_h


def _fetch_pubchem_3d_mol(cid: int) -> Chem.Mol | None:
    """
    Attempt to download a PubChem 3D SDF conformer and return it as an RDKit Mol.
    Returns None if not available or on any error.
    """
    try:
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
            "/SDF?record_type=3d"
        )
        resp = requests.get(url, timeout=15)
        if resp.status_code != 200:
            return None
        mol = Chem.MolFromMolBlock(resp.text, removeHs=False)
        return mol
    except Exception:
        return None


# ---------------------------------------------------------------------------
# get_3d_structure — main Python helper for 3D file generation
# ---------------------------------------------------------------------------

def get_3d_structure(
    query: str | int,
    formats: list[str] = None,
    output_stem: str = None,
    output_dir: str = ".",
    source: str = "auto",
    cid: int | None = None,
) -> list[str]:
    """
    Obtain a 3D structure for a compound and write it in one or more formats.

    This is the recommended single-compound entry point. Call it in a loop
    to process lists of compounds (see Examples below).

    Parameters
    ----------
    query : str or int
        PubChem CID (int or numeric string), compound name, or SMILES string.
    formats : list of str, optional
        Subset of {'xyz', 'sdf', 'pdb', 'mol'}.  Default: ['xyz'].
    output_stem : str, optional
        Full output path stem (directory + base name, no extension).
        E.g. '/results/aspirin' → writes '/results/aspirin.xyz'.
        If None, the stem is built automatically from CID and/or name
        in output_dir.
    output_dir : str, optional
        Directory for auto-named files (ignored when output_stem is given).
        Default: current directory.
    source : str, optional
        '3d_pubchem' — try PubChem 3D conformer first, fall back to RDKit.
        'rdkit'      — always use RDKit ETKDGv3 + MMFF94.
        'auto'       — same as '3d_pubchem'.
        Default: 'auto'.
    cid : int, optional
        If already known (e.g. from a prior PubChem lookup), pass it here to
        skip re-resolving.  Only used when source is 'auto' or '3d_pubchem'.

    Returns
    -------
    List of file paths that were successfully written.

    Examples
    --------
    >>> from mof_toolkit.molecule_manager import get_3d_structure

    # Single compound by CID
    >>> get_3d_structure(3033, formats=["xyz", "sdf"], output_dir="./structures/")

    # Single compound by name
    >>> get_3d_structure("aspirin", formats=["xyz"], output_dir="./structures/")

    # Single compound by SMILES
    >>> get_3d_structure("CC(=O)Oc1ccccc1C(=O)O", formats=["sdf", "pdb"],
    ...                  output_stem="./structures/aspirin")

    # Loop over a list of CIDs
    >>> cids = [3033, 3672, 644019]
    >>> for cid in cids:
    ...     get_3d_structure(cid, formats=["xyz"], output_dir="./structures/")

    # Loop over a list of SMILES with custom names
    >>> entries = [("CC(=O)O", "acetic_acid"), ("c1ccccc1", "benzene")]
    >>> for smiles, name in entries:
    ...     get_3d_structure(smiles, formats=["xyz", "sdf"],
    ...                      output_stem=f"./structures/{name}")

    # Loop over a mixed list
    >>> queries = [3033, "ibuprofen", "CC(=O)Oc1ccccc1C(=O)O"]
    >>> for q in queries:
    ...     get_3d_structure(q, formats=["xyz"], output_dir="./structures/")
    """
    if formats is None:
        formats = ["xyz"]

    # Validate formats
    valid_fmts = {"xyz", "sdf", "pdb", "mol"}
    bad = [f for f in formats if f not in valid_fmts]
    if bad:
        raise ValueError(f"Unknown format(s): {bad}. Choose from {sorted(valid_fmts)}.")

    # Resolve query to get SMILES, CID, names
    from mof_toolkit.rdkit_properties import resolve_compound_input

    try:
        resolved = resolve_compound_input(str(query))
    except Exception as e:
        print(f"  Error resolving '{query}': {e}")
        return []

    mol_2d   = resolved["mol"]
    smiles   = resolved["smiles"] or ""
    res_cid  = cid if cid is not None else resolved["cid"]
    common   = resolved["common_name"] or resolved["iupac_name"] or ""

    if mol_2d is None and not smiles:
        print(f"  Error: could not obtain a valid molecule for '{query}'.")
        return []

    # Build output stem if not provided
    if output_stem is None:
        safe_name = (
            common.replace(" ", "_").replace("/", "-").replace("\\", "-")[:40]
            if common else ""
        )
        if res_cid and safe_name:
            auto_stem = f"{res_cid}_{safe_name}"
        elif res_cid:
            auto_stem = str(res_cid)
        elif safe_name:
            auto_stem = safe_name
        else:
            from rdkit.Chem import rdMolDescriptors as _rmd
            formula = _rmd.CalcMolFormula(mol_2d) if mol_2d else "unknown"
            auto_stem = formula
        os.makedirs(output_dir, exist_ok=True)
        full_stem = os.path.join(output_dir, auto_stem)
    else:
        os.makedirs(os.path.dirname(output_stem) or ".", exist_ok=True)
        full_stem = output_stem

    # Build XYZ title string
    title_parts = []
    if res_cid:
        title_parts.append(f"CID={res_cid}")
    if common:
        title_parts.append(common)
    if smiles:
        title_parts.append(f"SMILES={smiles}")
    title = "  ".join(title_parts)

    # Try PubChem 3D conformer first (if CID is known and source allows it)
    mol_3d = None
    if source in ("auto", "3d_pubchem") and res_cid:
        print(f"  Trying PubChem 3D conformer for CID={res_cid} ...")
        mol_3d = _fetch_pubchem_3d_mol(res_cid)
        if mol_3d is not None:
            print(f"  Using PubChem 3D conformer.")
        else:
            print(f"  PubChem 3D conformer not available — falling back to RDKit ETKDGv3.")

    # Fall back to RDKit conformer generation
    if mol_3d is None:
        if mol_2d is None:
            mol_2d = Chem.MolFromSmiles(smiles)
        if mol_2d is None:
            print(f"  Error: cannot generate 3D structure — invalid molecule.")
            return []
        print(f"  Generating RDKit 3D conformer ...")
        try:
            mol_3d = embed_3d(mol_2d)
        except RuntimeError as e:
            print(f"  Error: {e}")
            return []

    return _write_mol_formats(mol_3d, full_stem, formats, title=title)


# ---------------------------------------------------------------------------
# CLI: smiles_to_3d
# ---------------------------------------------------------------------------

def smiles_to_3d_cli():
    """smiles_to_3d — convert a SMILES string to 3D structure file(s).

    3D coordinates are generated with RDKit (ETKDGv3 + MMFF94 optimisation).
    No internet connection is required.

    Arguments
    ---------
    SMILES        SMILES string (quote in the shell if it contains parentheses)
    --format, -f  Output formats: xyz sdf pdb mol  (default: xyz; multiple allowed)
    --output, -o  Filename stem, e.g. 'aspirin' → aspirin.xyz  (default: mol. formula)

    Examples
    --------
      smiles_to_3d "CC(=O)Oc1ccccc1C(=O)O"
      smiles_to_3d "CC(=O)Oc1ccccc1C(=O)O" --format xyz sdf pdb mol
      smiles_to_3d "CC(=O)Oc1ccccc1C(=O)O" --output aspirin --format sdf pdb
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="smiles_to_3d",
        description="Convert a SMILES string to 3D structure file(s) using RDKit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "3D geometry is generated locally — no internet needed.\n\n"
            "Examples:\n"
            '  smiles_to_3d "CC(=O)O"\n'
            '  smiles_to_3d "CC(=O)O" --format xyz sdf pdb mol\n'
            '  smiles_to_3d "CC(=O)O" --output acetic_acid --format sdf'
        ),
    )
    parser.add_argument(
        "smiles",
        help="SMILES string (quote it in the shell if it contains special characters)",
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
        "--output", "-o",
        default=None,
        help="Filename stem (default: molecular formula, e.g. C9H8O4)",
    )
    args = parser.parse_args()

    mol = Chem.MolFromSmiles(args.smiles)
    if mol is None:
        print(f"Error: invalid SMILES string: '{args.smiles}'")
        sys.exit(1)

    formula = rdMolDescriptors.CalcMolFormula(mol)
    stem    = args.output if args.output else formula
    title   = f"{formula}  SMILES={args.smiles}"

    print(f"Generating 3D conformer for {formula} ...")
    try:
        mol_3d = embed_3d(mol)
    except RuntimeError as e:
        print(f"Error: {e}")
        sys.exit(1)

    written = _write_mol_formats(mol_3d, stem, args.format, title=title)
    if not written:
        print("No files were written.")
        sys.exit(1)
