"""
molecule_manager.py
-------------------
3D conformer generation, structure file writing, and format conversion for
the MOF-Guest Toolkit.

Python helpers
--------------
    embed_3d(mol)
        Add Hs, embed a 3D conformer (ETKDGv3), optimise with MMFF94 (or UFF).
        Returns the 3D Mol with Hs on success; raises RuntimeError on failure.

    get_smiles_from_coords(filepath)
        Read a 3D structure file (xyz, sdf, mol, pdb) and return the canonical
        SMILES of the molecule it contains.

    get_3d_structure(query, formats, output_stem, output_dir, source, cid)
        Obtain a 3D structure for a CID, name, formula, SMILES, or InChIKey
        and write it in one or more formats.  Single-compound helper — call
        in a loop for lists.

CLI commands
------------
    mol_get_xyz -input <cid|name|smiles|…> [-outputformat <xyz sdf pdb mol>]
                [-output <stem>]
        Convert any supported identifier to 3D structure file(s).

    mol_get_xyz -batch <csv> [-outputformat <…>] [-output <dir>]
        Batch-generate 3D structure files from a CSV file.

    mol_file_translate -input <file> -output <file>
        Convert between any two supported 3D structure formats.
"""

from __future__ import annotations

import csv
import os
import re
import sys

import requests
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


# Supported 3D format extensions
SUPPORTED_FORMATS: frozenset[str] = frozenset({"xyz", "sdf", "mol", "pdb"})


# ---------------------------------------------------------------------------
# Low-level format writers
# ---------------------------------------------------------------------------

def _write_xyz(mol: Chem.Mol, filepath: str, title: str = "") -> None:
    """
    Write a 3D-embedded RDKit Mol to XYZ format.

    The first line is the atom count, the second is the title/comment,
    and subsequent lines are: symbol  x  y  z.
    """
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), title]
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"{atom.GetSymbol():<3}  {pos.x:12.6f}  {pos.y:12.6f}  {pos.z:12.6f}"
        )
    with open(filepath, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def _write_mol_formats(
    mol: Chem.Mol,
    stem: str,
    formats: list[str],
    title: str = "",
) -> list[str]:
    """
    Write a 3D RDKit Mol to one or more format files.

    Parameters
    ----------
    mol     : 3D RDKit Mol with Hs and a conformer.
    stem    : Output path stem (no extension), e.g. '/output/aspirin'.
    formats : List of format strings, subset of SUPPORTED_FORMATS.
    title   : Comment line used in XYZ header and SDF title field.

    Returns
    -------
    List of file paths that were successfully written.
    """
    written: list[str] = []
    for fmt in formats:
        outfile = f"{stem}.{fmt}"
        try:
            if fmt == "xyz":
                _write_xyz(mol, outfile, title=title)
            elif fmt == "sdf":
                writer = Chem.SDWriter(outfile)
                if title:
                    mol.SetProp("_Name", title)
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
# Format readers — used by get_smiles_from_coords and mol_file_translate
# ---------------------------------------------------------------------------

_FORMAT_READERS: dict[str, object] = {
    "sdf": Chem.MolFromMolFile,
    "mol": Chem.MolFromMolFile,
    "pdb": Chem.MolFromPDBFile,
    "xyz": None,  # handled separately below
}


def _read_mol_from_file(filepath: str) -> Chem.Mol | None:
    """
    Read a molecule from a structure file (sdf, mol, pdb, or xyz).

    For XYZ files, RDKit has no native reader; we use Open Babel if available
    (via subprocess) and fall back to a pure-Python heavy-atom reader that
    omits bond orders (useful for SMILES generation via distance geometry).

    Returns an RDKit Mol on success, None on failure.
    """
    ext = os.path.splitext(filepath)[1].lstrip(".").lower()
    if ext not in SUPPORTED_FORMATS:
        print(f"  Error: unsupported format '.{ext}'. "
              f"Supported: {', '.join(sorted(SUPPORTED_FORMATS))}")
        return None

    # ---- SDF / MOL
    if ext in ("sdf", "mol"):
        mol = Chem.MolFromMolFile(filepath, removeHs=False, sanitize=True)
        return mol

    # ---- PDB
    if ext == "pdb":
        mol = Chem.MolFromPDBFile(filepath, removeHs=False, sanitize=True)
        return mol

    # ---- XYZ — try Open Babel first, then naive fallback
    if ext == "xyz":
        # Attempt via Open Babel (openbabel Python bindings)
        try:
            from openbabel import pybel  # type: ignore
            ob_mol = next(pybel.readfile("xyz", filepath))
            sdf_str = ob_mol.write("sdf")
            mol = Chem.MolFromMolBlock(sdf_str, removeHs=False, sanitize=True)
            if mol is not None:
                return mol
        except Exception:
            pass

        # Naive fallback: parse atom symbols + coordinates, build mol by
        # RDKit's AllChem.MolFromSmiles(AtomBlock) — returns heavy atoms only,
        # useful for re-generating SMILES when bond-order info is not needed.
        try:
            with open(filepath, "r", encoding="utf-8") as fh:
                lines = [l.strip() for l in fh if l.strip()]
            # Skip count line and title line
            start = 2
            emol = Chem.RWMol()
            conf_atoms: list[tuple[float, float, float]] = []
            for line in lines[start:]:
                parts = line.split()
                if len(parts) < 4:
                    continue
                sym, x, y, z = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
                try:
                    atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(sym)
                    emol.AddAtom(Chem.Atom(atomic_num))
                    conf_atoms.append((x, y, z))
                except Exception:
                    pass
            if conf_atoms:
                from rdkit.Geometry import rdGeometry
                conf = Chem.Conformer(len(conf_atoms))
                for idx, (x, y, z) in enumerate(conf_atoms):
                    conf.SetAtomPosition(idx, rdGeometry.Point3D(x, y, z))
                emol.AddConformer(conf, assignId=True)
                try:
                    Chem.SanitizeMol(emol)
                    return emol.GetMol()
                except Exception:
                    pass
        except Exception as e:
            print(f"  Error reading XYZ file '{filepath}': {e}")
        return None

    return None


# ---------------------------------------------------------------------------
# 3D conformer generation (RDKit)
# ---------------------------------------------------------------------------

def embed_3d(mol: Chem.Mol) -> Chem.Mol:
    """
    Add hydrogens, embed a 3D conformer (ETKDGv3), and optimise with MMFF94.

    Falls back to UFF optimisation when MMFF94 is unavailable for the molecule
    (e.g. organometallics or unusual atom types).

    Parameters
    ----------
    mol : Chem.Mol
        2D or 3D RDKit molecule (Hs are added internally).

    Returns
    -------
    Chem.Mol
        New Mol object with explicit Hs and a 3D conformer.

    Raises
    ------
    RuntimeError
        If ETKDGv3 fails to generate a conformer.  This typically indicates
        that the molecule is too complex, highly strained, or the SMILES is
        malformed.

    Notes
    -----
    MMFF94 return codes: 0 = converged, 1 = not fully converged, -1 = unavailable.
    A warning is printed when convergence is incomplete (code 1).
    """
    mol_h  = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42

    result = AllChem.EmbedMolecule(mol_h, params)
    if result == -1:
        raise RuntimeError(
            "ETKDGv3 failed to generate a 3D conformer.  "
            "The molecule may be too complex, highly strained, or the SMILES invalid."
        )

    ff_code = AllChem.MMFFOptimizeMolecule(mol_h)
    if ff_code == -1:
        # MMFF94 not available — try UFF
        uff_code = AllChem.UFFOptimizeMolecule(mol_h)
        if uff_code == 1:
            print("  Warning: UFF optimisation did not fully converge.")
    elif ff_code == 1:
        print("  Warning: MMFF94 optimisation did not fully converge.")

    return mol_h


def _fetch_pubchem_3d_sdf(cid: int) -> str | None:
    """
    Download a PubChem 3D SDF conformer for *cid*.
    Returns the SDF string on success, None if unavailable or on error.
    """
    try:
        resp = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
            "/SDF?record_type=3d",
            timeout=15,
        )
        return resp.text if resp.status_code == 200 else None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# get_smiles_from_coords — read SMILES from a 3D file
# ---------------------------------------------------------------------------

def get_smiles_from_coords(filepath: str) -> str | None:
    """
    Read a 3D structure file and return its canonical SMILES string.

    Supported formats: xyz, sdf, mol, pdb.

    For XYZ files, Open Babel is used when available (provides bond order
    perception); otherwise a distance-geometry fallback is used which may
    give approximate results for complex molecules.

    Parameters
    ----------
    filepath : str
        Path to the 3D structure file.

    Returns
    -------
    Canonical SMILES string, or None if the file cannot be read or converted.

    Examples
    --------
    >>> get_smiles_from_coords("aspirin.sdf")
    'CC(=O)Oc1ccccc1C(=O)O'
    >>> get_smiles_from_coords("molecule.xyz")
    'CC(=O)O'
    """
    mol = _read_mol_from_file(filepath)
    if mol is None:
        return None
    mol_no_h = Chem.RemoveHs(mol)
    return Chem.MolToSmiles(mol_no_h) if mol_no_h else None


# ---------------------------------------------------------------------------
# get_3d_structure — main Python helper
# ---------------------------------------------------------------------------

def get_3d_structure(
    query: str | int,
    formats: list[str] | None = None,
    output_stem: str | None = None,
    output_dir: str = ".",
    source: str = "auto",
    cid: int | None = None,
) -> list[str]:
    """
    Obtain a 3D structure for a compound and write it in one or more formats.

    Parameters
    ----------
    query : str or int
        CID (int or numeric string), compound name, formula, SMILES, or InChIKey.
    formats : list of str, optional
        Subset of {'xyz', 'sdf', 'pdb', 'mol'}.  Default: ['xyz'].
    output_stem : str, optional
        Full output path stem (directory + base name, no extension).
        E.g. '/results/aspirin' → writes '/results/aspirin.xyz'.
        When None, the stem is built automatically from CID and/or name
        inside *output_dir*.
    output_dir : str, optional
        Directory for auto-named files (ignored when *output_stem* is given).
        Created automatically if it does not exist.  Default: current directory.
    source : str, optional
        '3d_pubchem' — try PubChem 3D conformer first, fall back to RDKit.
        'rdkit'      — always use RDKit ETKDGv3 + MMFF94 locally.
        'auto'       — same as '3d_pubchem'.
        Default: 'auto'.
    cid : int, optional
        Already-known PubChem CID; skips re-resolution when provided.

    Returns
    -------
    List of file paths that were successfully written.

    Examples
    --------
    >>> from mof_toolkit.molecule_manager import get_3d_structure

    # By CID — downloads PubChem 3D conformer when available
    >>> get_3d_structure(3033, formats=["xyz", "sdf"], output_dir="./structures/")

    # By name
    >>> get_3d_structure("aspirin", formats=["xyz"], output_dir="./structures/")

    # By SMILES — always uses RDKit (no CID known)
    >>> get_3d_structure("CC(=O)Oc1ccccc1C(=O)O", formats=["sdf", "pdb"],
    ...                  output_stem="./structures/aspirin")

    # Force RDKit even when a CID is available
    >>> get_3d_structure(3033, formats=["xyz"], output_dir="./out/", source="rdkit")

    # Loop over a list of CIDs
    >>> for c in [3033, 3672, 644019]:
    ...     get_3d_structure(c, formats=["xyz"], output_dir="./structures/")

    # Loop with custom output names
    >>> for smiles, name in [("CC(=O)O", "acetic_acid"), ("c1ccccc1", "benzene")]:
    ...     get_3d_structure(smiles, formats=["xyz", "sdf"],
    ...                      output_stem=f"./structures/{name}")
    """
    from mof_toolkit.rdkit_descriptors import resolve_compound_input

    if formats is None:
        formats = ["xyz"]

    # Validate requested formats
    bad = [f for f in formats if f not in SUPPORTED_FORMATS]
    if bad:
        raise ValueError(
            f"Unknown format(s): {bad}. "
            f"Choose from: {sorted(SUPPORTED_FORMATS)}"
        )

    # Resolve the query to a unified metadata dict
    try:
        resolved = resolve_compound_input(str(query))
    except Exception as e:
        print(f"  Error resolving '{query}': {e}")
        return []

    mol_2d       = resolved["mol"]
    smiles       = resolved["smiles"] or ""
    res_cid      = cid if cid is not None else resolved["cid"]
    common_name  = resolved["common_name"] or resolved["iupac_name"] or ""
    formula      = resolved["formula"]

    if mol_2d is None and not smiles:
        print(f"  Error: could not obtain a valid molecule for '{query}'.")
        return []

    # Build the output path stem when not explicitly provided
    if output_stem is None:
        safe_name = (
            re.sub(r"[/\\]", "-", common_name).replace(" ", "_")[:40]
            if common_name else ""
        )
        if res_cid and safe_name:
            auto_stem = f"{res_cid}_{safe_name}"
        elif res_cid:
            auto_stem = str(res_cid)
        elif safe_name:
            auto_stem = safe_name
        else:
            auto_stem = formula or "unknown"
        os.makedirs(output_dir, exist_ok=True)
        full_stem = os.path.join(output_dir, auto_stem)
    else:
        parent = os.path.dirname(output_stem) or "."
        os.makedirs(parent, exist_ok=True)
        full_stem = output_stem

    # Build the XYZ/SDF title string
    title_parts: list[str] = []
    if res_cid:
        title_parts.append(f"CID={res_cid}")
    if common_name:
        title_parts.append(common_name)
    if smiles:
        title_parts.append(f"SMILES={smiles}")
    title = "  ".join(title_parts)

    # Try PubChem 3D conformer first when the CID is known
    mol_3d: Chem.Mol | None = None
    if source in ("auto", "3d_pubchem") and res_cid:
        print(f"  Trying PubChem 3D conformer for CID={res_cid} …")
        sdf_text = _fetch_pubchem_3d_sdf(res_cid)
        if sdf_text:
            mol_3d = Chem.MolFromMolBlock(sdf_text, removeHs=False)
            if mol_3d is not None:
                print("  Using PubChem 3D conformer.")
            else:
                print("  PubChem returned SDF but RDKit could not parse it — "
                      "falling back to RDKit ETKDGv3.")
        else:
            print("  PubChem 3D conformer not available — "
                  "falling back to RDKit ETKDGv3.")

    # Fall back to local RDKit conformer generation
    if mol_3d is None:
        if mol_2d is None:
            mol_2d = Chem.MolFromSmiles(smiles)
        if mol_2d is None:
            print(f"  Error: cannot generate 3D structure — invalid molecule.")
            return []
        print("  Generating RDKit 3D conformer (ETKDGv3 + MMFF94) …")
        try:
            mol_3d = embed_3d(mol_2d)
        except RuntimeError as e:
            print(f"  Error: {e}")
            return []

    return _write_mol_formats(mol_3d, full_stem, formats, title=title)


# ---------------------------------------------------------------------------
# mol_file_translate — format conversion between structure files
# ---------------------------------------------------------------------------

def mol_file_translate_cli() -> None:
    """
    mol_file_translate — convert a 3D structure file between supported formats.

    Usage
    -----
      mol_file_translate -input molecule.sdf -output molecule.xyz
      mol_file_translate -input molecule.xyz -output molecule.pdb
      mol_file_translate -input molecule.pdb -output molecule.mol

    Supported formats: xyz, sdf, mol, pdb

    The input format is inferred from the input file extension.
    The output format is inferred from the output file extension.
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="mol_file_translate",
        description="Convert a 3D structure file between xyz/sdf/mol/pdb formats",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  mol_file_translate -input aspirin.sdf -output aspirin.xyz\n"
            "  mol_file_translate -input molecule.xyz -output molecule.pdb\n"
            "  mol_file_translate -input compound.pdb -output compound.mol\n"
        ),
    )
    parser.add_argument(
        "-input", dest="input", required=True, metavar="FILE",
        help="Input structure file (xyz, sdf, mol, pdb)",
    )
    parser.add_argument(
        "-output", dest="output", required=True, metavar="FILE",
        help="Output structure file; format inferred from extension",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: input file not found: {args.input}")
        sys.exit(1)

    in_ext  = os.path.splitext(args.input)[1].lstrip(".").lower()
    out_ext = os.path.splitext(args.output)[1].lstrip(".").lower()

    if in_ext not in SUPPORTED_FORMATS:
        print(f"Error: unsupported input format '.{in_ext}'. "
              f"Supported: {', '.join(sorted(SUPPORTED_FORMATS))}")
        sys.exit(1)
    if out_ext not in SUPPORTED_FORMATS:
        print(f"Error: unsupported output format '.{out_ext}'. "
              f"Supported: {', '.join(sorted(SUPPORTED_FORMATS))}")
        sys.exit(1)

    print(f"Reading {args.input} …")
    mol = _read_mol_from_file(args.input)
    if mol is None:
        print("Error: could not read the input file.")
        sys.exit(1)

    # Use the input filename (without extension) as the title
    title = os.path.splitext(os.path.basename(args.input))[0]
    stem  = os.path.splitext(args.output)[0]

    written = _write_mol_formats(mol, stem, [out_ext], title=title)
    if not written:
        print("Error: output file could not be written.")
        sys.exit(1)
    print(f"Converted: {args.input}  →  {args.output}")


# ---------------------------------------------------------------------------
# CLI: mol_get_xyz
# ---------------------------------------------------------------------------

def mol_get_xyz_cli() -> None:
    """
    mol_get_xyz — generate 3D structure file(s) from any compound identifier.

    Single compound
    ---------------
      mol_get_xyz -input <cid|name|smiles|formula|inchikey>
          Generate an XYZ file (default format).

      mol_get_xyz -input aspirin -outputformat xyz sdf pdb mol
          Generate files in multiple formats.

      mol_get_xyz -input aspirin -outputformat sdf -output ./structures/asp
          Write to a custom output stem (./structures/asp.sdf).

    Batch from CSV
    --------------
      mol_get_xyz -batch molecules.csv -output ./structures/
          Generate XYZ files for every compound in the CSV.

      mol_get_xyz -batch molecules.csv -outputformat xyz sdf -output ./out/
          Generate multiple format files per compound.

    The CSV must have at least one of: CID, Name, SMILES, InChIKey.
    3D structures are fetched from PubChem (when a CID is available) or
    generated locally with RDKit ETKDGv3 + MMFF94.
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="mol_get_xyz",
        description="Generate 3D structure file(s) from any compound identifier",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  mol_get_xyz -input 3033\n"
            "  mol_get_xyz -input aspirin -outputformat xyz sdf\n"
            "  mol_get_xyz -input \"CC(=O)O\" -outputformat sdf -output acetic_acid\n"
            "  mol_get_xyz -batch mols.csv -outputformat xyz -output ./structures/\n"
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
        "-outputformat", dest="outputformat", nargs="+",
        choices=["xyz", "sdf", "pdb", "mol"],
        default=["xyz"],
        metavar="FORMAT",
        help="Output format(s): xyz sdf pdb mol  (default: xyz)",
    )
    parser.add_argument(
        "-output", dest="output", metavar="PATH", default=None,
        help=(
            "Single-compound: output stem (e.g. 'aspirin' → aspirin.xyz). "
            "Batch: output directory (default: current directory)."
        ),
    )
    parser.add_argument(
        "--source", choices=["auto", "3d_pubchem", "rdkit"], default="auto",
        help="3D source: auto (PubChem first), 3d_pubchem, or rdkit (default: auto)",
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Single-compound mode
    # ------------------------------------------------------------------
    if args.input:
        output_stem = args.output or None
        output_dir  = "."
        if output_stem and os.path.isdir(output_stem):
            output_dir  = output_stem
            output_stem = None

        written = get_3d_structure(
            query=args.input,
            formats=args.outputformat,
            output_stem=output_stem,
            output_dir=output_dir,
            source=args.source,
        )
        if not written:
            print("No files were written.")
            sys.exit(1)
        return

    # ------------------------------------------------------------------
    # Batch mode
    # ------------------------------------------------------------------
    if not os.path.isfile(args.batch):
        print(f"Error: input file not found: {args.batch}")
        sys.exit(1)

    output_dir = args.output or "."
    os.makedirs(output_dir, exist_ok=True)

    with open(args.batch, newline="", encoding="utf-8") as fh:
        reader  = csv.DictReader(fh)
        headers = list(reader.fieldnames or [])
        rows    = list(reader)

    from mof_toolkit.rdkit_descriptors import _resolve_batch_row, _BATCH_COL_ALIASES

    accepted_cols = {alias for aliases in _BATCH_COL_ALIASES.values() for alias in aliases}
    if not accepted_cols.intersection(headers):
        print(
            "Error: input CSV must contain at least one of: "
            "CID, Name, SMILES, InChIKey"
        )
        sys.exit(1)

    total           = len(rows)
    missing_counter = 0

    for i, row in enumerate(rows, 1):
        res = _resolve_batch_row(row, i, total)
        if res is None:
            continue

        cid         = res["cid"]
        common_name = res["common_name"] or res["iupac_name"]
        safe_name   = (
            re.sub(r"[/\\]", "-", common_name).replace(" ", "_")[:40]
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

        print(f"    Writing {', '.join(args.outputformat)} → {stem}.*")
        get_3d_structure(
            query=res["smiles"] or (str(cid) if cid else ""),
            formats=args.outputformat,
            output_stem=os.path.join(output_dir, stem),
            source=args.source,
            cid=cid,
        )