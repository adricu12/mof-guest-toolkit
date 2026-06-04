"""
molecule_manager.py
-------------------
Local molecule manipulation — no PubChem dependency.

CLI commands
------------
  smiles_to_3d <SMILES> [--format xyz sdf pdb mol] [--output stem]
    Convert a SMILES string to 3D coordinates using RDKit's ETKDG conformer
    generator and MMFF94 geometry optimisation.  Writes one file per format
    requested.  Supported formats: xyz, sdf, pdb, mol.
"""

import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


# ---------------------------------------------------------------------------
# Helpers
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


def embed_3d(mol: Chem.Mol) -> Chem.Mol:
    """
    Add hydrogens, embed a 3D conformer (ETKDGv3), and run MMFF94 optimisation.
    Falls back to UFF if MMFF is unavailable for the molecule.
    Returns the 3D molecule (with Hs) on success, raises RuntimeError on failure.
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
        AllChem.UFFOptimizeMolecule(mol_h)
    return mol_h


# ---------------------------------------------------------------------------
# CLI: smiles_to_3d
# ---------------------------------------------------------------------------

def smiles_to_3d_cli():
    """smiles_to_3d <SMILES> [--format xyz sdf pdb mol] [--output stem]

    Convert a SMILES string to one or more 3D structure files.
    3D coordinates are generated with RDKit (ETKDGv3 + MMFF94).

    Arguments
    ---------
    SMILES      SMILES string (use quotes in the shell if it contains parentheses)
    --format    One or more of: xyz sdf pdb mol  (default: xyz)
    --output    Filename stem, e.g. 'aspirin' → aspirin.xyz  (default: mol. formula)

    Examples
    --------
      smiles_to_3d "CC(=O)Oc1ccccc1C(=O)O"
      smiles_to_3d "CC(=O)Oc1ccccc1C(=O)O" --format xyz sdf pdb mol
      smiles_to_3d "CC(=O)Oc1ccccc1C(=O)O" --output aspirin --format sdf
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert SMILES to 3D structure file(s) using RDKit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
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

    # Parse SMILES
    mol = Chem.MolFromSmiles(args.smiles)
    if mol is None:
        print(f"Error: invalid SMILES string: '{args.smiles}'")
        sys.exit(1)

    # Determine output stem
    formula = rdMolDescriptors.CalcMolFormula(mol)
    stem = args.output if args.output else formula

    # Generate 3D conformer
    print(f"Generating 3D conformer for {formula} ...")
    try:
        mol_3d = embed_3d(mol)
    except RuntimeError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Write requested formats
    written = []
    for fmt in args.format:
        outfile = f"{stem}.{fmt}"
        try:
            if fmt == "sdf":
                writer = Chem.SDWriter(outfile)
                writer.write(mol_3d)
                writer.close()
            elif fmt == "mol":
                Chem.MolToMolFile(mol_3d, outfile)
            elif fmt == "pdb":
                Chem.MolToPDBFile(mol_3d, outfile)
            elif fmt == "xyz":
                _write_xyz(mol_3d, outfile, title=f"{formula}  SMILES={args.smiles}")
            written.append(outfile)
            print(f"  Saved: {outfile}")
        except Exception as e:
            print(f"  Error writing {outfile}: {e}")

    if not written:
        print("No files were written.")
        sys.exit(1)
