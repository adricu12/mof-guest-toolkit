"""
mof_toolkit
-----------
MOF-Guest Toolkit — utilities for computational chemistry workflows with
an emphasis on MOF–guest systems, cheminformatics, and research automation.

Quick-start
-----------
>>> from mof_toolkit.rdkit_descriptors import (
...     resolve_compound_input,
...     get_rdkit_dict,
...     get_descriptor,
...     get_smiles,
...     display_table,
... )
>>> from mof_toolkit.molecule_manager import (
...     get_3d_structure,
...     get_smiles_from_coords,
...     embed_3d,
... )

See the individual module docstrings for full API documentation.
"""

from mof_toolkit.rdkit_descriptors import (
    resolve_compound_input,
    get_rdkit_dict,
    get_descriptor,
    get_smiles,
    get_all_rdkit_descriptors,
    display_table,
    DEFAULT_DESCRIPTORS,
    VIEWER_DESCRIPTORS,
)

from mof_toolkit.molecule_manager import (
    get_3d_structure,
    get_smiles_from_coords,
    embed_3d,
)

__all__ = [
    "resolve_compound_input",
    "get_rdkit_dict",
    "get_descriptor",
    "get_all_rdkit_descriptors",
    "display_table",
    "DEFAULT_DESCRIPTORS",
    "VIEWER_DESCRIPTORS",
    "get_smiles",
    "get_3d_structure",
    "get_smiles_from_coords",
    "embed_3d",
]