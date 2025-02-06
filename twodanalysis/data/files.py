"""
Location of data files
======================

Use as ::

    from twodanalysis.data.files import *

"""

__all__ = [
    "MDANALYSIS_LOGO",  # example file of MDAnalysis logo
]

import importlib.resources

data_directory = importlib.resources.files("twodanalysis") / "data"

MDANALYSIS_LOGO = data_directory / "mda.txt"
MEMBRANE_XTC = data_directory / "membrane_new.xtc"
MEMBRANE_TPR = data_directory / "veamos.tpr"
MD_NOWATER_TPR = data_directory / "md_nowater.tpr"
MD_NOWATER_PDB = data_directory / "md_nowater.pdb"
MD_TRAJ = data_directory / "trajout.xtc"
