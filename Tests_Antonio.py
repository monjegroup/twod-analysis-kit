import MDAnalysis as mda
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from protein2D_analysis import protein2D_analysis

trj_path='/home/antonio/Desktop/VIRMAT/Paper_PB_KDE/SIMs/RBD-PBLs_wGlyc_closed_layed/glyc_head/rep1/omicron_0/'

u=mda.Universe(f"{trj_path}md_0_1.gro",f"{trj_path}md_short_compact.xtc")

sel = u.select_atoms("protein")
handler_from_atomgroup = protein2D_analysis(sel,)
print(handler_from_atomgroup.atom_group.atoms)
pos=handler_from_atomgroup.getPositions()
print(pos.shape)