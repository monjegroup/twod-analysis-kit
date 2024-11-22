import MDAnalysis as mda
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from protein2D_analysis import protein2D_analysis

trj_path='/home/antonio/Desktop/VIRMAT/Paper_PB_KDE/SIMs/RBD-PBLs_wGlyc_closed_layed/glyc_head/rep1/omicron_0/'

u=mda.Universe(f"{trj_path}md_0_1.gro",f"{trj_path}md_short_compact.xtc")

sel = u.select_atoms("protein and resid 34-45")
ag_analysis = protein2D_analysis(sel)
print(np.shape(ag_analysis.atom_group.atoms))
pos=ag_analysis.getPositions()
print(ag_analysis.pos.shape)
pos=ag_analysis.getCOMs()
print(ag_analysis.com.shape)

########### TEST GENERAL MODULES #############
# zlim=40
# ag_analysis.FilterMinFrames(zlim=zlim, Nframes=200, control_plots=False)
# pos=handler_from_atomgroup.getPositions(inplace=False)
# print(handler_from_atomgroup.pos.shape)
############# TEST POLAR ANAÑYSIS ############
# hist_arr,pos_hist=ag_analysis.PolarAnalysis('resid 193-200 or resid 12',900, sort=[1,2,3,4,5,6,7,8,0],zlim=zlim,control_plots=False,plot=True)
# print(hist_arr.shape,pos_hist.shape)
############# TEST RADII of GYRATION ANALYSIS ########
# prot=ag_analysis.atom_group
# masses=prot.atoms.masses
# total_mass=np.sum(masses)

# rgs=ag_analysis.RG2D(masses, total_mass)
# print(rgs)
# ag_analysis.getRgs2D()
rgs=ag_analysis.RGs2D()
print(rgs)

