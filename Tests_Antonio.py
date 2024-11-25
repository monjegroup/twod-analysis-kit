import MDAnalysis as mda
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from protein2D_analysis import protein2D_analysis

trj_path='/home/antonio/Desktop/VIRMAT/Paper_PB_KDE/SIMs/RBD-PBLs_wGlyc_closed_layed/glyc_head/rep1/omicron_0/'
u=mda.Universe(f"{trj_path}md_0_1.tpr",f"{trj_path}md_short_compact.xtc")

trj_path2='/home/antonio/Desktop/VIRMAT/Paper_PB_KDE/SIMs/RBD-PBLs_wGlyc_closed_layed/glyc_head/rep1/omicron_0/'
u2=mda.Universe(f"{trj_path}md_0_1.tpr",f"{trj_path}md_short_compact.xtc")
print(u.trajectory.timeseries())
sel = u.select_atoms("resid 193-200 or protein")
ag_analysis = protein2D_analysis(sel)
print(np.shape(ag_analysis.atom_group.atoms))
ag_analysis.getPositions()
print(ag_analysis.pos.shape)



sel2 = u2.select_atoms("resid 193-200 or protein")
ag_analysis2 = protein2D_analysis(sel2)
print(np.shape(ag_analysis2.atom_group.atoms))
ag_analysis2.getPositions()
print(ag_analysis2.pos.shape)

ag_analysis_sum=ag_analysis+ag_analysis2
ag_analysis.system_name='Omicron PBL0'
print(ag_analysis_sum.universe)
########### TEST GENERAL MODULES #############
zlim=40
Nframes=200
# ag_analysis.FilterMinFrames(zlim=zlim, Nframes=Nframes, control_plots=False)
# pos=handler_from_atomgroup.getPositions(inplace=False)
# print(handler_from_atomgroup.pos.shape)
############# TEST POLAR ANALYSIS ############
# hist_arr,pos_hist=ag_analysis.PolarAnalysis('resid 193-200 or resid 12',900, sort=[1,2,3,4,5,6,7,8,0],zlim=zlim,control_plots=False,plot=True)
# print(hist_arr.shape,pos_hist.shape)
############# TEST RADII of GYRATION ANALYSIS ########

# rgs=ag_analysis.getRgs2D()
# print(rgs.shape)
# ag_analysis.RgPerpvsRgsPar(rgs, 'tab:green',show=True)

# ##########TEST Contour PLOTS ################

# paths=ag_analysis.getKDEAnalysis(zlim,Nframes)
# ag_analysis.plotPathsInLvl(1)
# areas=ag_analysis.getAreas(2,getTotal=True)
# print(areas)

##### TEST HBONDS PLOTS #####
ag_analysis.getHbonds('protein','resid 192-193', update_selections=False,trj_plot=True)



