import MDAnalysis as mda
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from protein2D_analysis import protein2D_analysis

trj_path='/home/antonio/Desktop/VIRMAT/Paper_PB_KDE/SIMs/RBD-PBLs_wGlyc_closed_layed/glyc_head/rep1/omicron_10/'
u=mda.Universe(f"{trj_path}md_0_1.tpr",f"{trj_path}md_short_compact.xtc")
sel = u.select_atoms("resid 193-200 or protein")
ag_analysis = protein2D_analysis(sel)
ag_analysis.getPositions()
print(ag_analysis.pos.shape)


########### TEST GENERAL MODULES #############
zlim=60
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

paths=ag_analysis.getKDEAnalysis(zlim,Nframes)
# ag_analysis.plotPathsInLvl(1)
# areas=ag_analysis.getAreas(2,getTotal=True)
# print(areas)

##### TEST HBONDS PLOTS #####
ag_analysis.getHbonds('protein','resid 193-200', update_selections=False,trj_plot=False)
df=ag_analysis.HbondsPerResidues()
print(df)

max_val=df['Count'].max()
str_resids=' '.join(np.array(df['ResIDs'],dtype=str))
print(str_resids)
res_w_hbonds=u.select_atoms(f'resid {str_resids}')
ag_w_hbonds=protein2D_analysis(res_w_hbonds)
ag_w_hbonds.getPositions()
ag_analysis.plotPathsInLvl(0)
ag_analysis.plotPathsInLvl(5)
ag_analysis.plotPathsInLvl(8)
for i in range(ag_w_hbonds.pos.shape[1]):
    norm_val=df['Count'].iloc[i]/len(ag_w_hbonds.universe.trajectory) #max_val
    norm_val_plot=df['Count'].iloc[i]/max_val
    pos=ag_w_hbonds.pos[:,i].mean(axis=0)
    plt.plot(pos[1],pos[3], 'o',
            label='%s-%s (%.2f)'%(ag_w_hbonds.atom_group.residues.resids[i],
                                  ag_w_hbonds.atom_group.residues.resnames[i], norm_val*100),)
    plt.scatter(pos[1],pos[3], s=(8*20*norm_val_plot)**2, alpha=.5)
plt.legend()
plt.show()



