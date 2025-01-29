import MDAnalysis as mda
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
sys.path.insert(0,'/home/antonio/Desktop/VIRMAT/2Danalysis/twod-analysis-kit/twodanalysis')
from BioPolymer2D import BioPolymer2D # type: ignore
# from twodanalysis import BioPolymer2D
from twodanalysis.data.files import MD_NOWATER_TPR, MD_TRAJ



u=mda.Universe(MD_NOWATER_TPR,MD_TRAJ)

sel = u.select_atoms("resid 193-200 or protein")
ag_analysis = BioPolymer2D(sel, surf_selection='resname DOL and name O1 and prop z > 16')
ag_analysis.system_name='Omicron PBL1'
# ag_analysis.startT=100
ag_analysis.getPositions()

zlim=15
Nframes=300

print('########### TEST GENERAL MODULES #############')

pos=ag_analysis.getPositions(inplace=False)
pos_selected=BioPolymer2D.FilterMinFrames(pos,zlim=zlim, Nframes=Nframes, control_plots=False)
# print(ag_analysis.pos.shape)
# print(pos_selected)
print('############# TEST POLAR ANALYSIS ############')
select_res='resid 198 200 12 8 40 45 111 115 173'
fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'}, figsize=(12, 12))
hist_arr,pos_hist=ag_analysis.PolarAnalysis(select_res,Nframes,# sort=[1,2,3,4,5,6,7,8,0],
                                            zlim=zlim,control_plots=False,plot=True)
# print(hist_arr.shape,pos_hist.shape)

print('############# TEST RADII OF GYRATION ANALYSIS ########')

rgs=ag_analysis.getRgs2D(plot=True)
plt.show()
# print(rgs.shape)
ag_analysis.RgPerpvsRgsPar(rgs, 'tab:green',show=True)

print('# ##########TEST CONTOUR PLOTS ################')
### CHECK
paths=ag_analysis.getKDEAnalysis(zlim,Nframes)
print(ag_analysis.kdeanalysis.kde)
fig, ax = plt.subplots(figsize=(6.5,6.5))
ag_analysis.plotPathsInLevel(paths,1,show=False,ax=ax)
areas=BioPolymer2D.getAreas(paths,2,getTotal=True)
print(areas)
ag_analysis.KDEAnalysisSelection('resid 198 200 12 8 40 45 111 115 173',Nframes,zlim,legend=True,ax=ax)
# plt.legend()
plt.show()
print('##### TEST HBONDS PLOTS #####')

sel_for_path = u.select_atoms("resid 193-200")
ag_for_path = BioPolymer2D(sel_for_path)
ag_for_path.getPositions()
paths=ag_for_path.getKDEAnalysis(zlim,Nframes)
ag_analysis.getHbonds('resname DOL','protein or resid 193-200', update_selections=False,trj_plot=False)
ag_analysis.plotHbondsPerResidues(paths,contour_lvls_to_plot=[0,5,8],top=5, print_table=True,filter=['DOL'])
plt.show()




