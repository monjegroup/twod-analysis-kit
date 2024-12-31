
import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from twodanalysis import BioPolymer2D
from twodanalysis.data.files import MD_NOWATER_TPR, MD_TRAJ


print('########### TEST GENERAL METHODS #############')

def test_getPositions(obj,pos_type='COM', inplace=True, select=None,getselection=False):
    pos=obj.getPositions(pos_type, inplace, select,getselection)
    if inplace:
        assert pos is None and isinstance(obj.pos, (list, np.ndarray))
    elif getselection is True:
        assert isinstance(pos, (list, np.ndarray)), isinstance(obj.atom_group, MDAnalysis.AtomGroup)  # type: ignore
    else:
        assert isinstance(pos, (list, np.ndarray))

def test_getCOMs(obj, inplace=True, select=None):
    com=obj.getCOMs(inplace, select)
    if inplace:
        assert com is None and isinstance(obj.com, (list, np.ndarray))
    else:
        assert isinstance(com, (list, np.ndarray))
        
def test_FilterMinFrames(obj, pos, zlim,Nframes,control_plots=False):
    selected_pos=obj.FilterMinFrames(pos, zlim,Nframes,control_plots)
    assert isinstance(selected_pos, (list, np.ndarray))

def test_PolarAnalysis(obj,select_res,Nframes,max_hist=None,sort=None,plot=False,control_plots=False, zlim=14,Nbins=1000,resolution=5):
    hist,ordered_selected_pos=obj.PolarAnalysis(select_res,Nframes,max_hist,sort,plot,control_plots, zlim,Nbins,resolution)
    assert isinstance(hist, (list, np.ndarray)) and np.shape(hist)==(np.shape(ordered_selected_pos)[1],2,Nbins-1) and isinstance(ordered_selected_pos, (list, np.ndarray))

def test_computeRG2D(obj):
    masses=obj.atom_group.atoms.masses
    total_mass=np.sum(masses)
    rg=obj.computeRG2D(masses,total_mass)
    print(np.shape(rg),'computeRG2D')
    assert isinstance(rg, (list, np.ndarray)) and np.shape(rg)[0]==3

def test_getRgs2D(obj):
    rgs=obj.getRgs2D()
    assert np.shape(rgs)==(obj.universe.trajectory[-1].frame,4)
def test_RgPerpvsRgsPar(obj,rgs,color, marker='s',plot=False,show=False):
    rg_ratio=obj.RgPerpvsRgsPar(rgs,color, marker,plot,show)
    assert isinstance(rg_ratio,(float))




u=mda.Universe(MD_NOWATER_TPR,MD_TRAJ) 

sel = u.select_atoms("resid 193-200 or protein")
ag = BioPolymer2D(sel)

test_getPositions(ag,pos_type='COM', inplace=True, select=None,getselection=False)
test_getCOMs(ag, inplace=True, select=None)
zlim,Nframes=60,300
pos=ag.getPositions(inplace=False)
test_FilterMinFrames(ag, pos, zlim,Nframes)
select_res="resid 111 200 100"
print('############# TEST POLAR ANALYSIS ############')
test_PolarAnalysis(ag,select_res,Nframes,zlim=zlim)

print('############# TEST RADII OF GYRATION ANALYSIS ########')

test_computeRG2D(ag)
test_getRgs2D(ag)

rgs=ag.getRgs2D(plot=False)
# plt.show()
# # print(rgs.shape)
test_RgPerpvsRgsPar(ag,rgs, color='green')


# print('# ##########TEST CONTOUR PLOTS ################')

# paths=ag_analysis.getKDEAnalysis(zlim,Nframes)
# print(ag_analysis.kdeanalysis.kde)
# ag_analysis.plotPathsInLevel(paths,1,show=False)
# areas=BioPolymer2D.getAreas(paths,2,getTotal=True)
# print(areas)
# ag_analysis.KDEAnalysisSelection('resid 198 200 12 8 40 45 111 115 173',Nframes,zlim)
# # plt.legend()
# plt.show()
# print('##### TEST HBONDS PLOTS #####')

# sel_for_path = u.select_atoms("resid 193-200")
# ag_for_path = BioPolymer2D(sel_for_path)
# ag_for_path.getPositions()
# paths=ag_for_path.getKDEAnalysis(zlim,Nframes)
# ag_analysis.getHbonds('resname DOL','resid 193-200', update_selections=False,trj_plot=False)
# ag_analysis.plotHbondsPerResidues(paths,contour_lvls_to_plot=[0,5,8],top=5, print_table=True,filter=['DOL'])



