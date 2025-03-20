
import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from twodanalysis import BioPolymer2D
from twodanalysis.data.files import MD_NOWATER_TPR, MD_TRAJ
import pytest

surf_selection='resname DOL and name O1 and prop z > 20' ## This is dependent on the system used for testing

def test_import():
    assert BioPolymer2D is not None

def test_init():
    u=mda.Universe(MD_NOWATER_TPR, MD_TRAJ)
    obj = BioPolymer2D(u)
    assert obj is not None
    assert isinstance(obj.universe, mda.Universe)
    assert isinstance(obj.atom_group, mda.AtomGroup)
    

@pytest.mark.parametrize("pos_type", ['COM', 'all'])
@pytest.mark.parametrize("inplace", [True, False])
@pytest.mark.parametrize("select", [None, 'name CA'])###
@pytest.mark.parametrize("getselection", [False, True])
def test_get_positions(pos_type, inplace, select,getselection):
    u=mda.Universe(MD_NOWATER_TPR, MD_TRAJ)
    obj = BioPolymer2D(u)
    pos=obj.getPositions(pos_type, inplace=inplace, select=select, getselection=getselection)
    if inplace:
        assert pos is None and isinstance(obj.pos, (list, np.ndarray))
    elif getselection is True:
        pos, ag = pos
        assert isinstance(pos, (list, np.ndarray))
        assert isinstance(ag, mda.AtomGroup)
    else:
        assert isinstance(pos, (list, np.ndarray))



@pytest.mark.parametrize("inplace", [True, False])
@pytest.mark.parametrize("select", [None, 'name CA'])
def test_get_COMs(inplace, select):
    u=mda.Universe(MD_NOWATER_TPR, MD_TRAJ)
    obj = BioPolymer2D(u)
    com=obj.getCOMs(inplace, select)
    if inplace:
        assert com is None and isinstance(obj.com, (list, np.ndarray))
    else:
        assert isinstance(com, (list, np.ndarray))

@pytest.mark.parametrize("pos_type", ['COM', 'all'])
@pytest.mark.parametrize("inplace", [True, False])
@pytest.mark.parametrize("select", [None, 'name CA'])
@pytest.mark.parametrize("Nframes", [None, 10])
@pytest.mark.parametrize("zlim", [14]) # Add test for 2 and 14, to validated interanl error management of the function
def test_FilterMinFrames(pos_type, inplace, select, Nframes,zlim,surf_selection=surf_selection):
    u=mda.Universe(MD_NOWATER_TPR, MD_TRAJ,)
    obj = BioPolymer2D(u,surf_selection=surf_selection) 
    if inplace:
        obj.getPositions(pos_type, inplace=inplace, select=select)
        pos=obj.pos
    else:
        pos=obj.getPositions(pos_type, inplace=inplace, select=select)
    # if Nframes is different from None, the function should pass if raise a ValueError or not
    # if not Nframes is None:
    #     with pytest.raises(ValueError):
    #         selected_pos=obj.FilterMinFrames(pos, zlim=zlim, Nframes=Nframes, control_plots=False)
    # else:
    selected_pos=obj.FilterMinFrames(pos, zlim=zlim, Nframes=Nframes, control_plots=False)
    assert isinstance(selected_pos, (list, np.ndarray))
def test_PolarAnalysis(surf_selection=surf_selection):
    u=mda.Universe(MD_NOWATER_TPR, MD_TRAJ)
    obj = BioPolymer2D(u,surf_selection=surf_selection)
    select_res='resid 198 200 12 8 40 45 111 115 173'
    Nframes=100
    zlim=14
    nbins=1000
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'}, figsize=(12, 12))
    hist_arr,pos_hist=obj.PolarAnalysis(select_res,Nframes,# sort=[1,2,3,4,5,6,7,8,0],
                                        zlim=zlim,control_plots=False,plot=True,ax=ax,Nbins=nbins)
    plt.close()
    assert isinstance(hist_arr, (list, np.ndarray)) 
    assert np.shape(hist_arr)==(np.shape(pos_hist)[1],2,nbins-1)
    assert isinstance(pos_hist, (list, np.ndarray))

def test_computeRG2D(surf_selection=surf_selection):
    u=mda.Universe(MD_NOWATER_TPR, MD_TRAJ)
    obj = BioPolymer2D(u,surf_selection=surf_selection)
    masses=obj.atom_group.atoms.masses
    total_mass=np.sum(masses)
    rgs=obj.computeRG2D(masses,total_mass)
    assert isinstance(rgs, (list, np.ndarray))

def test_getRgs2D(surf_selection=surf_selection):
    u=mda.Universe(MD_NOWATER_TPR, MD_TRAJ)
    obj = BioPolymer2D(u,surf_selection=surf_selection)
    rgs=obj.getRgs2D()
    assert np.shape(rgs)==(obj.universe.trajectory[-1].frame,4)

@pytest.mark.parametrize("plot", ['data', 'mean', 'both'])
@pytest.mark.parametrize("system_name", ['test_system',None])
@pytest.mark.parametrize("system_label", ['test_system1',None])
@pytest.mark.parametrize("legend", [True, False])
def test_RgPerpvsRgsPar(plot,system_name,system_label,legend,surf_selection=surf_selection):
    u=mda.Universe(MD_NOWATER_TPR, MD_TRAJ)
    obj = BioPolymer2D(u,surf_selection=surf_selection)
    rgs=obj.getRgs2D()
    color='tab:green'
    obj.system_name=system_name
    #if obj.system_name is None and systen_label is None, the function should pass if raise a KeyError
    if legend is True:
        if obj.system_name is None and system_label is None:
            with pytest.raises(KeyError):
                rg_ratio = obj.RgPerpvsRgsPar(rgs, color, plot=plot, legend=legend)
        else:
            # obj.system_name=system_name
            # assert obj.system_name!=system_label
            rg_ratio=obj.RgPerpvsRgsPar(rgs,color,plot=plot,system_label=system_label,legend=legend)
            assert isinstance(rg_ratio,(float))
    else:
        rg_ratio=obj.RgPerpvsRgsPar(rgs,color,plot=plot,legend=legend)
        assert isinstance(rg_ratio,(float))
    plt.close()

def test_getKDEAnalysis(surf_selection=surf_selection):
    u=mda.Universe(MD_NOWATER_TPR, MD_TRAJ)
    obj = BioPolymer2D(u,surf_selection=surf_selection)
    obj.getPositions()
    zlim=14
    Nframes=100
    paths=obj.getKDEAnalysis(zlim,Nframes)
    assert isinstance(paths, list,)

@pytest.mark.parametrize("lvl", [0,-1])
@pytest.mark.parametrize("getTotal", [True, False])
def test_getAreas(lvl,getTotal,surf_selection=surf_selection):
    u=mda.Universe(MD_NOWATER_TPR, MD_TRAJ)
    obj = BioPolymer2D(u,surf_selection=surf_selection)
    obj.getPositions()
    zlim=14
    Nframes=100
    paths=obj.getKDEAnalysis(zlim,Nframes)
    areas=obj.getAreas(paths,contour_lvl=lvl,getTotal=getTotal)
    if getTotal:
        assert isinstance(areas, (float))
    else:
        assert isinstance(areas, (list, np.ndarray))

# Create a matplotlib axis for testing
fig, ax = plt.subplots()
@pytest.mark.parametrize("ax", [ax,None])
@pytest.mark.parametrize("getNframes", [True, False])
def test_KDEAnalysisSelection(ax,getNframes,):
    u = mda.Universe(MD_NOWATER_TPR, MD_TRAJ)
    obj = BioPolymer2D(u, surf_selection=surf_selection)
    
    # Define test parameters
    select_res = 'resid 198 200 12 8 40 45 111 115 173'
    Nframes = 100
    zlim = 14
    legend = True
    plot_COM = True
    # Call the function with test parameters
    result = obj.KDEAnalysisSelection(
        select_res=select_res,
        Nframes=Nframes,
        zlim=zlim,
        ax=ax,
        show=False,
        legend=legend,
        plot_COM=plot_COM,
        getNframes=getNframes
    )

    # Validate the output
    if getNframes:
        paths_arr_arr, res, n_used_frames = result
        assert isinstance(paths_arr_arr, list)
        assert isinstance(res, mda.AtomGroup)
        assert isinstance(n_used_frames, int)
    else:
        paths_arr_arr, res = result
        assert isinstance(paths_arr_arr, list)
        assert isinstance(res, mda.AtomGroup)

    # Ensure the paths array is consistent with the residues
    assert len(paths_arr_arr) == len(res.residues)

    # Close the plot to avoid resource warnings
    plt.close(fig)

from MDAnalysis.exceptions import SelectionError

@pytest.mark.parametrize("selection1, selection2, expected_exception", [
    ("resname DOL", "protein", None),  # Valid selections, no exception expected
    ("name CA", "protein", SelectionError),  # Invalid selection2, should raise SelectionError
    ("protein", "name CA", SelectionError),  # Invalid selection1, should raise SelectionError
    ("name CA", "resname DOL and prop z <16 and prot z > 5", SelectionError),  # Both selections invalid, should raise SelectionError
])
@pytest.mark.parametrize("inplace", [True, False])
def test_getHbonds(selection1, selection2, expected_exception, inplace):
    u = mda.Universe(MD_NOWATER_TPR, MD_TRAJ)
    obj = BioPolymer2D(u, surf_selection=surf_selection)

    if expected_exception:
        with pytest.raises(expected_exception):
            obj.getHbonds(selection1=selection1, selection2=selection2)
    else:
        if inplace:
            obj.getHbonds(selection1=selection1, selection2=selection2)
            assert hasattr(obj, "hbonds")  # Ensure the results object has expected attributes
            assert hasattr(obj.hbonds, "hbonds")  # Ensure the results object has expected attributes
        else:
            results = obj.getHbonds(selection1=selection1, selection2=selection2,inplace=inplace)
            assert results is not None


def test_HbondsPerResidues():
    u = mda.Universe(MD_NOWATER_TPR, MD_TRAJ)
    obj = BioPolymer2D(u, surf_selection=surf_selection)

    # Ensure that getHbonds is called before HbondsPerResidues
    obj.getHbonds(selection1="resname DOL", selection2="protein", inplace=True)
    df = obj.HbondsPerResidues()
    assert isinstance(df, pd.DataFrame)  # Ensure the result is a DataFrame
    assert "ResIDs" in df.columns  # Check for expected columns
    assert "ResNames" in df.columns
    assert "Count" in df.columns
fig,ax=plt.subplots()
@pytest.mark.parametrize("ax", [ax,None])
def test_plotHbondsPerResidues(ax):
    u = mda.Universe(MD_NOWATER_TPR, MD_TRAJ)
    obj = BioPolymer2D(u, surf_selection=surf_selection)
    obj.getPositions()
    obj.getKDEAnalysis(zlim=14, Nframes=100)
    # Ensure that getHbonds is called before HbondsPerResidues
    obj.getHbonds(selection1="resname DOL", selection2="protein", inplace=True)

    
    sorted_df = obj.plotHbondsPerResidues(paths_for_contour=obj.kdeanalysis.paths,ax=ax)
    assert isinstance(sorted_df, pd.DataFrame)  # Ensure the result is a DataFrame
    assert "ResIDs" in sorted_df.columns  # Check for expected columns
    assert "ResNames" in sorted_df.columns
    assert "Count" in sorted_df.columns
    # ax = obj.plotHbondsPerResidues()
    # assert ax is not None  # Ensure the plot is created

    # Close the plot to avoid resource warnings
    plt.close()

# fig,ax=plt.subplots()
# @pytest.mark.parametrize("ax", [ax,None])
# def test_plotPathsInLevel(ax):
#     u = mda.Universe(MD_NOWATER_TPR, MD_TRAJ)
#     obj = BioPolymer2D(u, surf_selection=surf_selection)
#     obj.getPositions()
#     obj.getKDEAnalysis(zlim=14, Nframes=100)
#     contour_lvl=0

#     # Call the function without providing an axis
#     obj.plotPathsInLevel(paths=obj.kdeanalysis.paths, contour_lvl=contour_lvl, ax=ax)

#     # Ensure the function creates a new axis if none is provided
#     assert ax is None, "The function did not create an axis when ax=None."
#     assert isinstance(ax, plt.Axes), "The created object is not a valid matplotlib Axes."

#     # Ensure the axis contains plotted elements
#     assert len(ax.lines) > 0 or len(ax.collections) > 0, "No paths were plotted on the axis."

#     # Close the plot to avoid resource warnings
#     plt.close(ax.figure)
