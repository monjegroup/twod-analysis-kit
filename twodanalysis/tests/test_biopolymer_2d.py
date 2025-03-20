
import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from twodanalysis import BioPolymer2D
from twodanalysis.data.files import MD_NOWATER_TPR, MD_TRAJ
import pytest

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
def test_FilterMinFrames(pos_type, inplace, select, Nframes,zlim):
    u=mda.Universe(MD_NOWATER_TPR, MD_TRAJ,)
    obj = BioPolymer2D(u,surf_selection='resname DOL and name O1 and prop z > 20') ## This is dependent on the system used for testing
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
