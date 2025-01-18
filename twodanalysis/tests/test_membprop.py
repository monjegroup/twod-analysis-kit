import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pytest



from twodanalysis import MembProp


from twodanalysis.data.files import MEMBRANE_TPR, MEMBRANE_XTC





tpr = MEMBRANE_TPR
traj = MEMBRANE_XTC

# Creating the class

@pytest.fixture
def universe():
    return mda.Universe(tpr, traj)

def test_chain_length():
    u = mda.Universe(tpr, traj)
    membrane = MembProp(u)
    membrane.guess_chain_lenght()
    chains = membrane.chain_info
    print(membrane.chain_info)

    expected_dict = {"DODMA" : [17,17],
                     "CHL1" : [-1,7],
                     "DSPC" : [17,17],
                     "POPE" : [15,17],}
    assert expected_dict == chains


def test_last_carbon():
    u = mda.Universe(tpr, traj)
    membrane = MembProp(u)
    membrane.guess_last_cs()
    membrane.guess_chain_lenght()

    lipids = membrane.lipid_list
    last_cs = [membrane.working_lip[lipid]["last_c"] for lipid in lipids if lipid != "CHL1" ]
    chains = membrane.chain_info
    expected = [[f"C3{chains[lipid][0] + 1}", f"C2{chains[lipid][1] + 1}"] for lipid in lipids if lipid != "CHL1" ]

    assert last_cs == expected



