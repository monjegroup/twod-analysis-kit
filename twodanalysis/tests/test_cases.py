import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pytest
import time


from twodanalysis import Memb2D
#import imageio
import os
import time
from twodanalysis.data.files import MEMBRANE_TPR, MEMBRANE_XTC





tpr = MEMBRANE_TPR
traj = MEMBRANE_XTC

# Creating the class

@pytest.fixture
def universe():
    return mda.Universe(tpr, traj)







#lipid_list = list(membrane.lipid_list)
#first_lipids = membrane.first_lipids




######### Lipid order 2d code related ########
"""
layers = ["top", "bot", "both"]
lipid_list.remove("CHL1")
nbins = 50
lipids = membrane.chain_info

for layer in layers:
    for key in lipid_list:
        start_time = time.time()
        H, edges = membrane.order_histogram(key, layer, nbins, lipids[key], start=60, final=100)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time: {elapsed_time:.6f} seconds")


        print(key, layer, nbins, lipids[key], 0, 180)
        plt.imshow(H,cmap = "Spectral", extent = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]])
        plt.colorbar(cmap = "Spectral")
        plt.savefig(f"{key}_test1_{layer}.png")
        plt.show()
        plt.close()
        start_time = time.time()
        H1, edges = membrane.order_histogram(key, layer, nbins, lipids[key], start=60, final=100, method = "numpy")
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time: {elapsed_time:.6f} seconds")

        print(key, layer, nbins, lipids[key], 0, 180)
        plt.imshow(H1,cmap = "Spectral", extent = [edges[0], edges[-1], edges[0], edges[-1]])
        plt.colorbar(cmap = "Spectral")
        plt.savefig(f"{key}_test2_{layer}.png")

        plt.show()
        plt.close()

        plt.scatter(H.flatten(), H1.flatten())
        plt.show()




layer = "top"
"""
nbins = 50
def test_all_order(universe):

    membrane = Memb2D(universe,
                        verbose = True,
                        add_radii = True)

    mat_top, edges = membrane.all_lip_order("top", nbins,
                        start = 0,
                        final = 100,
                        step = 1)#


    assert isinstance(mat_top, np.ndarray)

#plt.imshow(mat_top,extent=edges, cmap = "Spectral")
#plt.show()
#plt.imshow(mat_bot, extent=edges,cmap = "Spectral")
#plt.show()
#plt.imshow(mat_both,extent=edges, cmap = "Spectral")
#plt.show()



##### Packing deffects related ######


# Adding POPE lipids that are not taken into account
#membrane.non_polar_dict["POPE"].append("H101")
#membrane.non_polar_dict["POPE"].append("H91")

""" Print selections to test packing deffects with VMD
lipid_polar = {}
for key in first_lipids.keys():
    polar_dict = {"polar":[], "nonpolar": []}
    #print(membrane.non_polar_dict[key])
    atoms = membrane.u.select_atoms(f"resid {first_lipids[key]}")
    names = atoms.names
    for name in names:
        #print(membrane.non_polar_dict[key])
        if name in membrane.non_polar_dict[key]:
            polar_dict["nonpolar"].append(name)
        else:
            polar_dict["polar"].append(name)

    lipid_polar[key] = polar_dict


for key in lipid_polar.keys():
    print(key)
    print(f"###### resname {key} and {membrane.build_name(  lipid_polar[key]['polar'])}    ######\n\n")
    print(f"###### resname {key} and {membrane.build_name(lipid_polar[key]['nonpolar'])}    ######\n\n")

"""
# Data for RNA membrane
traj = "membrane.xtc"
tpr = "../data/veamos.tpr"
universe = mda.Universe(tpr, traj)
membrane = Memb2D(universe,
                        verbose = True,
                        add_radii = True)

defects, defects_dict = membrane.packing_defects(layer = "top", nbins = 200, periodic=True)

print(membrane.print_dict(defects_dict))

# Plot defects
#plt.imshow(defects, cmap = "viridis")
#plt.xlabel("x  $[\AA]$")
#plt.ylabel("y  $[\AA]$")
#plt.show()
#membrane.visualize_polarity()

#print(membrane.non_polar_dict["POPE"])

#membrane.visualize_polarity()

#membrane.non_polar_dict["POPE"].append("H91")
#membrane.non_polar_dict["POPE"].append("H101")
#membrane.non_polar_dict["POPI15"].append("H91")
#membrane.non_polar_dict["POPI15"].append("H101")
#membrane.non_polar_dict["POPI24"].append("H91")
#membrane.non_polar_dict["POPI24"].append("H101")
#membrane.visualize_polarity()
"""
vmin, vmax = membrane.guess_minmax_space()
print(vmin, vmax)
#membrane.v_max = vmax
#membrane.v_min = vmin

count = 0
filenames = []
positions = []
protein = membrane.u.select_atoms("resid 1-180")
for ts in membrane.u.trajectory[51::10]:
    positions.append(protein.positions)

positions = np.concatenate(positions, axis = 0)



print("fdsfdsfd################################3", positions.shape)


plt.plot(positions[:,0], positions[:,1])
plt.xlim(0, 180)
plt.ylim(0, 180)
plt.show()
plt.close()


percentage = []
for ts in membrane.u.trajectory[51:500:2]:

    start_time = time.time()
    matrix, return_dict = membrane.packing_defects(
                                    nbins = 150,
                                    layer = "bot",
                                    height = False,
                                    periodic=True,
                                    area = True,
                                    edges = [80,120,40,80]
                                    )
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.6f} seconds")


    edges = return_dict["edges"]



    area = f" deffects: {return_dict['area']['deffects']}, total {return_dict['area']['total']} "
    area_p = f"percentage : {return_dict['area']['deffects']/return_dict['area']['total']}"
    percentage.append(return_dict['area']['deffects']/return_dict['area']['total'])
    #fig,ax = plt.subplots(1,2)
    #ax[0].imshow(np.rot90(matrix), extent = [edges[0], edges[1], edges[0], edges[1]])
    #ax[1].imshow(np.rot90(matrix_height), extent = [edges[0], edges[1], edges[0], edges[1]])
    #ax[0].set_title(f"Frame {count}, area {area_p}")
    #plt.savefig(f"frame_{count}_periodic.png")
    #plt.close()

    filenames.append(f"frame_{count}_periodic.png")
    count += 1



with imageio.get_writer("gif_packing.gif", mode = "I", duration = 2) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
    plt.close()

for filename in filenames:
    os.remove(filename)

plt.close()
plt.clf()
data = pd.DataFrame()
x = list(range(len(percentage)))
data["percentage"] = percentage
data["frame"] = x
data["frame"] = data["frame"]*2
data["rolling"] = data["percentage"].rolling(50, center=True).mean()

plt.plot(data["frame"], data["percentage"])
plt.plot(data["frame"], data["rolling"])
plt.show()
plt.close()

"""


#### Membrane thickness related code ##########

#mat_both, edges = membrane.all_lip_order("both", nbins,
#                        start = 0,
#                        final = 100,
#                        step = 1)

#mat_thi, edges = membrane.thickness(50, start = 0, final = 100, step = 1)
#plt.scatter(mat_thi.flatten(), mat_both.flatten())
#plt.colorbar(cmap = "Spectral")
#plt.show()
#plt.close()
#plt.scatter(mat_both.flatten(), mat_thi.flatten(), alpha = 0.5)
#plt.savefig("corr_thilip.png")
#plt.close()
#print(membrane.lipid_list)




###### Code to test apl ######


voronoi_dict = membrane.voronoi_apl(layer = "top")


colors = {"DODMA": "blue",
        "POPE": "red",
        "DSPC": "green",
        "CHL1": "yellow"
}

print(len(voronoi_dict["points"]), len(voronoi_dict["resnames"]))
print(membrane.lipid_list)
n_lipids = voronoi_dict["orig_len"]

plt.scatter(*voronoi_dict["points"][:n_lipids].T, c =[colors[lipid] for lipid in voronoi_dict["resnames"][:n_lipids]])
count = 0
for vertices in voronoi_dict["vertices"][:n_lipids]:

    temp = vertices
    temp = np.concatenate([temp, [temp[0]]])
    plt.plot(*temp.T, c = colors[voronoi_dict["resnames"][count]])
    count+=1
print(membrane.print_dict(voronoi_dict["apl"]))
plt.show()

xmin = membrane.v_min
xmax = membrane.v_max
ymin = membrane.v_min
ymax = membrane.v_max

print(voronoi_dict["points"].shape, len(voronoi_dict["areas"]))
apl, edges = membrane.map_voronoi(voronoi_dict["points"], voronoi_dict["areas"], 180, [xmin, xmax, ymin, ymax])

plt.imshow(apl, extent = edges, cmap = "Spectral")
plt.xlabel("$x [\AA]$")
plt.ylabel("$y [\AA]$")
plt.colorbar()
plt.show()


#membrane.map_voronoitest(voronoi_dict["vertices"], voronoi_dict["areas"], 100, [membrane.v_min, membrane.v_max, membrane.v_min, membrane.v_max])
#membrane.map_voronoi(voronoi_dict["points"], voronoi_dict["areas"], 300, [membrane.v_min, membrane.v_max, membrane.v_min, membrane.v_max])






resu = membrane.grid_apl(layer = "top", start = 100, final = 101, step = 1, lipid_list = None)
print(resu)

matrices = membrane.windows_apl(layer = "top",
                                start = 260,
                                final = 290,
                                step = 1,
                                w_size = 2,
                                lipid_list = None,
                                nbins = 150,
                                )

fig, ax = plt.subplots(1,11, figsize = (5*11, 5), sharex = True, sharey = True)
x_ref = matrices[0].flatten()
indices = []
for i, matrix in enumerate(matrices):
    print(len(matrices))
    y = matrix.flatten()
    if i <= 10:
        print(i)
        y = matrix.flatten()
        ax[i].imshow(matrix, cmap = "Spectral")
    indices.append(np.corrcoef(x_ref, y)[0,1])
plt.show()
plt.close()

#index = 4*(np.array(list(range(len(indices))))/5)
#lista = [0.9999999999999998, 0.6025176760541637, 0.3472441438579375, 0.26539688826701563, 0.2850464370677984, 0.1480777929465581, 0.04842089982626965, 0.08878377145501166, 0.17719296655312355, 0.2072503188050951, 0.11521773655176341, 0.08574894206908254]
#plt.plot(index,indices)
#plt.plot(list(range(len(lista))),lista)
#print(indices)
#plt.show()

