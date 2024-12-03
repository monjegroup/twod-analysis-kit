import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.integrate import simps
import time
from matplotlib.patches import Patch


class analysis:
    def __init__(self, top, traj, start = 0,
                                final = -1,
                                step =1,
                                membrane = False,
                                protein = False,
                                rna = False,
                                lipid_list = None,
                                protein_resids = None,
                                rna_resids = None,
                                time_scale = None,
                                ddrange = None,
                                ):

        self.top = top # Set the topology
        self.traj = traj # Set the traj

        self.color = ['tab:blue', 'tab:orange', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink'] # Set a palette color if needed

        if lipid_list:
            membrane = True
        if membrane:
            self.lipid_list = lipid_list # Lipids in study


        self.working_lip = {
                                "CHL1" : {"head" :"O3", "charge" : 0},
                                "DODMA" : {"head" :"N1", "charge" : -0.21},
                                "DSPC" : {"head" :"P", "charge" : 1.1},
                                "POPE" : {"head" :"P", "charge" : 1.1},
                                "DOPS" : {"head" :"P", "charge" : 0.1},
                                "POPS" : {"head" :"P", "charge" : 0.1},
                                "DSPE" : {"head" :"P", "charge" : 1.3},
                                "DOPC" : {"head" :"P", "charge" : 1.3},
                                "DOPE" : {"head" :"P", "charge" : 1.3},
                                "POPI1" : {"head" :"P", "charge" : 1.3},
                                "POPI2" : {"head" :"P", "charge" : 1.3},
                            } #List of known lipids and lipids head people usually use to work


        self.test = True
#        self.rna_at = rna_at
#        self.rna_resid = rna_resid
        self.start = start
        self.final = final
        self.step = step
        self.lower_lim = None
        self.upper_lim = None
        self.ha_table = None
        self.p_table = None

        self.range = ddrange
#        self.working_lip = {
#            "CHL1" : "O3",
#            "DODMA" : "N1",
#            "DSPC" : "P",
#            "POPE" : "P",
#            "DOPS" : "P",
#            "POPS" : "P",
#            "DSPE" : "P",
#        }

#        self.charge_li = {'DSPC' : 1.1,
#		                    'POPE' : 1.1,
#		                    'DODMA' : -0.21,
#		                    'POPS' : 0.1,
#		                    'DOPS' : 0.1,
#		                    'DSPE' : 1.3
#                            }
        print("###### Correctly initialized class rna ##########################")
        print(f"With \ntopology: {self.top},  trajectory: {self.traj}")
        print(f"Start of the trajectory: {self.start}")
        print(f"Final of the trajectory: {self.final}")
        print(f"Step of the trajectory: {self.step}")

        print(f"Lipids under consideration to find the middle of the membrane {self.lipid_list}")
 #       print(f"Heavy atoms of RNA: {self.rna_at}")
        #print(f"Resids corresponding to RNA: {self.rna_resid}")
        self.u = mda.Universe(self.top, self.traj)
        print(f"Lenght of the trajectory:{len(self.u.trajectory)}")
        self.lentraj = len(self.u.trajectory)
        self.all_p = self.u.select_atoms(self.build_resname(self.lipid_list) + " and name P")



        if time_scale != None:
            print(f"Start of the trajectory: {self.start*time_scale}ns")
            print(f"Final of the trajectory: {self.final*time_scale}ns")
            print(f"Step of the trajectory: {self.step*time_scale}ns")
            print(f"Lenght of the trajectory:{time_scale*len(self.u.trajectory)}ns")






    @staticmethod
    def build_resname(resnames_list):
        string = " (resname " + resnames_list[0]

        for resname in resnames_list[1:]:
            string = string + " or resname " + resname

        string = string + ") "
        return string

    @staticmethod
    def build_name(resnames_list):
        string = " (name " + resnames_list[0]

        for resname in resnames_list[1:]:
            string = string + " or name " + resname

        string = string + ") "
        return string

    # Return a dictionary with the positions of the lipids heads over a period of time
    def surface(self,
                start = None,
                final = None,
                step = None,
                lipid = "DSPC",
                layer = 'top',
                filename = None, include_charge = False):

        if start == None:
            start = self.start
        if final == None:
            final = self.final
        if step == None:
            step = self.step

        if filename == None:
            filename = f"{lipid}_{layer}_{start}_{final}.dat"

        lipid_list = self.lipid_list


        print("######### Running surface function ########## ")
        print(f"We will compute the surface files for a membrane with there lipids {lipid_list}")
        print(f"Currently working on: {lipid}")
        print(f"Layer: {layer}")
        print(f"Writing under the name of {filename}")


        if layer == "top":
            sign = " > "
        elif layer == "bot":
            sign = " < "


        ##### Select all the P atoms to find the middle of the membrane
        all_p = self.all_p



        #### Loop over trajectory to find the lipids in the requested membrane
        pos_data = []
        for ts in self.u.trajectory[start:final:step]:
            positions = all_p.positions[:,2]
            mean_z = positions.mean()

            # Selects the lipid head and of the working lipid
            if layer == "both":
                selection_string = f"(resname {lipid} and name {self.working_lip[lipid]['head']})"
            else:
                selection_string = f"(resname {lipid} and name {self.working_lip[lipid]['head']}) and prop z {sign} {str(mean_z)}"


            # Find the positions of the P atoms
            atoms = self.u.select_atoms(selection_string)


            #### Check number of atoms
            n_at = atoms.n_atoms

            ### Get positions
            atom_pos = atoms.positions

            ### Get resids
            atom_resid = atoms.resids
            atom_resid = atom_resid[:,np.newaxis]

            atom_pos = np.concatenate((atom_pos, atom_resid), axis = 1)
            atom_pos[:,2] = np.abs(atom_pos[:,2]-z_mean)

            pos_data.append(atom_pos)



        pos_data = np.concatenate(pos_data, axis = 0)
        df_data = pd.DataFrame(pos_data, columns = ["x", "y", "z", "id"])
        df_data["id"] = df_data["id"].astype(int)
        if include_charge:
            df_data["charge"] = self.charge_li[lipid]
            df_data.to_csv(f"pd_{filename}", index = False)
        df_data.to_csv(f"{filename}", index = False)

        return df_data   # Maybe have to change, it does not make sense to return this



    # Return a dictionary with the data of surface
    def surface_list(self,
                start = None,
                final = None,
                step = None,
                lipids = ["DSPC", "DOPC"],
                layer = 'top',
                filename = None, include_charge = False):
        if start == None:
            start = self.start
        if final == None:
            final = self.final
        if step == None:
            step = self.step
        lipid_data_dict = {}
        for lipid in lipids:
            #try:
            #    filename = f"{lipid}__{layer}_{start}_{final}.dat"
            #    lipid_data_dict[lipid] = pd.read_csv(filename)
            #except:
            filename = f"{lipid}_{layer}_{start}_{final}.dat"
            lipid_data_dict[lipid] = self.surface(start = start,
                    final = final,
                    step = step,
                    lipid = lipid,
                    layer = layer,
                    filename = filename, include_charge = include_charge)
        return lipid_data_dict

################## Code related to 2D order parameters ##########################

    # Given a list of the shape [C, Hx, Hy, Hz], gets the vectos contecting teh hydrogen ti the C and compute the indivifual cos^2theta) for each carbon for each lipid in the original selection (must be called after individual_order_sn1 or ns2
    @staticmethod
    def get_individual(lista):
        angles = [] # Store the angles
        for i in (range(len(lista)-1)): # Accounts for variable number of list (Change if the carbon has or not double bonds)
            vectores = lista[i+1].positions - lista[0].positions # Hidrogen - Carbons; output of shape (n_lipids, 3)
            costheta = vectores[:,2]**2/np.linalg.norm(vectores, axis = 1)**2 # Compute the costheta^2
            angles.append(costheta) # dim (n_lipids,)
        angles = np.array(angles) # dim ((2 or 3),n_lipids)
        #print("angles", angles.shape)
        angles = np.mean(angles, axis = 0) # output is dim n_lipids, it means the cos^2(theta) or the Carbon passed for each lipid
        return angles




    # Get the cos^2(theta) for each carbon in the selection, for sn1
    def individual_order_sn1(self, sel, lipid, n_chain):
        # Define list to store the chain cos^2(theta)
        chains = []

        # Loop over carbons
        for i in range(n_chain):
            # Define selections for H and C in the chain
            print(f"Value of the chain {i} sn1")
            selections = [
                            f"name C3{i+2}",
                            f"name H{i+2}X and not name HX",
                            f"name H{i+2}Y and not name HY",
                            f"name H{i+2}Z and not name HZ"
                        ]
            #print(selections)


            # Define a list to store atoms
            lista = []

            for selection in selections:
                atoms = sel.select_atoms(selection)


                if atoms.n_atoms != 0:
                    lista.append(atoms)
            # Call get_individual that computes the cos^2(theta) for each carbon.
            chains.append(self.get_individual(lista))
            #print(i, self.get_individual(lista).shape, self.get_individual(lista))
        chains = np.array(chains) # Expect array of dim (n_chain, n_lipids)
        return chains





    # Get the cos^2(theta) for each carbon in the selection, for sn2

    def individual_order_sn2(self, sel, lipid, n_chain):
        # Define list to store the chain cos^2(theta)
        chains = []

        # Loop over carbons
        max_v = 0
        for i in range(n_chain):
            # Define selections for H and C in the chain
            #print(f"Value of the chain {i} sn2")
            selections = [
                            f"name C2{i+2}",
                            f"name H{i+2}R and not name HR",
                            f"name H{i+2}S and not name HS",
                            f"name H{i+2}T and not name HT"
                        ]
            if lipid == "POPE" or lipid == "POPS" or lipid == "POPI1" or lipid == "POPI2":
                if selections[0] == "name C29":
                    selections[1] = "name H91"
                if selections[0] == "name C210":
                    selections[1] = "name H101"
            # Define a list to store atoms
            lista = []

            for selection in selections:
                atoms = sel.select_atoms(selection)


                if atoms.n_atoms != 0:
                    lista.append(atoms)
            # Call get_individual that computes the cos^2(theta) for each carbon.
            angles = self.get_individual(lista)
            if len(angles) > max_v:
                max_v = len(angles)
            chains.append(angles)

        #for i in range(len(chains)):
        #    print(len(chains[i]))
        #    while len(chains[i]) < max_v:
        #        chain[i].append(np.nan)

        chains = np.array(chains) # Expect array of dim (n_chain, n_lipids)
        #print("Individual_order_sn2:",chains.shape )
        return chains







    @staticmethod
    def count_order(data, min_lenght, n_chain):
        columns = ["index"]
        carbons_sn2 = False
        try:
            carbons_sn1 = [f"sn1-{i+2}" for i in range(n_chain[0]) ]
            carbons_sn2 = [f"sn2-{i+2}" for i in range(n_chain[1])]
            columns = columns + carbons_sn1 + carbons_sn2
        except:
            carbons_sn1 = [f"sn1-{i+2}" for i in range(n_chain) ]
            columns = columns + carbons_sn1


        df = pd.DataFrame(data, columns = columns)
        result = []

        for i in range(min_lenght):
            temp = df[df["index"] == i]
            if len(temp) > 0 :
                sn1 = temp[carbons_sn1]
                sn1 = sn1.mean()
                sn1 = 1.5 * sn1 - 0.5
                sn1 = sn1.abs()
                sn1 = sn1.mean()
                final = sn1
                if carbons_sn2:
                    sn2 = temp[carbons_sn2]
                    sn2 = sn2.mean()
                    sn2 = 1.5 * sn2 - 0.5
                    sn2 = sn2.abs()
                    sn2 = sn2.mean()

                    final = (final + sn2)*0.5
                result.append([i,final])
            else:
                result.append([i,0])
        result = np.array(result)
        result = result[:,1]
        return result


    # Method to average vector to pseudovector program
    @staticmethod
    def average_vector(data, min_lenght):
        columns = ["index", "x", "y", "z"] # Data expected is an np array with columns ["index", "x", "y", "z"]

        df = pd.DataFrame(data, columns = columns)
        result = []

        for i in range(min_lenght):
            temp = df[df["index"] == i]
            if len(temp) > 0 :
                bin_vect = temp[columns[1:]]
                bin_vect = bin_vect.mean()
                result.append(bin_vect.to_list())
            else:
                result.append([np.nan, np.nan, np.nan])
        result = np.array(result)

        return result


















    # Computes the average vector for each bin, sample are the raw x,y positions and weights are the vectors related to the head
    def pseudohistogram2D(self,sample1, weights, bins = 10, v_min = None, v_max = None):
        if v_min == None:
            v_min = np.min(sample1)
        if v_max == None:
            v_max = np.max(sample1)

        #print(v_min, v_max)
        nbin = np.empty(2,np.intp)
        edges = 2*[None]

        for i in range(2):
            edges[i] = np.linspace(v_min, v_max, bins +1)
            nbin[i] = len(edges[i]) + 1

        Ncount = (tuple(np.searchsorted(edges[i], sample1[:,i], side = "right") for i in range(2)))

        for i in range(2):
            on_edge = (sample1[:,i] == edges[i][-1])
            Ncount[i][on_edge] -= 1


        xy = np.ravel_multi_index(Ncount, nbin)
        xy_test = xy.reshape(-1,1)

        xy_test = np.concatenate((xy_test, weights), axis = 1)
        hist = self.average_vector(xy_test, nbin.prod())
        nbin = (nbin[0], nbin[1], 3)
        hist = hist.reshape(nbin)
        hist = hist.astype(float, casting = "safe")
        core = 2*(slice(1,-1),)
        hist = hist[core]

        return hist, edges










    # Computes teh histogram of the average order parameters in each bin
    def histogram2D(self,sample1, weights, n_chain, bins = 10, v_min = None, v_max = None):
        if v_min == None:
            v_min = np.min(sample1)
        if v_max == None:
            v_max = np.max(sample1)

        #print(v_min, v_max)
        nbin = np.empty(2,np.intp)
        edges = 2*[None]

        for i in range(2):
            edges[i] = np.linspace(v_min, v_max, bins +1)
            nbin[i] = len(edges[i]) + 1

        Ncount = (tuple(np.searchsorted(edges[i], sample1[:,i], side = "right") for i in range(2)))

        for i in range(2):
            on_edge = (sample1[:,i] == edges[i][-1])
            Ncount[i][on_edge] -= 1


        xy = np.ravel_multi_index(Ncount, nbin)
        xy_test = xy.reshape(-1,1)

        xy_test = np.concatenate((xy_test, weights), axis = 1)

        hist = self.count_order(xy_test, nbin.prod(), n_chain)
        hist = hist.reshape(nbin)
        hist = hist.astype(float, casting = "safe")
        core = 2*(slice(1,-1),)
        hist = hist[core]

        return hist, edges







    def order_histogram(self, lipid, layer, n_grid,
                        n_chain,
                        v_min = None,
                        v_max = None,
                        all_p = None,
                        start = None,
                        final = None,
                        step = 1):

        if all_p == None:
            all_p = self.all_p
        if start == None:
            start = self.start
        if final == None:
            final = self.final

        if layer == "top":
            sign = " > "
        elif layer == "bot":
            sign = " < "

        try:
            n_chain1 = n_chain[0]
            n_chain2 = n_chain[1]
        except:
            n_chain1 = n_chain
            n_chain2 = 0

        matrix = [] # this will store a matrix of the shape (2+n_chain,

        for ts in self.u.trajectory[start:final:step]:
            z = all_p.positions[:,2]
            z_mean = z.mean() # get middel of the membrane
            #Pick atoms in the layer
            if layer == "both":
                layer = self.u.select_atoms(f"byres ((resname {lipid} and name {self.working_lip[lipid]['head']}))")
            else:
                layer = self.u.select_atoms(f"byres ((resname {lipid} and name {self.working_lip[lipid]['head']}) and prop z {sign} {z_mean})")
            #print("Info:", all_p.n_atoms, z_mean, layer.n_atoms)

            only_p = layer.select_atoms(f"name {self.working_lip[lipid]['head']}")
            positions = only_p.positions[:,:2]
            angles_sn1 = self.individual_order_sn1(layer, lipid, n_chain1)
            angles_sn1 = angles_sn1.T

            #print(angles_sn1.T.shape, positions.shape)
            #print(angles_sn1.shape, positions.shape)
            to_write = np.concatenate([positions, angles_sn1], axis = 1)
            if n_chain2 != 0:
                angles_sn2 = self.individual_order_sn2(layer, lipid, n_chain2)
                angles_sn2 = angles_sn2.T
                to_write = np.concatenate([to_write, angles_sn2], axis = 1)

            matrix.append(to_write) # Expect dim (n_lipids, 2+n_chain1+n_chain2)
            #print("Frame:",to_write.shape)

        #matrix = np.array(matrix) # Expect dim (frames, n_lipids, 2+n_chain1+n_chain2)
        matrix = np.concatenate(matrix, axis = 0) # Expect dim (n_lipids*frames, 2+n_chain1+n_chain2)


        H, edges = self.histogram2D(matrix[:,:2], matrix[:,2:], n_chain, bins = n_grid, v_min = v_min, v_max = v_max)
        H = np.rot90(H)
        H[H==0] = np.nan

        return H, edges






"""



    def surface(self, start = None, final = None, step = None,
                                    lipid = "DSPC",
                                    layer = 'top',
                                    filename = 'test.dat', include_charge = False):
#        print(lipid_list)

        if start == None:
            start = self.start
        if final == None:
            final = self.final
        if step == None:
            step = self.step



        lipid_list = self.lipid_list
        print("######### Running surface function ########## ")
        print(f"We will compute the surface files for {lipid_list}")
        print(f"Currently working on: {lipid}")
        print(f"Layer: {layer}")
        print(f"Writing under the name of {filename}")


        sign = " > "
        if layer != "top":
            sign = " < "
        ##### Select all the P atoms to find the middle of the membrane
        all_p = self.all_p



        #### Loop over trajectory to find the lipids in the requested membrane
        pos_data = []
        for ts in self.u.trajectory[start:final:step]:
            positions = all_p.positions[:,2]
            mean_z = positions.mean()
            selection_string = f"(resname {lipid} and name {self.working_lip[lipid]}) and prop z {sign} {str(mean_z)}"

            #print(selection_string)
            atoms = self.u.select_atoms(selection_string)


            #### Check number of atoms
            n_at = atoms.n_atoms

            ### Get positions
            atom_pos = atoms.positions

            ### Get resids
            atom_resid = atoms.resids
            atom_resid = atom_resid[:,np.newaxis]

            atom_pos = np.concatenate((atom_pos, atom_resid), axis = 1)
            atom_pos[:,2] = np.abs(atom_pos[:,2] - mean_z)

            pos_data.append(atom_pos)



        pos_data = np.concatenate(pos_data, axis = 0)
        if include_charge:
            print(pos_data.shape)
            df_data = pd.DataFrame(pos_data, columns = ["x", "y", "z", "id"])
            df_data["charge"] = self.charge_li[lipid]
            df_data.to_csv(f"pd_{filename}", index = False)
        np.savetxt(filename, pos_data, delimiter = ',')

        return self.build_resname(lipid_list)   # Maybe have to change, it does not make sense to return this



    def surface_vectors(self, start = None, final = None, step = None,
                                    lipid = "DSPC",
                                    layer = 'top',
                                    multiple_step=None,
                                    n_frames = 100,
                                    filename = None, include_charge = False):
#        print(lipid_list)

        if multiple_step:
            step1 = multiple_step[0]
            step2 = multiple_step[1]
            frames_list = []
            frames_list.append(start)
            for i in range(n_frames-1):
                if i % 2 == 0:
                    frames_list.append(frames_list[i] + step1)
                else:
                    frames_list.append(frames_list[i] + step2)



        if start == None:
            start = self.start
        if final == None:
            final = self.final
        if step == None:
            step = self.step

        sel_for_vec = {
                        "DOPS": ["name P", "name C13"],
                        "POPS": ["name P", "name N"], #May change N to C12
                        "DODMA": ["name C2", "name N1"],
                        "POPE": ["name P", "name N"],
                        "DSPC": ["name P", "name N"],
                        "DSPE": ["name P", "name N"], # May topdate to C2 from the PEGA that is 5 resid far from the original (44->49)
                        "CHL1": ["name C3", "name C17"],
                        }

        lipid_list = self.lipid_list


        sign = " > "
        if layer != "top":
            sign = " < "
        ##### Select all the P atoms to find the middle of the membrane
        all_p = self.all_p



        #### Loop over trajectory to find the lipids in the requested membrane
        pos_data = []
        count = start
        frames_check = 0

        for ts in self.u.trajectory[start:final:step]:
            positions = all_p.positions[:,2]
            mean_z = positions.mean()
            selection_string = f"byres (resname {lipid} and name {self.working_lip[lipid]}) and prop z {sign} {str(mean_z)}"

            #print(selection_string)
            atoms = self.u.select_atoms(selection_string)

            origin = atoms.select_atoms(sel_for_vec[lipid][0])
            fin = atoms.select_atoms(sel_for_vec[lipid][1])
            vectors = fin.positions - origin.positions # Expected dim n_lip X 3
            #print(vectors.shape)
            atoms = origin

            atoms_pos = atoms.positions

            vect_df = pd.DataFrame(vectors, columns = ["x", "y", "z"])
            #### Check number of atoms
            n_at = atoms.n_atoms


            ### Get resids
            atom_resid = atoms.resids
            vect_df["id"] = atom_resid
            vect_df["x_0"] = atoms_pos[:,0]
            vect_df["y_0"] = atoms_pos[:,1]
            vect_df["z_0"] = atoms_pos[:,2]

            if multiple_step:
                #print(count, frames_list)
                if count in frames_list:
                    pos_data.append(vect_df)
                    frames_check += 1
            else:
                pos_data.append(vect_df)
            count += 1

        if frames_check != 0:
            print(f"Frames has been written with double steps and the number of frames is: {frames_check}")
            print(f"The start frame is {frames_list[0]} and the last is {frames_list[-1]}, which includes a lenght of {len(frames_list)}")

        pos_data = pd.concat(pos_data, axis = 0)

        if filename == None:
            filename = f"{lipid}_{start}.dat"
        pos_data.to_csv(f"pd_{filename}", index = False)

        return pos_data   # Maybe have to change, it does not make sense to return this

    def plot_2dvectors(self, lipid = "DOPS",stage = "touchbot", start = 0, final = -1, layer = "top",filename = None):
        if filename == None:
            filename = f"pd_{lipid}_{start}.dat"
        try:
            vect = pd.read_csv(filename)
        except:
            vect = self.surface_vectors(lipid = lipid, start = start, final = final, layer = layer, filename = filename.replace("pd_", ""))


        #ha_com, p_pos = self.rna_positions(start = start, final = final)
        nts, ha_com = self.get_ordered_nts(f"../{stage}/")





        print(vect)
        plt.close()
        plt.quiver(vect["x_0"],
                    vect["y_0"],
                    vect["x"],
                    vect["y"],
                    angles = "xy",
                    scale_units = "xy",
                    color = "blue",
                    scale = 3)
        nts = nts.iloc[:6]
        if len(nts) > 0:
            plt.plot(nts['x'], nts['y'],marker = 'o', markerfacecolor = 'white',
					 markeredgecolor = 'black',
					 markeredgewidth = 2,
					 linestyle='None')

        sns.kdeplot(x = ha_com["x"], y = ha_com["y"], color = "black", levels = 1)
        plt.xlabel("x [$\\AA$]")
        plt.ylabel("y [$\\AA$]")
        plt.savefig(f"vectors_{lipid}_{start}.png")
        plt.close()





        vect1 = vect.groupby(by="id").mean()
        vect1["norm"] = np.linalg.norm(vect1[["x", "y", "z"]].values, axis =1)
        print("Grotoped by :", vect1)
        vect["norm"] = np.linalg.norm(vect[["x", "y", "z"]].values, axis =1)
        print("Descride vect:", vect.describe())
        print("Describe", vect1.describe())
        v = -1
        print("Plotted vectors ",vect1.iloc[:v])
        #plt.scatter(vect1["x_0"].iloc[:v],vect1["y_0"].iloc[:v], color = "red" )
        #plt.scatter(vect1["x_0"].iloc[:v]+vect1["x"].iloc[:v],vect1["y_0"].iloc[:v]+vect1["y"].iloc[:v], color = "blue" )
        plt.quiver(vect1["x_0"].iloc[:v],
                    vect1["y_0"].iloc[:v],
                    vect1["x"].iloc[:v],
                    vect1["y"].iloc[:v],
                    angles = "xy",
                    scale_units = "xy",
                    scale = 1,
                    color = "blue",
                    )
        nts = nts.iloc[:6]
        if len(nts) > 0:
            plt.plot(nts['x'], nts['y'],marker = 'o', markerfacecolor = 'white',
					 markeredgecolor = 'black',
					 markeredgewidth = 2,
					 linestyle='None')

        sns.kdeplot(x = ha_com["x"], y = ha_com["y"], color = "black", levels = 1)
        plt.xlabel("x [$\\AA$]")
        plt.xlim(0,175)
        plt.ylim(0,175)
        plt.ylabel("y [$\\AA$]")
        plt.savefig(f"avg_vectors_{lipid}_{start}.png")
        plt.close()





        sample1 = vect[["x_0", "y_0"]].values
        weights = vect[["x", "y", "z"]].values
        initial_time = time.time()
        f_hist, X, Y = np.histogram2d(x = vect["x_0"], y = vect["y_0"], weights = vect["x"], range =[[0,175],[0,175]], bins = 50 )
        f_hist_v, X, Y = np.histogram2d(x = vect["x_0"], y = vect["y_0"], range =[[0,175],[0,175]], bins = 50 )
        f_hist_v[f_hist_v == 0] = 1
        f_hist_x = f_hist/f_hist_v
        f_hist, X, Y = np.histogram2d(x = vect["x_0"], y = vect["y_0"], weights = vect["y"], range =[[0,175],[0,175]], bins = 50 )
        f_hist_v, X, Y = np.histogram2d(x = vect["x_0"], y = vect["y_0"], range =[[0,175],[0,175]], bins = 50 )
        f_hist_v[f_hist_v == 0] = 1
        f_hist_y = f_hist/f_hist_v
        f_hist, X, Y = np.histogram2d(x = vect["x_0"], y = vect["y_0"], weights = vect["z"], range =[[0,175],[0,175]], bins = 50 )
        f_hist_v, X, Y = np.histogram2d(x = vect["x_0"], y = vect["y_0"], range =[[0,175],[0,175]], bins = 50 )
        f_hist_v[f_hist_v == 0] = 1
        f_hist = f_hist/f_hist_v
        print("hist",f_hist)
        final_time = time.time()
        exect = final_time - initial_time
        print(exect)

        initial_time = time.time()
        hist, edges = self.pseudohistogram2D(sample1, weights, bins = 50, v_min = 0, v_max = 175)
        final_time = time.time()
        exect = final_time - initial_time





        print("pse:",exect)

        print("pseudohist",hist[:,:,2])
        edges = np.array(edges[0][:-1])
        print(hist.shape, len(edges))


        dx = edges[1]-edges[0]
        edges = edges + 0.5 * dx
        x, y = np.meshgrid(edges, edges)
        print(x.shape, y.shape)
        plt.close()
        plt.quiver(x,
                    y,
                    hist[:,:,0].T,
                    hist[:,:,1].T,
                    angles = "xy",
                    color = "blue",
                    scale_units = "xy",
                    scale = 1)
        nts = nts.iloc[:6]
        if len(nts) > 0:
            plt.plot(nts['x'], nts['y'],marker = 'o', markerfacecolor = 'white',
					 markeredgecolor = 'black',
                     alpha = 0.,
					 markeredgewidth = 2,
					 linestyle='None')

        sns.kdeplot(x = ha_com["x"], y = ha_com["y"], color = "black", levels = 1)
        plt.xlabel("x [$\\AA$]")
        plt.ylabel("y [$\\AA$]")
        plt.savefig(f"{lipid}_field_{start}.png")
        return vect



##### Plot vector fields and vector for vatious lipids at a time
    def plot_vector_lipids(self, lipids = ["DOPS"],stage = "touchbot", start = 0, final = -1,output = "field.png", layer = "top", zoom = False, multiple_step = None, n_frames = 100, axis1 = None, axis2 = None, axis3 = None):

        dict_data = {}
        print(multiple_step, n_frames)

        for lipid in lipids:
            try:
                dict_data[lipid] = pd.read_csv(f"pd_{lipid}_{start}.dat")
                dict_data[lipid]["norm"] = np.linalg.norm(dict_data[lipid][["x", "y", "z"]].values, axis =1)
                dict_data[lipid]["x"] = dict_data[lipid]["x"]/dict_data[lipid]["norm"]
                dict_data[lipid]["y"] = dict_data[lipid]["y"]/dict_data[lipid]["norm"]
                dict_data[lipid]["z"] = dict_data[lipid]["z"]/dict_data[lipid]["norm"]
            except:
                if multiple_step:
                    dict_data[lipid] = self.surface_vectors(lipid = lipid,
                                                            start = start,
                                                            final = final,
                                                            layer = layer,
                                                            multiple_step=multiple_step,
                                                            n_frames = n_frames)
                else:
                    dict_data[lipid] = self.surface_vectors(lipid = lipid, start = start, final = final, layer = layer)
                dict_data[lipid]["norm"] = np.linalg.norm(dict_data[lipid][["x", "y", "z"]].values, axis =1)
                dict_data[lipid]["x"] = dict_data[lipid]["x"]/dict_data[lipid]["norm"]
                dict_data[lipid]["y"] = dict_data[lipid]["y"]/dict_data[lipid]["norm"]
                dict_data[lipid]["z"] = dict_data[lipid]["z"]/dict_data[lipid]["norm"]
        nts, ha_com = self.get_ordered_nts(f"../{stage}/", start = start, final = final)
        color_lip = {"DSPC": "blue", "POPE":"red", "DOPS":"green", "POPS":"black", "DSPE":"orange", "DODMA":"gray"}
        colors = ["blue", "red", "green", "orange"]
        v_min = None
        v_max = None
        markerwidth = 2
        markersize = 5
        if zoom:
            v_min = min(ha_com["x"].min(), ha_com["y"].min()) - 5
            v_max = max(ha_com["x"].max(), ha_com["y"].max()) + 5
            markerwidth = 4
            markersize = 15
        fig, ax_array = plt.subplots(1,3,figsize=(30,10), sharey = True)

        ##### Plot raw postitions with vectors ####
        plt.close()
        count = 0

        for lipid in lipids:
            vect = dict_data[lipid].copy()

            if axis1:
                print("imporimiendo en el axis")
                axis1.quiver(vect["x_0"],
                    vect["y_0"],
                    vect["x"],
                    vect["y"],
                    angles = "xy",
                    scale_units = "xy",

                    label = lipid,
                    color = color_lip[lipid],
                    scale = 1)
            else:

                plt.quiver(vect["x_0"],
                    vect["y_0"],
                    vect["x"],
                    vect["y"],
                    angles = "xy",
                    scale_units = "xy",
                    color = color_lip[lipid],
                    label = lipid,
                    scale = 1)
                ax_array[0].quiver(vect["x_0"],
                    vect["y_0"],
                    vect["x"],
                    vect["y"],
                    angles = "xy",
                    scale_units = "xy",
                    label = lipid,
                    color = color_lip[lipid],
                    scale = 1)
            count += 1

        nts = nts.iloc[:6]
        if len(nts) > 0:
            if not axis1:
                plt.plot(nts['x'], nts['y'],marker = 'o', markerfacecolor = 'white',
					 markeredgecolor = 'black',
                     markersize = markersize,
					 markeredgewidth = markerwidth,
					 linestyle='None')
                ax_array[0].plot(nts['x'], nts['y'],marker = 'o', markerfacecolor = 'white',
					 markeredgecolor = 'black',
                     markersize = markersize,
					 markeredgewidth = markerwidth,
					 linestyle='None')
            if axis1:
                axis1.plot(nts['x'], nts['y'],marker = 'o', markerfacecolor = (0., 0., 0., 0.),
					 markeredgecolor = 'black',
                     markersize = markersize/4,
					 markeredgewidth = markerwidth/4,
                     linewidth = 0.5,
					 linestyle='None')


        if axis1:
            sns.kdeplot(x = ha_com["x"], y = ha_com["y"],ax = axis1, color = "black", levels = 1)
            if v_min != None and v_max != None:
                axis1.set_xlim(v_min, v_max)
                axis1.set_ylim(v_min, v_max)
                axis1.set_xlim(v_min,v_max)
                axis1.set_ylim(v_min,v_max)
        else:
            sns.kdeplot(x = ha_com["x"], y = ha_com["y"], color = "black", levels = 1)
            sns.kdeplot(x = ha_com["x"], y = ha_com["y"],ax = ax_array[0], color = "black", levels = 1)
        plt.xlabel("x [$\\AA$]")
        plt.ylabel("y [$\\AA$]")
        ax_array[0].set_xlabel("x [$\\AA$]")
        ax_array[0].set_ylabel("y [$\\AA$]")
        ax_array[0].set_title("Raw points")
        plt.title("Raw points")
        if v_min != None and v_max != None:
            plt.xlim(v_min, v_max)
            plt.ylim(v_min, v_max)
            ax_array[0].set_xlim(v_min,v_max)
            ax_array[0].set_ylim(v_min,v_max)

        #plt.legend(bbox_to_anchor=(0., 1.02,1, .102), loc='lower left', borderaxespad=0.,mode = "expand",frameon=False, ncols = len(lipids))
        plt.savefig(f"vectors_{start}.png")
        plt.show()
        plt.close()

        return axis1
        ###### Plot mean plot ######
        '''
        plt.close()
        count = 0
        for lipid in lipids:
            vect = dict_data[lipid].copy()
            vect = vect.groupby(by="id").mean()
            plt.quiver(vect["x_0"],
                    vect["y_0"],
                    vect["x"],
                    vect["y"],
                    angles = "xy",
                    scale_units = "xy",
                    color = colors[count],
                    label = lipid,
                    scale = 0.25)
            ax_array[1].quiver(vect["x_0"],
                    vect["y_0"],
                    vect["x"],
                    vect["y"],
                    angles = "xy",
                    scale_units = "xy",
                    color = colors[count],
                    label = lipid,
                    scale = 0.25)
            count += 1
        nts = nts.iloc[:6]
        if len(nts) > 0:
            plt.plot(nts['x'], nts['y'],marker = 'o', markerfacecolor = 'white',
					 markeredgecolor = 'black',
                     markersize = markersize,
					 markeredgewidth = markerwidth,
					 linestyle='None')
            ax_array[1].plot(nts['x'], nts['y'],marker = 'o', markerfacecolor = 'white',
					 markeredgecolor = 'black',
                     markersize = markersize,
					 markeredgewidth = markerwidth,
					 linestyle='None')

        sns.kdeplot(x = ha_com["x"], y = ha_com["y"], color = "black", levels = 1)
        sns.kdeplot(x = ha_com["x"], y = ha_com["y"],ax =ax_array[1] , color = "black", levels = 1)
        plt.xlabel("x [$\\AA$]")
        plt.ylabel("y [$\\AA$]")
        plt.title("Mean")
        ax_array[1].set_xlabel("x [$\\AA$]")
        ax_array[1].set_ylabel("y [$\\AA$]")
        ax_array[1].set_title("Mean")
        if v_min != None and v_max != None:
            plt.xlim(v_min, v_max)
            plt.ylim(v_min, v_max)
            ax_array[1].set_xlim(v_min,v_max)
            ax_array[1].set_ylim(v_min,v_max)
        plt.legend(bbox_to_anchor=(0., 1.02,1, .102), loc='lower left', borderaxespad=0.,mode = "expand",frameon=False, ncols = len(lipids))
        plt.savefig(f"avg_{start}.png")
        plt.show()
        plt.close()

        ###### Plot field #####
        count = 0
        for lipid in lipids:
            vect = dict_data[lipid].copy()
            sample1 = vect[["x_0", "y_0"]].values
            weights = vect[["x", "y", "z"]].values
            hist, edges = self.pseudohistogram2D(sample1, weights, bins = 50, v_min = 0, v_max = 175)



            edges = np.array(edges[0][:-1])


            dx = edges[1]-edges[0]
            edges = edges + 0.5 * dx
            x, y = np.meshgrid(edges, edges)
            plt.quiver(x,
                    y,
                    hist[:,:,0].T,
                    hist[:,:,1].T,
                    angles = "xy",
                    color = colors[count],
                    label = lipid,
                    scale_units = "xy",
                    scale = 0.25)
            ax_array[2].quiver(x,
                    y,
                    hist[:,:,0].T,
                    hist[:,:,1].T,
                    angles = "xy",
                    color = colors[count],
                    label = lipid,
                    scale_units = "xy",
                    scale = 0.25)
            count += 1
        plt.legend(bbox_to_anchor=(0., 1.02,1, .102), loc='lower left', borderaxespad=0.,mode = "expand",frameon=False, ncols = len(lipids))
        nts = nts.iloc[:6]
        if len(nts) > 0:
            plt.plot(nts['x'], nts['y'],marker = 'o', markerfacecolor = 'white',
					 markeredgecolor = 'black',
                     markersize = markersize,
					 markeredgewidth = markerwidth,
					 linestyle='None')
            ax_array[2].plot(nts['x'], nts['y'],marker = 'o', markerfacecolor = 'white',
					 markeredgecolor = 'black',
                     markersize = markersize,
					 markeredgewidth = markerwidth,
					 linestyle='None')

        sns.kdeplot(x = ha_com["x"], y = ha_com["y"], color = "black", levels = 1)
        sns.kdeplot(x = ha_com["x"], y = ha_com["y"], ax = ax_array[2], color = "black", levels = 1)
        plt.xlabel("x [$\\AA$]")
        plt.ylabel("y [$\\AA$]")
        plt.title("Field (cumulative)")
        ax_array[2].set_xlabel("x [$\\AA$]")
        ax_array[2].set_ylabel("y [$\\AA$]")
        ax_array[2].set_title("Field (cumulative)")
        if v_min != None and v_max != None:
            plt.xlim(v_min, v_max)
            plt.ylim(v_min, v_max)
            ax_array[2].set_xlim(v_min,v_max)
            ax_array[2].set_ylim(v_min,v_max)
        plt.savefig(f"field_{start}.png")
        plt.show()

        fig.savefig(f"{start}_output")
        fig.show()

        return axis1
        '''

    def test_func(self):
        print(type(self.u.trajectory))









    def height_matrix(self, lipids, layer,start = None, final = None, step = None, nbins = 50):

        if start == None:
            start = self.start
        if final == None:
            final = self.final
        if step == None:
            step = self.step


        print(f"Computing matrix for {layer} in frames {start}-{final}")
        data = []
        for lipid in lipids:
            filename = f"{lipid}_{layer}_{start}-{final}.dat"
            try:
                df_data = pd.read_csv(f"pd_{filename}")
            except:
                self.surface(lipid = lipid, layer = layer, filename = filename, include_charge = True, start = start, final = final)
                df_data = pd.read_csv(f"pd_{filename}")
            data.append(df_data)

        data = pd.concat(data, axis = 0)
        print(data)
        xmin = data["x"].min()
        xmax = data["x"].max()
        ymin = data["y"].min()
        ymax = data["y"].max()

        xmin = 0
        xmax = 175
        ymin = 0
        ymax = 175

        H_height, x_edges, y_edges = np.histogram2d(x = data["x"], y = data["y"], weights = data["z"], bins = nbins, range = [[xmin,xmax], [ymin,ymax]])
        H_count, x_edges, y_edges = np.histogram2d(x = data["x"], y = data["y"], bins = nbins, range = [[xmin, xmax], [ymin,ymax]])

        H_count[H_count == 0] = 1.

        H_avg = H_height/H_count

        H_avg[H_avg == 0] =  np.nan

        H_avg = np.rot90(H_avg)

        np.savetxt(f'Height_{layer}_{start}_{final}.dat', H_avg, fmt = '%.2f')
        np.savetxt(f"edges_{layer}_{start}_{final}.dat", x_edges, fmt = "%.2f")
        return H_avg, x_edges, y_edges

    def charge_matrix(self, lipids, layer,start = None, final = None, step = None, nbins = 50):

        if start == None:
            start = self.start
        if final == None:
            final = self.final
        if step == None:
            step = self.step


        print(f"Computing matrix for {layer} in frames {start}-{final}")
        data = []
        for lipid in lipids:
            filename = f"{lipid}_{layer}_{start}-{final}.dat"
            try:
                df_data = pd.read_csv(f"pd_{filename}")
            except:
                self.surface(lipid = lipid, layer = layer, filename = filename, include_charge = True, start = start, final = final)
                df_data = pd.read_csv(f"pd_{filename}")
            data.append(df_data)

        data = pd.concat(data, axis = 0)
        print(data)
        xmin = data["x"].min()
        xmax = data["x"].max()
        ymin = data["y"].min()
        ymax = data["y"].max()

        xmin = 0
        xmax = 175
        ymin = 0
        ymax = 175

        H_height, x_edges, y_edges = np.histogram2d(x = data["x"], y = data["y"], weights = data["charge"], bins = nbins, range = [[xmin,xmax], [ymin,ymax]])
        H_count, x_edges, y_edges = np.histogram2d(x = data["x"], y = data["y"], bins = nbins, range = [[xmin, xmax], [ymin,ymax]])

        H_count[H_count == 0] = 1.

        H_avg = H_height/H_count

        H_avg[H_avg == 0] =  np.nan

        H_avg = np.rot90(H_avg)

        np.savetxt(f'H{layer}_{start}_{final}.dat', H_avg, fmt = '%.2f')
        return H_avg, x_edges, y_edges

    def plot_matrix(self, H, x_edges, y_edges,stage = "highest",colorbarlabel = "Height $\\AA$", filename = "deltaplot.png"):
        #rcParams['font.family']='serif'
        bigsize = 17
        #plt.figure(figsize = (5,2))
        plt.rc('axes', labelsize=bigsize)
        plt.rc('xtick', labelsize=bigsize)
        plt.rc('ytick', labelsize=bigsize)
        plt.rc('legend', fontsize=13)
        #plt.rc('xtick', labelsize = 28)
        vmin = np.nanmin(H)
        vmax = np.nanmax(H)
        fig, ax = plt.subplots()

        xmin = x_edges.min()
        xmax = x_edges.max()
        ymin = y_edges.min()
        ymax = y_edges.max()

        plott = ax.imshow(H, vmin = vmin,
                vmax=vmax,
            extent = [xmin,xmax,ymin,ymax],
            cmap = 'Spectral')

        cbar = fig.colorbar(plott, ax=ax)
        ticklabs = cbar.ax.get_yticklabels()

        try:
            nts, ha_com = self.get_ordered_nts(f"../{stage}/")
            print(f" File plotting ../{stage}/full_ha_table.dat")
        except:
            ha_com, p_pos  = self.rna_positions(ha_file = f"ha_com.dat", p_file = f"p_pos.dat")

        try:
            nts = nts.iloc[:6]
            if len(nts) > 0:
                plt.plot(nts['x'], nts['y'],marker = 'o', markerfacecolor = 'white',
					 markeredgecolor = 'black',
					 markeredgewidth = 2,
					 linestyle='None')
        except:
            pass

        sns.kdeplot(x = ha_com["x"], y = ha_com["y"], ax = ax, color = "blue", levels = 1)
        cbar.set_label(colorbarlabel)
        plt.xlabel('x [$\\AA$]')
        plt.ylabel('y [$\\AA$]')
        plt.tight_layout()
        plt.show()
        plt.savefig(filename)


    def get_ordered_nts(self, filename, start = None, final = None):
        try:
            data_df = pd.read_csv(f"{filename}full_ha_table.dat")
        except:
            print(f"There is not data in {filename} or {filename} does not exist. In this case, nts does not exist. We will compute rna_positions and return this value")
            ha_com, p_pos = self.rna_positions(start = start, final = final)
            return pd.DataFrame(), ha_com
            #raise ValueError('A very specific bad thing happened.')

        data_df = data_df.iloc[:6]
        resids = data_df["resname"].tolist()
        ha_com, p_pos  = self.rna_positions(ha_file = f"{filename}ha_com.dat", p_file = f"{filename}p_pos.dat")
        #print(ha_com.columns)
        means = ha_com.groupby("resid").mean()
        return means.loc[resids], ha_com










    def min_dist(self,
                ha_filename = "closest_nts_ha.dat",
                p_filename = "closest_nts_p.dat"):

        print("######### Currently working on min distance computations ####################")

        ####### If data exists provided, data is already there, in that case, read the data and return the dataframe

        try:

            ha_data = pd.read_csv(ha_filename)
            p_data = pd.read_csv(p_filename)

            self.ha_closest = ha_data
            self.p_closest = p_data
            print("Warning: you are using files that where already there")
            return [self.ha_closest, self.p_closest]
        except:
            print("There is not files containing min_dist, well compute them")


        print("##### Computing tables after data_existe  ########")


        ######## Define lipids and its ehavy atoms
        ## Example of lipid dicts
        # lipids_dict = {"DODMA": ["CN1", "CN2", "N1"]}

        #######
        keys = list(self.lipids_dict.keys())
        if ("DSPE" in keys):
            lipids_at = {"DSPE":"((resname DSPE and (" + self.build_name(self.lipids_dict["DSPE"]) +")) or ( resname PEGM and (" +self.build_name(self.pegm_ha)+ " )))"}
            del self.lipids_dict["DSPE"]
        else:
            lipids_at = {}




        for lipid in self.lipids_dict.keys():
            lipids_at[lipid] = "(resname " + lipid + " and "+ self.build_name(self.lipids_dict[lipid]) + ")"

        keys_lip_at = list(lipids_at.keys())
        h_atoms = lipids_at[keys_lip_at[0]]
        for lipids in keys_lip_at[1:]:
            h_atoms = h_atoms + " or " + lipids_at[lipids]

        ## rna heacy atoms
        rna_ha = self.build_name(self.rna_at)

        ##### This cutoff helps to optimize the computations, this way we do not take into account all the lipids but we asure that the lipids uder consideration are within a cutoff of 10A

        cutoff = 20

        ha_atoms = []
        z_ha_dist = []
        p_atoms = []
        z_p_dist = []
        #print(ha_atoms)

        for ts in self.u.trajectory[self.start:self.final:self.step]:

            ha_temp = []
            z_temp = []
            p_temp = []
            z_p_temp = []

            for id in range(self.rna_resid[0], self.rna_resid[1] + 1):
                # For each nucleotide choose the possible lipids that are around cutoff distance
                possible_lipids = self.u.select_atoms("(byres (around " + str(cutoff) + " (resid "+str(id)+"))) and not resid 1-43")

                # Further choose the heavy atoms of the lipids
                li_ha = possible_lipids.select_atoms(h_atoms)

                # Choose the residue and its heavy atoms
                nts_ha = self.u.select_atoms("resid " + str(id) + " and "+ rna_ha)

                # Choose only the atom P of the possible lipids
                li_p =possible_lipids.select_atoms("name P")


                # Compute the center of mass of the nucleotide heavy atoms
                com_nt = nts_ha.center_of_mass()
                # Choose the P atoms of the nucleotide
                p = self.u.select_atoms("resid " + str(id) + " and name P")
                min_v = np.nan
                z_mean = np.nan
                # Only do if there is possible lipids, if not, no needed
                if len(li_ha) > 0:                                                  # Compute center of mass of the ha of the lipids by residue: got a NX3 matrix
                    com_lip = li_ha.center_of_mass(compound = 'fragments')                                                                                                                                           # Cast the matrix with the COM of the nucleotide
                    distances_mat = com_lip - com_nt                                                                                                # Get the actual distances
                    distances = np.linalg.norm(distances_mat, axis = 1)                                                                             #Get the min distance
                    min_v = np.min(distances)                                                                                                             # Get the lipids within the cutoff of 10A
                    dist_to_filter = pd.DataFrame({'z' : distances_mat[:,2], 'dist' : distances})
                    dist_to_filter = dist_to_filter[dist_to_filter['dist'] <= 10]
                    z_mean = dist_to_filter['z'].mean()
                    #print(li_ha)


                ha_temp.append(min_v) # Store the min value for the nucleotides

                z_temp.append(np.abs(z_mean)) # Store the z mean distance

        # For the first nucleotide append NAN since it does not have P

                if id == 1:

                    z_p_temp.append(np.nan)

                    p_temp.append(np.nan)
                else:
                    if len(li_p) > 0:

                # Get center of m:ass of the P atoms (Not actually needed)
                        com_p = p.center_of_mass(compound = 'fragments')

                # Get the center of mass of the P atoms by residue (Not needed either)
                        com_lip_p = li_p.center_of_mass(compound = 'fragments')

                # Get distances
                        distances_p_mat = com_lip_p - com_p
                        distances_p = np.linalg.norm(distances_p_mat, axis = 1)

                # Append the min distance
                        p_temp.append(np.min(distances_p))

                # Get the lipids within the cutoff of 10A, to get the average z distance
                        dist_to_filter_p = pd.DataFrame({'z' : distances_p_mat[:,2], 'dist' : distances_p})
                        dist_to_filter_p = dist_to_filter_p[dist_to_filter_p['dist'] <= 10]
                        z_mean_p = dist_to_filter_p['z'].mean()
                        z_p_temp.append(np.abs(z_mean_p))

                    else: # If there is not lipids within the cutoff, we append NAN
                        p_temp.append(np.nan)
                        z_p_temp.append(np.nan)

            ha_atoms.append(ha_temp)
            z_ha_dist.append(z_temp)
            p_atoms.append(p_temp)
            z_p_dist.append(z_p_temp)

        ha_atoms = pd.DataFrame(ha_atoms)
        z_ha_dist = pd.DataFrame(z_ha_dist)
        p_atoms = pd.DataFrame(p_atoms)
        z_p_dist = pd.DataFrame(z_p_dist)


        ha_closest = pd.DataFrame(ha_atoms.mean(), columns = ['Mean'])
        ha_closest = ha_closest.reset_index()
        ha_closest['index'] += 1
        ha_closest.columns = ['resid', 'ha-min-distance']
        ha_closest['z-distace'] = z_ha_dist.mean()

        print('##### closes_ha written ############')
        ha_closest.to_csv(ha_filename, index = False)



        p_closest = pd.DataFrame(p_atoms.mean(), columns = ['Mean'])
        p_closest = p_closest.reset_index()
        p_closest['index'] += 1
        p_closest.columns = ['resid', 'p-min-distance']
        p_closest['z-p-distace'] = z_p_dist.mean()

        print('##### closes_p written ############')

        p_closest.to_csv(p_filename, index=False)

        #print(h_atoms)
        self.ha_closest = ha_closest
        self.p_closest = p_closest
        return [ha_closest, p_closest]


    def rna_positions(self, ha_file = "ha_com.dat", p_file = "p_pos.dat", start = None, final = None, step = None, use_default = True):

        if start == None:
            start = self.start
        if final == None:
            final = self.final
        if step == None:
            step = self.step
        if use_default:
            try:
                ha_com = pd.read_csv(ha_file, index_col = 0)
                p_pos = pd.read_csv(p_file, index_col = 0)
                self.ha_com = ha_com
                self.p_pos = p_pos
                return [ha_com, p_pos]
            except:
                print("Files are not yet present in the folder, we will compute them")


        rna_sel_ha = self.u.select_atoms("resid " + str(self.rna_resid[0])+"-"+str(self.rna_resid[1]) + " and " + self.build_name(self.rna_at))
        rna_sel_p = self.u.select_atoms("resid " + str(self.rna_resid[0])+"-"+str(self.rna_resid[1]) + " and name P")



        rna_resname = rna_sel_ha.residues.resnames
        for i in range(1,44):
            rna_resname[i-1] = rna_resname[i-1][0]+str(i)
        ha_com = []
        p_pos = []
        for ts in self.u.trajectory[start:final:step]:

                # Store rna has' positions
            rna_com = rna_sel_ha.center_of_mass(compound = 'residues')
            rna_ha_df = pd.DataFrame(rna_com, columns = ['x', 'y', 'z'])
            rna_ha_df['resid'] = rna_resname
            ha_com.append(rna_ha_df)
                # Store rna P positions
            rna_p = rna_sel_p.positions
            rna_p_df = pd.DataFrame(rna_p, columns = ['x', 'y', 'z'])
            rna_p_df['resid'] = rna_resname[1:]
            p_pos.append(rna_p_df)
                #print(rna_p_df)

        ha_com = pd.concat(ha_com)
        p_pos = pd.concat(p_pos)

        ha_com.to_csv('ha_com.dat')
        p_pos.to_csv('p_pos.dat')
        self.ha_com = ha_com
        self.p_pos = p_pos

        return [ha_com, p_pos]

   #def compute_areas(self):


    # Function to get the area
    def get_area(self, df,resid):
        temp = df[df['resid'] == resid]
    #   print(temp)
        kde_plot = sns.kdeplot(data=temp, x = 'x', y = 'y', levels = 1)
        if resid == 12:
            plt.savefig('res.png')
        #print(temp)
        contour_lines = kde_plot.collections
        area = np.nan
        for contour_line in contour_lines[:1]:
            try:
                x_values = contour_line.get_paths()[0].vertices[:,0]
                y_values = contour_line.get_paths()[0].vertices[:,1]
                #print(x_values, y_values)
                vertices = pd.DataFrame({'x' : x_values, 'y' : y_values, 'code': contour_line.get_paths()[0].codes})
                area = simps(vertices['x'], vertices['y'])
            except:
                print("Problem with contour")
        plt.close()
        return area

    def compute_areas(self, ha_file = "ha_com.dat",
                            p_file = 'p_pos.dat',
                            out_ha_file = "full_ha_table.dat",
                            out_p_file = "full_p_table.dat"

                            ):

        try:
            self.ha_table = pd.read_csv(out_ha_file)
            self.p_table = pd.read_csv(out_p_file)
            print("Files full_*_table.dat already exist, we will use those to do further computations")

            return [self.ha_table, self.p_table]
        except:
            print("Files full_*_table.dat does not exist, we will compute them now...")



        try:
            ha_com = self.ha_com
            p_pos = self.p_pos
            print(ha_com, p_pos)
            print("We are using preexistent h_com, p_pos to compute areas")
        except:
            ha_com, p_pos = self.rna_positions(ha_file, p_file)
            print("no files related withh_com, p_pos to compute areas, we will computhem now ...")


        print("ha", ha_com,"##########p", self.ha_com)

        rna_resname = ha_com['resid'].unique().tolist()
        print(rna_resname)
        # Store the area computed
        areas_ha = []
        areas_p = []

        for resid in rna_resname:
            areas_ha.append(self.get_area(ha_com, resid))
            if resid != rna_resname[0]:
                areas_p.append(self.get_area(p_pos, resid))

        # Read min-distance and z distance data
        try:
            ha_data = self.ha_closest
            p_data = self.p_closest
            print("We will use preexisting closest_nts_*.dat files")
        except:
            self.min_dist()
            ha_data = self.ha_closest
            p_data = self.p_closest

        #ha_data = pd.read_csv('closest_nts_ha.dat')
        #p_data = pd.read_csv('closest_nts_p.dat')

        ha_data['area'] = areas_ha
        ha_data['resname'] = rna_resname
        ha_data  = ha_data[ha_data['ha-min-distance'] <= 10]
        ha_data = ha_data.sort_values(by = 'ha-min-distance')
        self.ha_table = ha_data
        ha_data.to_csv(out_ha_file)
        #print(ha_data)

        p_data['area'] = [np.nan] + areas_p
        p_data['resname'] = rna_resname
        p_data  = p_data[p_data['p-min-distance'] <= 10]
        p_data = p_data.sort_values(by = 'p-min-distance')
        p_data.to_csv(out_p_file)
        self.p_table = p_data

        return [self.ha_table, self.p_table]



    def nose():
        print("hs")

    @staticmethod
    def activate_plot():
        import matplotlib as mpl
        label_size=25
        tick_width=3
        tick_size=5
        ax_label_size=30
        #plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["figure.dpi"] = 400
        plt.rcParams['figure.figsize'] = (10,10)
        plt.rcParams['font.size'] = 20
        plt.rcParams['axes.titlesize'] = 20
        plt.rcParams['figure.titlesize'] = 23
        mpl.rcParams['xtick.major.size'] = tick_size
        mpl.rcParams['xtick.major.width'] = tick_width
        mpl.rcParams['ytick.major.size'] = tick_size
        mpl.rcParams['ytick.major.width'] = tick_width
        plt.rcParams['axes.labelsize'] = ax_label_size
        plt.rcParams['xtick.labelsize'] = label_size
        plt.rcParams['ytick.labelsize'] = label_size
        mpl.rcParams['xtick.direction'] = 'in'
        mpl.rcParams['ytick.direction'] = 'in'
        plt.rcParams['axes.linewidth'] = 3
        plt.rcParams['lines.linewidth'] = 2
        #mpl.rcParams['axes.prop_cycle'] = mpl.cycler('color', ['#2ca02c','#1f77b4', '#ff7f0e', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])
        plt.rcParams["figure.autolayout"] = True





############## Plot the contour of the RNA and the 6 first nucleotides ##########################
    def full_contour(self, data_df, filename, nucl):
        #print(data_df)

        plt.close()
        kde_plot_1=sns.kdeplot(data=data_df, x='x',y='y', levels = 1, color = 'black')
	    #df_1 = df.sort_values(by='Mean')
        nts_list = nucl['resname'].iloc[:7].tolist()
        self.color = ['tab:blue', 'tab:orange', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink']
        thresh_list = [0.25, 0.35, 0.5, 0.65, 0.8]
        kde_plot_1=sns.kdeplot(data=data_df, x='x',y='y', levels = thresh_list , color='black')


	    #print(nts_list[:7])
        count = 0
        for nt in nts_list[:6]:
	        temp = data_df[data_df['resid'] == nt]
	        pos = temp[['x', 'y']].mean()
	        plt.text(pos[0], pos[1], str(nt))

	        sns.kdeplot(data=temp, x='x',y='y', levels = thresh_list, fill=False, color=self.color[count])
	        sns.kdeplot(data = temp, x = 'x', y = 'y', fill = True, color=self.color[count], alpha = 0.5)
	        count += 1
        plt.xlabel('x $[\\AA]$')
        plt.xlim(self.lower_lim,self.upper_lim)
        plt.ylim(self.lower_lim,self.upper_lim)

        legend_handles = []
        for i in range(len(nts_list[:6])):
	        legend_handles.append(Patch(color = self.color[i], label = nts_list[i]))
        plt.legend(handles=legend_handles, labelspacing = 0, handletextpad=0.1)
        plt.ylabel('y $[\\AA]$')
        plt.savefig(filename+ ".png")
        plt.show()


###################################################

    def plot_p(self, x , y, co):
        fig, ax1 = plt.subplots()
        ax1.plot(x,y)
        print(x,y)
        plt.savefig(co + 'testing.png')

    def get_lims(self, data_df):
        fig, temp = plt.subplots()

        kde_plot=sns.kdeplot(data=data_df,
                                        x='x',
                                        y='y',
                                        fill=False,
                                        levels = 1,
                                        color = 'black',
                                        alpha = 0.8)

        contour_lines_1 = kde_plot.collections
        x_values = contour_lines_1[0].get_paths()[0].vertices[:,0]
        y_values = contour_lines_1[0].get_paths()[0].vertices[:,1]
        print(data_df, len(contour_lines_1))
        xlim_low  = np.min(x_values)
        xlim_top = np.max(x_values)
        ylim_low  = np.min(y_values)
        ylim_top = np.max(y_values)
        #if self.test:
        self.lower_lim = min(xlim_low, ylim_low) - 4
        self.upper_lim = max(xlim_top, ylim_top) + 4
        return [self.lower_lim, self.upper_lim]



############## Plot the contour of the RNA and the 6 first nucleotides adding a function for axis #########################
    def full_contour_ax(self, data_df, nucl, ax):
        #print(data_df)
        #plt.close()
        kde_plot_3=sns.kdeplot(data=data_df,
                                        x='x',
                                        y='y',
                                        fill=False,
                                        levels = 1,
                                        color = 'black',
                                        alpha = 0.8,
                                        ax = ax)
        print(f"Here I write the origin of the limits {self.lower_lim}, {self.upper_lim}")

	    #df_1 = df.sort_values(by='Mean')
        nts_list = nucl['resname'].iloc[:7].tolist()
        self.color = ['tab:blue', 'tab:orange', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink']
        thresh_list = [0.25, 0.35, 0.5, 0.65, 0.8]
        kde_plot_1=sns.kdeplot(data=data_df,
                                        x='x',
                                        y='y',
                                        levels = thresh_list ,
                                        color='black',
                                        alpha = 0.8,
                                        ax = ax)


	    #print(nts_list[:7])

        count = 0
        for nt in nts_list[:6]:
	        temp = data_df[data_df['resid'] == nt]
	        pos = temp[['x', 'y']].mean()
	        ax.text(pos[0], pos[1], str(nt), fontsize = 12)

	        sns.kdeplot(data=temp,
                                x='x',
                                y='y',
                                levels = thresh_list,
                                fill=False,
                                color=self.color[count],
                                alpha = 0.9,
                                ax = ax)


	        sns.kdeplot(data = temp,
                                x = 'x',
                                y = 'y',
                                fill = True,color=self.color[count],
                                alpha = 0.5,
                                ax = ax)
	        count += 1

        ax.set_xlabel('x $[\\AA]$')
        ax.set_xlim(self.lower_lim,self.upper_lim)
        ax.set_ylim(self.lower_lim,self.upper_lim)

        plt.xlim(self.lower_lim, self.upper_lim)
        plt.ylim(self.lower_lim, self.upper_lim)

        ax.set_ylabel('y $[\\AA]$')


###################################################


###### Plot the lipid backgound with green  for each axis ################


    def lipid_background(self, filename, ax):
        lipids_data = np.loadtxt(filename, delimiter = ",")
        lipids_data = pd.DataFrame(lipids_data, columns = ['x', 'y', 'z', 'resid'])
        lipids_data['resid'] = lipids_data['resid'].astype(int)
        #print(lipids_data)
        resids_u = lipids_data['resid'].unique()

        for id in resids_u:
            temp = lipids_data[lipids_data['resid']==id]
            dis_temp = temp.diff(axis = 0).abs()

        # Choose those that cross the periodic box
            val_x_50 = dis_temp[dis_temp['x']>50]
            val_y_50 = dis_temp[dis_temp['y']>50]

        # Filter lipids that cross the barrier
            if len(val_x_50) > 0 or len(val_y_50) > 0:
                boolean = True
            else:
            # If lipids has not cross the periodic box, plot it
                kde = sns.kdeplot(data = temp, x =     'x', y = 'y', fill = True, alpha = 0.5, color='green', ax = ax)
                #print("Added to ax")
        #ax.set_ylim(self.lower_lim, self.upper_lim)
        #ax.set_xlim(self.lower_lim, self.upper_lim)
        #plt.xlim(self.lower_lim, self.upper_lim)
        #plt.ylim(self.lower_lim, self.upper_lim)
        title = filename.split("/")
        title = title[-1]
        title = title.replace("cumulative", "")
        title = title.replace("top.dat", "")
        title = title.replace("top.dat", "")
        title = title.replace("bot.dat", "")
        ax.set_title(title)

##### Plot the final plor for 4 or 3 lipids passed as argument #############3



    def lipid_contour_plot_p(self, lipids_paths, out_file):

        #print(lipids_paths)
        plt.rcParams['lines.linewidth'] = 1.5

        if len(lipids_paths) == 4:
            fig, ax = plt.subplots(2,2,sharex=True, sharey = True, layout='constrained')


            i = 0
            j = 0
            for filename in lipids_paths:
                self.lipid_background(filename, ax[i][j])
                self.full_contour_ax(self.p_pos, self.p_table, ax[i][j])
                if i == 1:
                    i = -1
                    j = 1
                i += 1




        else:
            plt.rcParams['figure.figsize'] = (17,6)
            fig, ax = plt.subplots(1,3,sharex=True, sharey = True, layout='constrained')
            i = 0
            for filename in lipids_paths:
                self.lipid_background(filename, ax[i])
                self.full_contour_ax(self.p_pos, self.p_table, ax[i])
                ax[i].set_ylim(self.lower_lim,self.upper_lim)
                ax[i].set_xlim(self.lower_lim,self.upper_lim)
                i += 1






        nts_list = self.p_table['resname'].iloc[:7].tolist()
        color = self.color
        legend_handles = []
        for i in range(len(nts_list[:6])):
            legend_handles.append(Patch(color = color[i], label = nts_list[i]))
        fig.legend(handles=legend_handles,
                        labelspacing = 0, handletextpad=0.1,
                        loc = "outside lower center",
                        frameon = False,
                        ncol = len(color))


        plt.savefig("plt" + out_file)
        plt.show()
        fig.savefig(out_file)
        plt.close()
        plt.rcParams['figure.figsize'] = (10,10)







    def lipid_contour_plot_ha(self, lipids_paths, out_file):

        plt.rcParams['lines.linewidth'] = 1.5


        if len(lipids_paths) == 4:
            fig, ax = plt.subplots(2,2,sharex=True, sharey = True, layout='constrained')


            i = 0
            j = 0
            for filename in lipids_paths:
                self.lipid_background(filename, ax[i][j])
                self.full_contour_ax(self.ha_com, self.ha_table, ax[i][j])
                if i == 1:
                    i = -1
                    j = 1
                i += 1





        else:
            plt.rcParams['figure.figsize'] = (17,6)
            fig, ax = plt.subplots(1,3,sharex=True, sharey = True, layout='constrained')
            i = 0
            for filename in lipids_paths:
                self.lipid_background(filename, ax[i])
                self.full_contour_ax(self.ha_com, self.ha_table, ax[i])
                print(f"############ Here I put the lipids {self.lower_lim}, {self.upper_lim}")
                ax[i].set_ylim(self.lower_lim, self.upper_lim)
                ax[i].set_xlim(self.lower_lim, self.upper_lim)
                i += 1





        nts_list = self.ha_table['resname'].iloc[:7].tolist()
        color = self.color
        legend_handles = []
        for i in range(len(nts_list[:6])):
            legend_handles.append(Patch(color = color[i], label = nts_list[i]))
        fig.legend(handles=legend_handles,
                        labelspacing = 0, handletextpad=0.1,
                        loc = "outside lower center",
                frameon = False,
                ncol = len(color))


        plt.savefig("plt" + out_file)
        plt.show()
        fig.savefig(out_file)
        plt.rcParams['figure.figsize'] = (10,10)
        plt.close()



##################### Plot to compare the value of the areas and dist ####################


    def plot_area_dis(self,filename1, filename2, cosa):
        identity = 1
        if cosa == "area":
            identity = 2
        data = pd.read_csv(filename1)
        data_2 = pd.read_csv(filename2)
        columns = list(data.columns)
        columns = columns[2:]
        data = data[columns]
        data_2 = data_2[columns]
        data = data.loc[:9]
        data_2 = data_2.loc[:9]
        nts_u = ["A", "G", "C", "U"]

        count = 0
        count1 = 0
        fig,ax1 = plt.subplots(2,2)
        for nt in nts_u:
            temp = data[data["resname"].str.contains(nt)]
            temp1 = data_2[data_2["resname"].str.contains(nt)]
            #print('####',temp[columns[2]], "####3", temp, "####", temp[columns[2]].idxmax())
            temp = temp.sort_values(by=columns[identity])
            temp1 = temp1.sort_values(by=columns[identity])
        #print("Ordered_columns", temp, temp1)

    #max_table = pd.DataFrame(max_table, columns = ['Nucleotide', 'Distance', 'dist-resid', 'Area', 'area_resid'])


        #print(temp)
            ax1[count1][count].plot(temp['resname'], temp[columns[identity]], color= 'tab:blue', marker = 'o', label = 'touchbot')
            ax1[count1][count].plot(temp1['resname'], temp1[columns[identity]], color= 'tab:orange', marker = 'o', label = 'highest')
            ax1[count1][count].tick_params(axis='x', rotation=45)
            print(count1, count)
            if count == 1:
                count = -1
                count1 = 1
            count+=1
            plt.legend()




        print(columns, columns[identity], cosa)


        fig.text(0.02, 0.5, 'Area $\\AA^2$', va='center', rotation='vertical')
        fig.text(0.5, 0.02, 'Nucleotides', ha='center')
        name = filename2.split("/")[-1]
        plt.savefig("plot_maxplot_"+columns[identity]+name+".png")
        plt.show()
        plt.close()

















    @staticmethod
    def radial_tilt(file_closest = "full_ha_table.dat",layer = "top",ha_com = "ha_com.dat",  file_surf_vect = "pd_DSPC_240.dat", filename = "test.png", ax_p = None):
        table_closest = pd.read_csv(file_closest)
        data_tilt = pd.read_csv(file_surf_vect)
        ha_com = pd.read_csv(ha_com)
        columns_ha = list(ha_com.columns)
        ha_com = ha_com[columns_ha[1:]]
        ha_com = ha_com.groupby("resid").mean()

        table_closest = table_closest.iloc[:10]
        axis = ["x", "y", "z"]
        plt.close()
        # Get the vectors normalized
        data_tilt["norm"] = np.linalg.norm(data_tilt[axis].values, axis = 1)
        for ax in axis:
            data_tilt[ax] = data_tilt[ax]/data_tilt["norm"]
        if layer == "bot":
            data_tilt["z"] = -data_tilt["z"]
        # Get the z axis
        data_tilt["angle"] = np.arccos(data_tilt["z"])

        print(file_surf_vect,"data tilt:",data_tilt)
        print(table_closest)
        print(ha_com)
        temp_pos0 = data_tilt[["x_0", "y_0"]].values


        df_data_nts = pd.DataFrame()
        for index, row in table_closest.iterrows():
            ha_pos = ha_com.loc[row["resname"]].to_numpy()[:2]
            data_tilt["dist_temp"] = np.linalg.norm(temp_pos0 - ha_pos, axis = 1)
            threshold = np.sqrt(row["area"]/np.pi)

            temp_df = data_tilt[data_tilt["dist_temp"] < threshold]
            print(temp_df)
            temp_df = temp_df.sort_values(by = "dist_temp")
            hist, bins = np.histogram(temp_df["angle"], bins = 180, range=(0,2*np.pi))
            print("largo del ######",len(temp_df))
            if len(temp_df) > 0:
                hist = hist/len(temp_df)
            #hist_w, bins = np.histogram(temp_df["dist_temp"],weights = temp_df["angle"], density = True, bins = 10, range=(0,threshold))
            #means = hist_w
            width = bins[1]-bins[0]
            bins = bins + width/2
            bins = bins[:-1]
            df_data_nts["bins"] = bins
            df_data_nts[f"{row['resname']}"] = hist

            if ax_p != None:
                print(f"Max value bins {max(hist)}, min value {min(hist)}, width{width}")
                ax_p.bar(bins, hist,width = width, bottom = 0.,alpha = 0.5,label = row["resname"])
                print(ax_p.get_yticklabels())
            else:
                plt.plot(bins, hist, label = row["resname"])

        if ax_p != None:

            #ax_p.set_ylabel("Mean angle [Deg]")
            yticks = [0,0.2,0.4,0.6,0.8,1]
            ax_p.set_yticks(yticks)
            ax_p.set_ylim(bottom = -0.2,top = None)
            ax_p.legend(loc = "lower center", ncols = 2)
        else:
            plt.ylabel("Mean angle [Deg]")
            plt.xlabel("Radius [$\\AA$]")
            plt.legend(loc = "lower center", ncols = 2)

        plt.savefig(filename)

        return df_data_nts



















    #################### Order parameters ################################


    @staticmethod
    def get_vectors(lista):
        angles = [] # Store the angles
        for i in (range(len(lista)-1)): # Accounts for variable number of list (Change if the carbon has or not double bonds)
            vectores = lista[i+1].positions - lista[0].positions # Hidrogen - Carbons; output of shape (n_lipids, 3)
            #print(vectores)
            angles.append(vectores)


        angles = np.concatenate(angles, axis = 0)


        costheta = angles[:,2]**2/np.linalg.norm(angles, axis = 1)**2 # Compute the costheta^2
        #print(costheta)
        order = np.mean(costheta) # Expect dim = 1, order parameter computer or the corresponding carbon
        #print(order)
        return order





    ######### Computes order parameters para sn1 #############
    def order_sn1(self, sel, lipid, n_chain):
        # Loop over all the carbons in the chains
        chains = []
        for i in range(n_chain):
	    # Select the carbon C3, H*X, H*Y, H*Z, atoms for all the lipids presents in sel1
            selections = ['name C3' + str(i+2),
                     "name H" + str(i+2) + "X and not name HX",
                     "name H" + str(i+2) + "Y and not name HY",
                      "name H" + str(i+2) + "Z and not name HZ"]
            #print(selections)
	    # Creates a list to store carbons and hydrogens
            lista = []
	    # Loop over the carbons selections and hydrogens selectifor selection in selections:
            for selection in selections:
                atoms = sel.select_atoms(selection) # Here the selection is done
                #if (i == n_chain-1) and (count == 0):
                #   print('Selection:', atoms.n_atoms, n_chain,i,selections)

                if atoms.n_atoms != 0:
                    lista.append(atoms) # The atoms are appended in a list, [[C3i], [HiX], [HiY], [HiZ]]
            #print([atoms.n_atoms for atoms in lista], [selection for selection in selections])
            chains.append(self.get_vectors(lista))
        chains = np.array(chains)
        #print("an1",chains.shape, n_chain, [at.n_atoms for at in lista])
        mean_order = np.mean(chains)
        return chains






    def order_sn2(self, sel, lipid, n_chain):
        # Define list to store the chain cos^2(theta)
        chains = []

        # Loop over carbons
        for i in range(n_chain):
            # Define selections for H and C in the chain
            selections = [
                            f"name C2{i+2}",
                            f"name H{i+2}R and not name HR",
                            f"name H{i+2}S and not name HS",
                            f"name H{i+2}T and not name HT"
                        ]
            if lipid == "POPE" or lipid == "POPS":
                if selections[0] == "name C29":
                    selections[1] = "name H91"
                if selections[0] == "name C210":
                    selections[1] = "name H101"
            # Define a list to store atoms
            lista = []

            for selection in selections:
                atoms = sel.select_atoms(selection)


                if atoms.n_atoms != 0:
                    lista.append(atoms)
            # Call get_individual that computes the cos^2(theta) for each carbon.
            angles = self.get_vectors(lista)
            chains.append(angles)
            #print("sn2", angles, n_chain, [at.n_atoms for at in lista])


        chains = np.array(chains) # Expect array of dim (n_chain)
        return chains




"""
