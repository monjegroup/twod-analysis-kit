import MDAnalysis as mda
import numpy as np
import pandas as pd
from twodanalysis import MembProp
from twodanalysis.analysis import OrderParameters


class Cumulative2D(MembProp):
    def __init__(self,
                 universe,
                lipid_list = None,
                verbose = False,
                edges = None,
                ):

        super().__init__(universe,
                         lipid_list=lipid_list,
                         verbose=verbose)

        if edges is None:
            positions = self.all_head.positions[:,:2]
            self.v_min = np.min(positions)
            self.v_max = np.max(positions)
            self.edges = [self.v_min, self.v_max, self.v_min, self.v_max]
        else:
            self.edges = edges
            self.v_min = edges[0]
            self.v_max = edges[1]

        self.start = 0
        self.final = -1
        self.step = 1


        self.guess_chain_lenght()



    ############## Order parameters related code ####################
    def order_matrix(self,
                        lipid = "DOPC",
                        layer = "top",
                        nbins = 50,
                        n_chain = None,
                        edges = None,
                        start = None,
                        final = None,
                        step = 1,
                        method = "numpy",
                        ):


        """ Computes the 2D histogram for SCD
        ----------
        lipid : str
            Working lipid to compute the order parameters
        layer : str
            working layer, can be top, bot or both
        nbins : int
            number of divisions of the grid
        n_chain : int or list
            number of carbons in the first chain or in both chains, e.g., 16 or [16,18], defaults to None
        v_min : float, optional
            min value for the 2D grid, by default None
        v_max : float, optional
            min value for the 2D grid,, by default None
        all_head : AtomGroup, optional
            atoms considered to define the middle of the membrane (all p atoms used as default), by default None
        start : int, optional
            start frame, by default None
        final : int, optional
            final frame, by default None
        step : int, optional
            frames to skip, by default 1
        method : str, optional
            method to compute the 2d histogram, by default "numpy" which uses np.histogram2d for each carbon in the lipid tails.

        Returns
        -------
        ndarray(n_grid,n_grid), ndarray(4) (Still check the return)
            matrix containind the 2D SCD and the edges in the following disposition [v_min,v_max,v_min,_vmax] (Can be used to plot directly with extent)
        """

        all_head = self.all_head

        if start is None:
            start = self.start
        if final is None:
            final = self.final

        if layer == "top":
            sign = " > "
        elif layer == "bot":
            sign = " < "

        if n_chain is None:
            n_chain = self.chain_info[lipid]

        try:
            n_chain1 = n_chain[0]
            n_chain2 = n_chain[1]
        except:
            n_chain1 = n_chain
            n_chain2 = 0

        matrix = [] # this will store a matrix of the shape (2+n_chain,
        for _ in self.u.trajectory[start:final:step]:
            z = all_head.positions[:,2]
            z_mean = z.mean() # get middle of the membrane

            #Pick atoms in the layer
            if layer == "both":
                layer_at = self.memb.select_atoms(f"byres ((resname {lipid} and name {self.working_lip[lipid]['head']}))")
            else:
                layer_at = self.memb.select_atoms(f"byres ((resname {lipid} and name {self.working_lip[lipid]['head']}) and prop z {sign} {z_mean})")
            #print("Info:", all_head.n_atoms, z_mean, layer.n_atoms)

            only_p = layer_at.select_atoms(f"name {self.working_lip[lipid]['head']}")
            positions = only_p.positions[:,:2]
            angles_sn1 = OrderParameters.individual_order_sn1(layer_at, lipid, n_chain1)
            angles_sn1 = angles_sn1.T

            #print(angles_sn1.T.shape, positions.shape)
            #print(angles_sn1.shape, positions.shape)
            to_write = np.concatenate([positions, angles_sn1], axis = 1)
            if n_chain2 != 0:
                angles_sn2 = OrderParameters.individual_order_sn2(layer_at, lipid, n_chain2)
                angles_sn2 = angles_sn2.T
                to_write = np.concatenate([to_write, angles_sn2], axis = 1)

            matrix.append(to_write) # Expect dim (n_lipids, 2+n_chain1+n_chain2)
            #print("Frame:",to_write.shape)



        #matrix = np.array(matrix) # Expect dim (frames, n_lipids, 2+n_chain1+n_chain2)
        matrix = np.concatenate(matrix, axis = 0) # Expect dim (n_lipids*frames, 2+n_chain1+n_chain2)
        print(matrix.shape)
        v_min = self.v_min
        v_max = self.v_max
        if edges is None:
            edges = [v_min,v_max,v_min,v_max]

        if method == "numpy":
            H, edges = self.numpyhistogram2D(matrix[:,:2], matrix[:,2:], n_chain, bins = nbins, edges = edges)
        else:
            H, edges = self.histogram2D(matrix[:,:2], matrix[:,2:], n_chain, bins = nbins, edges = edges)

        H = np.rot90(H)
        H[H==0] = np.nan

        return H, edges



    # Computes teh histogram of the average order parameters in each bin
    def histogram2D(self,sample1, weights, n_chain, bins = 10, edges = None):
        """ Computes the 2D histogram of 2D data with various values taking an average of them

        Parameters
        ----------
        sample1 : np.array(n,2)
            2D data information
        weights : np.array(n,m)
            m values can be attached to the data (usually lenght of the tails)
        n_chain : int or list
            Number of carbons in each chain, e.g, 16 or [16,16] for both chains
        bins : int, optional
            Number of bins to split the space, by default 10
        edges : list(float)
            Edges for the 2D grid

        Returns
        -------
        np.array, list
            matrix containining the averaged 2D histogram, edges corresponding to te matrix
        """
        if edges is None:
            limits = [[self.v_min,self.v_max], [self.v_min, self.v_max]]
        else:
            limits = [[edges[0],edges[1]], [edges[2], edges[3]]]
        #print(v_min, v_max)
        nbin = np.empty(2,np.intp)
        edges = 2*[None]

        for i in range(2):
            edges[i] = np.linspace(limits[i][0], limits[i][1], bins +1)
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
        edges = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]]

        return hist, edges

    # Computes teh histogram of the average order parameters in each bin
    def numpyhistogram2D(self,sample1, weights, n_chain, bins = 10, edges = None):
        """ Computes the 2D histogram of 2D data with various values taking an average of them

        Parameters
        ----------
        sample1 : np.array(n,2)
            2D data information
        weights : np.array(n,m)
            m values can be attached to the data (usually lenght of the tails)
        n_chain : int or list
            Number of carbons in each chain, e.g, 16 or [16,16] for both chains
        bins : int, optional
            Number of bins to split the space, by default 10
        edges : list(float)
            Edges for the grid [xmin,xmax,ymin,ymax]
        Returns
        -------
        np.array, np.array
            matrix containining the averaged 2D histogram, edges corresponding to te matrix
        """
        if edges is None:
            xmin = self.v_min
            xmax = self.v_max
            ymin = self.v_min
            ymax = self.v_max
        else:
            xmin = edges[0]
            xmax = edges[1]
            ymin = edges[2]
            ymax = edges[3]

        #print(sample1.shape, weights.shape)
        if type(n_chain) == int:
            n_chain = [n_chain]

        hist, xedges,yedges = np.histogram2d(sample1[:,0], sample1[:,1], bins = bins, range = [[xmin, xmax], [ymin, ymax]])
        hist[hist == 0] = 1

        count = 0
        mat_chain = []

        n_feat = 0
        for chain in n_chain:
            n_feat += chain

        #print("here", n_feat, weights.shape[1])
        if n_feat == weights.shape[1]:
            "Las dimensiones don correctad"

        weights = 1.5*weights-0.5

        for chain in n_chain:
            matrix = np.zeros(hist.shape)
            for i in range(chain):
                temp, xedges, yedges = np.histogram2d(sample1[:,0], sample1[:,1], weights = weights[:,count * n_chain[0]+ i], bins = bins, range = [[xmin, xmax], [ymin, ymax]])
                matrix += np.abs(temp)
            count += 1
            matrix = matrix/(chain*hist)
            matrix[matrix == 0] = np.nan
            mat_chain.append(matrix)

        matrix = 0.5*(mat_chain[0] +mat_chain[1])


        return matrix, [xedges[0],xedges[-1], yedges[0], yedges[-1]]





    @staticmethod
    def count_order(data, min_lenght, n_chain):
        """ Function used to count and average the order parameter in each grid square

        Parameters
        ----------
        data : ndarray(n,2 + len(n_chain[0]) or 2 + len(n_chain[0]) +len(n_chain[1]))
            Array that contains the order parameters data to be averaged.
        min_lenght : int
            size of the data
        n_chain : int or list
            Sets the number of carbons in the fatty acida chain

        Returns
        -------
        _np.array
            Array containing the mean SCD for each grid square
        """
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


    def all_lip_order(self, layer, nbins,
                        edges = None,
                        start = None,
                        final = None,
                        step = 1,
                        plot = False):
        """all_lip_order Find the 2D order parameters for all lipids



        Parameters
        ----------
        layer : (str)
            Layer, can be top, bot, both
        nbins : (int)
            number of bins
        edges : list(float)
            Edges for the grid in the shape [xmin,xmax,ymin,ymax]
        all_head : (mda selection, optional), optional
            heads to be considered, by default None
        start : (int, optional), optional
            start frame, by default None
        final : (int, optional), optional
            final frame, by default None
        step : (int, optional), optional
            step, by default 1
        plot : bool, optional
            plot the resulting matrix, by default False

        Returns
        -------
        ndarray(n,n), ndarray(n+1)
            matrix containing the 2d order, edges of the matrix
        """





        lipid_list = list(self.lipid_list)
        if "CHL1" in lipid_list:
            lipid_list.remove("CHL1")

        lipids = self.chain_info






        matrices = []
        for key in lipid_list:
            print(key)
            H, edges = self.order_matrix(key, layer, nbins, lipids[key],edges,
                        start = start,
                        final = final,
                        step = step)

            matrices.append(H)
        matrices = np.array(matrices)
        matrices = np.nanmean(matrices, axis = 0)

        if plot:
            plt.close()
            plt.imshow(matrices[1:-1,1:-1] ,cmap = "Spectral", extent = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]])
            plt.colorbar(cmap = "Spectral")
            plt.savefig(f"all_lip1_{layer}.png")
            plt.close()

        return matrices, edges


    ############## End of order parameters related code ############################3



    def surface(self,
                start = None,
                final = None,
                step = None,
                lipid_list = "DSPC",
                layer = 'top',
                include_charge = False,
                splay = False):
        """Code to loop over the trajectory and print [x,y,z(referenced to zmean), charge] in a file.


        Parameters
        ----------
        start : int, optional
            Start Frame, by default None
        final : int, optional
            Final frame, by default None
        step : int, optional
            Frames to skip, by default None
        lipids : str or list, optional
             Lipid to work, by default "DSPC"
        layer : str, optional
            Layer to work, by default 'top'
        filename : str, optional
            filename to write data, by default None
        include_charge : bool, optional
            Include or not charge, by default False

        Returns
        -------
        _type_
            _description_
        """


        if start is None:
            start = self.start
        if final is None:
            final = self.final
        if step is None:
            step = self.step





        if layer == "top":
            sign = " > "
        elif layer == "bot":
            sign = " < "


        ##### Select all the P atoms to find the middle of the membrane
        all_p = self.all_head

        if isinstance(lipid_list, list):
            names = f" {self.build_name(self.working_lip[lipid_list[0]]['head'])}"
            selection_string = f"(resname {lipid_list[0]} and {names}) "
            for lipid in lipid_list[1:]:
                selection_string += f" or (resname {lipid} and {self.build_name(self.working_lip[lipid]['head'])}) "

        else:
            lipid_list = [lipid_list]
            selection_string = f"(resname {lipid_list} and {self.build_name(self.working_lip[lipid]['head'])}) "

        #print(self.working_lip, selection_string)
        columns = ["x", "y", "z", "id"]


        if splay:
            self.guess_last_cs()
            carbons1 = [self.working_lip[lipid]["last_c"][0] for lipid in lipid_list]
            carbons2 = [self.working_lip[lipid]["last_c"][1] for lipid in lipid_list]
            heads = [self.working_lip[lipid]["head"] for lipid in lipid_list]
            heads = self.build_resname_atom(lipid_list, heads)
            carbons1 = self.build_resname_atom(lipid_list, carbons1)
            carbons2 = self.build_resname_atom(lipid_list, carbons2)
            columns.append("splay")

        if self.verbose:
            print("######### Running surface function ########## ")
            print(f"We will compute the surface files for a membrane with there lipids {lipid_list}")
            print(f"Currently working on: {lipid_list}")
            print(f"Layer: {layer}")



        #### Loop over trajectory to find the lipids in the requested membrane
        pos_data = []
        for _ in self.u.trajectory[start:final:step]:
            positions = all_p.positions[:,2]
            mean_z = positions.mean()

            # Selects the lipid head and of the working lipid
            if layer == "both":
                final_selection_string = selection_string
            else:
                #print(selection_string)
                final_selection_string = f"({selection_string}) and prop z {sign} {str(mean_z)}"

            # Select the atoms in the head
            if splay:
                final_selection_byres = f"byres {final_selection_string}"
                lipid_ats = self.u.select_atoms(final_selection_byres)
                head_p = lipid_ats.select_atoms(heads)
                c1 = lipid_ats.select_atoms(carbons1)
                c2 = lipid_ats.select_atoms(carbons2)
                v1 = c1.positions - head_p.positions
                v2 = c2.positions - head_p.positions
                #print(heads,carbons2,carbons1, head_p.n_atoms, c1.n_atoms)
                costheta = np.sum(v1 * v2, axis=1)/(np.linalg.norm(v1, axis = 1)* np.linalg.norm(v2, axis = 1))# Compute the cos of splay angle, must bet lenmght nlipids
                costheta = np.rad2deg(costheta[:,np.newaxis])
                atoms = lipid_ats.select_atoms(selection_string)
            else:
                atoms = self.u.select_atoms(final_selection_string)

            ### Get positions # Maybe have to check what happens with masses 0
            atom_pos = atoms.center_of_mass(compound="residues")




            ### Get resids
            atom_resid = atoms.residues.resids
            #print(atom_resid.shape)
            #print(f"costheta {costheta.shape}")
            atom_resid = atom_resid[:,np.newaxis]
            #print(atom_resid.shape)

            atom_pos = np.concatenate((atom_pos, atom_resid), axis = 1)
            if splay:
                atom_pos = np.concatenate((atom_pos, costheta), axis = 1)
            atom_pos[:,2] = np.abs(atom_pos[:,2]-mean_z)





            pos_data.append(atom_pos)



        pos_data = np.concatenate(pos_data, axis = 0)
        df_data = pd.DataFrame(pos_data, columns = columns)
        df_data["id"] = df_data["id"].astype(int)
        #if include_charge:
        #    df_data["charge"] = self.charge_li[lipid]
        #    df_data.to_csv(f"pd_{filename}", index = False)
        #df_data.to_csv(f"pd_{filename}", index = False)

        #print(df_data)

        return df_data   # Maybe have to change, it does not make sense to return thi



    @staticmethod
    def build_resname_atom(resnames, atomsnames):
        resnames = list(resnames)
        string = f"( (resname {resnames[0]}  and name {atomsnames[0]} ) "

        for i, resname in enumerate(resnames[1:]):
            string += f" or (resname {resname}  and name {atomsnames[ i + 1]}) "

        string +=  " ) "
        return string


    #def individual_splay(self, atom_group, lipids):



    def height_matrix(self, lipid_list, layer, edges = None, start = None, final = None, step = None, nbins = 50):
        """ Code to divide the space in a 2D grid and compute the height referenced to zmean

        Parameters
        ----------
        lipids : list(str)
            Lipids to include in the height analysis
        layer : str
            Working layer for thickness
        edges : list
            Edges for the grid
        start : int, optional
            Frame to start analysis, by default None
        final : int, optional
            Final frame for the analysis, by default None
        step : int, optional
            Steps to skip, by default None
        nbins : int, optional
            Number of bins to divide the grid space, by default 50
        clean : bool, optional
            Decide if rerun and overwrite surface generated files, by default True

        Returns
        -------
        ndarray(nbins,nbins)
            Retun a matrix with the height information
        """


        if start is None:
            start = self.start
        if final is None:
            final = self.final
        if step is None:
            step = self.step


        print(f"Computing matrix for {layer} in frames {start}-{final}")
        data = []

        filename = f"{lipid_list[0]}_{layer}_{start}-{final}.dat"

        data = self.surface(lipid_list= lipid_list, layer = layer, include_charge = True, start = start, final = final)



        if edges is not None:
            xmin = edges[0]
            xmax = edges[1]
            ymin = edges[2]
            ymax = edges[3]
        else:
            xmin = self.v_min
            xmax = self.v_max
            ymin = self.v_min
            ymax = self.v_max


        H_height, x_edges, y_edges = np.histogram2d(x = data["x"], y = data["y"], weights = data["z"], bins = nbins, range = [[xmin,xmax], [ymin,ymax]])
        H_count, x_edges, y_edges = np.histogram2d(x = data["x"], y = data["y"], bins = nbins, range = [[xmin, xmax], [ymin,ymax]])

        H_count[H_count == 0] = 1.

        H_avg = H_height/H_count

        H_avg[H_avg == 0] =  np.nan

        H_avg = np.rot90(H_avg)

        np.savetxt(f'Height_{layer}_{start}_{final}.dat', H_avg, fmt = '%.2f')
        np.savetxt(f"edges_{layer}_{start}_{final}.dat", x_edges, fmt = "%.2f")
        return H_avg, [x_edges[0],x_edges[-1],y_edges[0], y_edges[-1]]


    def splay_matrix(self, lipid_list, layer, edges = None, start = None, final = None, step = None, nbins = 50):
        """ Code to divide the space in a 2D grid and compute the height referenced to zmean

        Parameters
        ----------
        lipids : list(str)
            Lipids to include in the height analysis
        layer : str
            Working layer for thickness
        edges : list
            Edges for the grid
        start : int, optional
            Frame to start analysis, by default None
        final : int, optional
            Final frame for the analysis, by default None
        step : int, optional
            Steps to skip, by default None
        nbins : int, optional
            Number of bins to divide the grid space, by default 50
        clean : bool, optional
            Decide if rerun and overwrite surface generated files, by default True

        Returns
        -------
        ndarray(nbins,nbins)
            Retun a matrix with the height information
        """


        if start is None:
            start = self.start
        if final is None:
            final = self.final
        if step is None:
            step = self.step




        print(f"Computing matrix for {layer} in frames {start}-{final}")
        data = []

        data = self.surface(lipid_list=lipid_list, layer = layer, include_charge = True, start = start, final = final, splay=True)


        #print(data)

        if edges is not None:
            xmin = edges[0]
            xmax = edges[1]
            ymin = edges[2]
            ymax = edges[3]
        else:
            xmin = self.v_min
            xmax = self.v_max
            ymin = self.v_min
            ymax = self.v_max


        H_splay, x_edges, y_edges = np.histogram2d(x = data["x"], y = data["y"], weights = data["splay"], bins = nbins, range = [[xmin,xmax], [ymin,ymax]])
        H_count, x_edges, y_edges = np.histogram2d(x = data["x"], y = data["y"], bins = nbins, range = [[xmin, xmax], [ymin,ymax]])

        H_count[H_count == 0] = 1.

        H_avg = H_splay/H_count

        H_avg[H_avg == 0] =  np.nan

        H_avg = np.rot90(H_avg)

        return H_avg, [x_edges[0],x_edges[-1],y_edges[0], y_edges[-1]]



    def thickness(self, nbins, edges = None,lipid_list = None, start = 0, final=-1, step = 1):
        """Find the thichness mapped in a 2d grid

        Parameters
        ----------
        nbins : int
            number of bins for thickness
        edges : list
            Edges for the grid
        lipids : list(str)
            Lipids to be considered in the thickness computation
        start : int, optional
            Start frame, by default 0
        final : int, optional
            Final frame, by default -1
        step : int, optional
            Step frame, by default 1

        Returns
        -------
        np.array, np.array
            Matrix with the thickness, edeges for the matrix
        """
        if lipid_list is None:
            lipid_list = list(self.lipid_list)
            if "CHL1" in lipids:
                lipid_list.remove("CHL1")

        matrix_up, edges = self.height_matrix(lipid_list,
                        "top",
                        edges = edges,
                        start = start,
                        final = final,
                        step = step,
                        nbins = nbins)
        matrix_bot, edges = self.height_matrix(lipid_list,
                        "bot",
                        edges = edges,
                        start = start,
                        final = final,
                        step = step,
                        nbins = nbins)

        mat_thickness = matrix_bot + matrix_up

        mat_thickness[mat_thickness == 0] = np.nan
        #print(mat_thickness,mat_thickness.shape,matrix_bot.shape,[edges[0], edges[-1], edges[0], edges[-1]])

        return mat_thickness, edges
