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
    def order_histogram(self,
                        lipid = "DOPC",
                        layer = "top",
                        nbins = 50,
                        n_chain = [17,17],
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
            number of carbons in the first chain or in both chains, e.g., 16 or [16,18]
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
                        all_head = None,
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
            H, edges = self.order_histogram(key, layer, nbins, lipids[key],edges,
                        all_head = all_head,
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
