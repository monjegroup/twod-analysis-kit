"""
Memb2D
=============
Class created mainly to analyze lipid membranes in different ways

Classes
-------

.. autoclass:: Memb2D
    :members:
    :undoc-members:
    :show-inheritance:


"""


import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt







class Memb2D:

    def __init__(
                self,
                obj,
                lipid_list = None,
                guess_chain_l = True,
                chain_info = None,
                v_min = None,
                v_max = None,
                add_radii = False,
                verbose = False,
            ):
        """__init__ _summary_

        _extended_summary_

        Parameters
        ----------
        top : _type_
            _description_
        traj : _type_
            _description_
        lipid_list : _type_, optional
            _description_, by default None
        tpr : _type_, optional
            _description_, by default None
        info : bool, optional
            _description_, by default False
        guess_chain_l : bool, optional
            _description_, by default True
        chain_info : _type_, optional
            _description_, by default None
        v_min : _type_, optional
            _description_, by default None
        v_max : _type_, optional
            _description_, by default None
        add_radii : bool, optional
            _description_, by default False
        verbose : bool, optional
            _description_, by default False
        """



        # Read trajectory depending if tpr is provided or not

        if isinstance(obj, mda.Universe):
            self.u = obj
        elif isinstance(obj,mda.core.groups.AtomGroup):
            self.u = obj.u
        else:
            raise TypeError("Input must be an MDAnalysis Universe or AtomGroup")




        # Select elements in the membrane (in principle only lipids)
        if not lipid_list: # Select only elements of the membrane
            self.memb = self.u.select_atoms("all and not protein and not (resname URA or resname GUA or resname ADE or resname CYT or resname THY)")
            self.lipid_list = set(self.memb.residues.resnames)
        else:
            self.memb = self.u.select_atoms(f"{self.build_resname(list(lipid_list))}")


        # Set ercentage for periodicity

        self.periodicity = 0.1


        # Set radius sizes of different elements
        self.radii_dict = 0
        if add_radii:
            self.radii_dict = {"H": 0.7,
                            "N": 1.85,
                            "C": 2.06,
                            "P": 2.15,
                            "O": 1.65,
                            }


            # Add radii as a topology attribute for Mdanalysis
            string_array = self.memb.elements
            radii_array = np.array([self.radii_dict[element] for element in string_array])
            self.u.add_TopologyAttr("radii")
            self.memb.radii = radii_array

            # May use to build a different way to select polar atoms
            #polar_motif = "N HN1 HN2 HN3 C12 H12A C13 O13A O13B C11 H11A H11B"
            #polar_PS = "N HN1 HN2 HN3 C12 H12A C13 O13A O13B C11 H11A H11B"
            #polar_PI = "C12 H2 O2 HO2 C13 H3 O3 HO3 C14 H4 O4 HO4 C15 H5 O5 HO5 C16 H6 O6 HO6 C11 H1"
            #polar_PA = "H12 "
            #polar_PC = "N C12 C13 C14 C15 H12A H12B H13A H13B H13C H14A H14B H14C H15A H15B H15C C11 H11A H11B"
            #polar_PE = "N HN1 HN2 HN3 C12 H12A H12B C11 H11A H11B"
            #polar_CHL = "O3 H3'"


            #polar_chains = [polar_motif, polar_PS, polar_PI, polar_PA, polar_PC, polar_PE]
            #polar_atoms = [chain.split() for chain in polar_atoms]
            #dspc = self.memb.select_atoms("(resname DSPC and not (name C3* or name H*X or name H*Y or name C2* or name H*R or name H*S)) or (resname DSPC and(name C3 or name HX or name HY or name C2 or name HR or name HS))")
            #print(set(dspc.atoms.names))








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
                                "POPI15" : {"head" :"P", "charge" : 1.3},
                                "POPI24" : {"head" :"P", "charge" : 1.3},
                            } #List of known lipids and lipids head people usually use to work

        self.chain_info = chain_info


        if guess_chain_l: # Guess the chain lenght of lipids. Chain sn2 start with C2 and chain sn1 start with C3
            self.chain_info = {}
            self.non_polar_dict = {}
            self.first_lipids = {}
            self.non_polar_visualize = {}
            for lipid in self.lipid_list:
                first_lipid = self.memb.select_atoms(f"resname {lipid}").resids[0]
                actual_sn1 = self.memb.select_atoms(f"resid {first_lipid} and name C3*")
                actual_sn2 = self.memb.select_atoms(f"resid {first_lipid} and name C2*")
                actual_sn1 = actual_sn1.names
                actual_sn2 = actual_sn2.names
                self.chain_info[lipid] = [len(actual_sn1) - 2, len(actual_sn2) - 2]
                self.first_lipids[lipid] = first_lipid

                if lipid == "CHL1":
                    non_polar = self.memb.select_atoms(f"resid {first_lipid} and not (name O3 or name H3')")
                    all_lip = self.memb.select_atoms(f"resid {first_lipid}")

                else:
                    non_polar = self.memb.select_atoms(f"resid {first_lipid} and (name *C3* or name H*Y or name H*X or name H*Z  or name *C2* or name H*R or name H*S or name H*T) and not (name C3 or name C31 or name HY or name HX or name HZ  or name C2 or name C21 or name HR or name HS or name HT)")
                    all_lip = self.memb.select_atoms(f"resid {first_lipid}")
                self.non_polar_dict[lipid] = list(non_polar.names)

                self.non_polar_visualize[lipid] = [all_lip, non_polar]



        self.all_head = self.u.select_atoms(self.build_resname(self.lipid_list) + " and name P")
        if v_min is None and v_max is None:
            positions = self.all_head.positions[:,:2]
            self.v_min = np.min(positions)
            self.v_max = np.max(positions)
        self.start = 0
        self.final = 100
        self.step = 1
        self.verbose = verbose

        if verbose:
            print(f"This system contains the following lipids : {self.lipid_list}\n\n")
            print(f"The chain lenght is : \n{self.print_dict(self.chain_info)}\n")
            print(f"We will use the following heads and charges for the following lipids. If the lipid is not here we will use P as head as default \n{self.print_dict(self.working_lip)}\n")
            print("Note: To compute the middle of the membrane we use only P heads\n\n")
            print(f"The default start frame is {self.start}, final {self.final}, step {self.step}\n\n")

    # Method to print dictionaries
    @staticmethod
    def print_dict(dictio):
        string = ""
        for key in dictio.keys():
            string += f"{key} : {dictio[key]}\n"
        return string


    def visualize_polarity(self, lipids = "all"):
        """visualize_polarity

        This function is used to visualize what atoms are considered in polarity

        Parameters
        ----------
        lipids : str, optional
            Lipids to show polarity, by default "all"
        """
        aspect_ratio = [1, 1, 1]
        # Get lipids to work
        if lipids == "all":
            lipids = self.lipid_list
        else:
            if isinstance(lipids, str):
                lipids = [lipids]
        #Guess bonds if needed
        try:
            self.u.bonds
        except:
            bonds = mda.topology.guessers.guess_bonds(self.u.atoms, self.u.atoms.positions)
            self.u.add_TopologyAttr("bonds", bonds)
            #print(mda.topology.guessers.guess_bonds(first_lipid, first_lipid.positions))

        # Make a 3D plot of the individual lipids and shows it
        if not self.non_polar_dict:
            print("Non polar atoms are not yet established, please set self.non_polar_dict first")
            return
        else:
            nplots = len(lipids)
            fig = plt.figure(figsize=(5 * nplots, 5))


            count = 0
            for lipid in lipids:
                first_lipid = self.memb.select_atoms(f"resid {self.first_lipids[lipid]} and resname {lipid}")
                lipid_pos = first_lipid.positions
                lipid_ats = first_lipid.names
                polarity = ["blue" if name in self.non_polar_dict[lipid] else "red" for name in lipid_ats]
                #print(polarity, lipid_ats)
                ax = fig.add_subplot(1,nplots, count+1, projection = "3d")
                ax.scatter(*lipid_pos.T, s = 10, c = polarity)
                dists = []
                mids = []
                for i in range(3):
                    value = np.max(lipid_pos[:,i]) - np.min(lipid_pos[:,i])
                    dists.append(value)
                    mids.append(np.min(lipid_pos[:,i] + value/2))
                dists = max(dists)
                dists = dists + 0.01*dists
                for atom in first_lipid:
                    bonds = atom.bonded_atoms
                    for bond in bonds:
                        vector = np.array([atom.position, bond.position])
                        plt.plot(*vector.T, color = "black")
                for i in range(len(lipid_ats)):
                    if lipid_ats[i] in self.non_polar_dict[lipid]:
                        ax.text(lipid_pos[i,0], lipid_pos[i,1],  lipid_pos[i,2], lipid_ats[i], color = "black", fontsize = 6)
                    else:
                        ax.text(lipid_pos[i,0], lipid_pos[i,1], lipid_pos[i,2], lipid_ats[i], color = "black", fontsize = 6)
                ax.set_box_aspect(aspect_ratio)
                ax.set_title(lipid)
                ax.set_xlim(mids[0]-dists/2, mids[0]+dists/2)
                ax.set_ylim(mids[1]-dists/2, mids[1]+dists/2)
                ax.set_zlim(mids[2]-dists/2, mids[2]+dists/2)
                count += 1
        plt.show()






    ############## Order parameters related code ####################




    def order_histogram(self, lipid, layer, n_grid,
                        n_chain,
                        edges = None,
                        all_head = None,
                        start = None,
                        final = None,
                        step = 1,
                        method = "numpy",
                        ):
        """Method that allows for the computation of order parameters in 2D fashion

        Parameters
        ----------
        lipid : str
            Working lipid to compute the order parameters
        layer : str
            working layer, can be top, bot or both
        n_grid : int
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

        if all_head is None:
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
        for ts in self.u.trajectory[start:final:step]:
            z = all_head.positions[:,2]
            z_mean = z.mean() # get middel of the membrane

            #Pick atoms in the layer
            if layer == "both":
                layer_at = self.memb.select_atoms(f"byres ((resname {lipid} and name {self.working_lip[lipid]['head']}))")
            else:
                layer_at = self.memb.select_atoms(f"byres ((resname {lipid} and name {self.working_lip[lipid]['head']}) and prop z {sign} {z_mean})")
            #print("Info:", all_head.n_atoms, z_mean, layer.n_atoms)

            only_p = layer_at.select_atoms(f"name {self.working_lip[lipid]['head']}")
            positions = only_p.positions[:,:2]
            angles_sn1 = self.individual_order_sn1(layer_at, lipid, n_chain1)
            angles_sn1 = angles_sn1.T

            #print(angles_sn1.T.shape, positions.shape)
            #print(angles_sn1.shape, positions.shape)
            to_write = np.concatenate([positions, angles_sn1], axis = 1)
            if n_chain2 != 0:
                angles_sn2 = self.individual_order_sn2(layer_at, lipid, n_chain2)
                angles_sn2 = angles_sn2.T
                to_write = np.concatenate([to_write, angles_sn2], axis = 1)

            matrix.append(to_write) # Expect dim (n_lipids, 2+n_chain1+n_chain2)
            #print("Frame:",to_write.shape)

        #matrix = np.array(matrix) # Expect dim (frames, n_lipids, 2+n_chain1+n_chain2)
        matrix = np.concatenate(matrix, axis = 0) # Expect dim (n_lipids*frames, 2+n_chain1+n_chain2)
        v_min = self.v_min
        v_max = self.v_max

        if method == "numpy":
            H, edges = self.numpyhistogram2D(matrix[:,:2], matrix[:,2:], n_chain, bins = n_grid, edges = edges)
        else:
            H, edges = self.histogram2D(matrix[:,:2], matrix[:,2:], n_chain, bins = n_grid, edges = edges)
        plt.imshow(H,cmap="Spectral")
        plt.colorbar()
        plt.close()
        H = np.rot90(H)
        H[H==0] = np.nan

        return H, edges





    @staticmethod
    def get_individual(lista
                    ):
        r"""This function gets a list with a specific carbon (e.g. C34 or C22)
        and its respective hidrogens (e.g. H4X, H4Y). It computes the vectors
        that connect the carbons and the hydrogens and computes the :math:`cos(\theta)^2`, where :math:`\theta` is
        the angle between each vector and the z-axis. Finally, this function returns a vector with the individual (per lipid)
        :math:`\braket{cos(\theta)^2}`, where the mean is computed over the hydrogens of each carbon.

        Parameters
        ----------
        lista : list
            Vector of the shape :math:`[C*i, HiX,HiY, HiZ]`, the minimun len is 2 (when the carbon
            only have one hydrogen) and the maximun is 4 (when there is three hydrogens)
            Note: If there is N lipids, there will be N carbons :math:`C*i`, and the i represents
            the position of the carbon in the lipid tail.

        Returns
        -------
        order : array(n_lipids)
            Float with the mean of :math:`\braket{cos(\theta)^2}`

        Notes
        -----
        The average of the angle of the i-th carbon for all the lipids in the selection is computed
        as follows:

        .. math:: \braket{cos(\theta_i)^2}

        where :math:`\theta_i` is the angle between the z- axis and the vector that connects the i-th carbon and the hydrogen.



        """

        angles = [] # Store the angles for the working carbon
        for i in (range(len(lista)-1)): # Accounts for variable number of list (Change if the carbon has or not double bonds)
            vectores = lista[i+1].positions - lista[0].positions # Hidrogen - Carbons; output of shape (n_lipids, 3)
            costheta = vectores[:,2]**2/np.linalg.norm(vectores, axis = 1)**2 # Compute the costheta^2
            angles.append(costheta) # dim (n_lipids,)
        angles = np.array(angles) # dim ((1, 2 or 3),n_lipids)
        #print("angles", angles.shape)
        angles = np.mean(angles, axis = 0) # output is dim n_lipids, it means the cos^2(theta) or the Carbon passed for each lipid
        return angles


    # Get the cos^2(theta) for each carbon in the selection, for sn1
    def individual_order_sn1(self, sel, lipid, n_chain):
        """

        Code to loop over the number of carbons_summary_ in the lipid tail and get a list with the carbon and its
        hydrogens for each carbon in the lipid tail: :math:`[C3i, HiX, HiY, ...]`. This list is passed to get_vectors which
        return the averages of each i-th carbon. This code returns an array of dim n_chain with the mean :math:`\braket{cos(\theta_i)^2}`

        Parameters
        ----------
        lipid : str
            Name of the lipid to compute the order parameters
        n_chain : int
            Number of carbons in the lipid tail sn1

        Returns
        -------
        chains : ndarray
            Vector dim n_chains with the mean :math:`\braket{cos(\theta_i)^2}`

        Notes
        -----
        The return is a vector containing the value :math:`\braket{cos^2(\theta_i}`. As follows:

        .. math:: [\braket{cos^2(\theta_2}, \braket{cos^2(\theta_3}, ..., \braket{cos^2(\theta_{n_chain}}]

        The index starts at 2 because that is the carbon the lipid tail starts with.

        """

        # Define list to store the chain cos^2(theta)
        chains = []


        # Loop over carbons
        for i in range(n_chain):
            # Define selections for H and C in the chain
            #print(f"Value of the chain {i} sn1")
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

        r"""

        Code to loop over the number of carbons in the lipid tail and get a list with the carbon and its
        hydrogens for each carbon in the lipid tail: :math:`[C2i, HiX, HiY, ...]`. This list is passed to get_vectors which
        return the averages of each i-th carbon. This code returns an array of dim n_chain with the mean :math:`\braket{cos(\theta_i)^2}`

        Parameters
        ----------
        lipid : str
            Name of the lipid to compute the order parameters
        n_chain : int
            Number of carbons in the lipid tail sn1

        Returns
        -------
        chains : ndarray
            Vector dim n_chains with the mean :math:`\braket{cos(\theta_i)^2}`

        Notes
        -----
        The return is a vector containing the value :math:`\braket{cos^2(\theta_i}`. As follows:

        .. math:: [\braket{cos^2(\theta_2}, \braket{cos^2(\theta_3}, ..., \braket{cos^2(\theta_{n_chain}}]

        The index starts at 2 because that is the carbon the lipid tail starts with.

        """
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
            if lipid == "POPE" or lipid == "POPS" or lipid == "POPI15" or lipid == "POPI24":
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
            angles = self.get_individual(lista)
            chains.append(angles)
        chains = np.array(chains) # Expect array of dim (n_chain, n_lipids)
        return chains



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

    """
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

    """







    @staticmethod
    def build_resname(resnames_list):
        resnames_list = list(resnames_list)
        string = " (resname " + resnames_list[0]

        for resname in resnames_list[1:]:
            string = string + " or resname " + resname

        string = string + ") "
        return string

    def build_resname_head(self,resnames_list):
        resnames_list = list(resnames_list)
        string = f"( (resname {resnames_list[0]}  and name {self.working_lip[resnames_list[0]]['head']}) "

        for resname in resnames_list[1:]:
            string += f" or (resname {resnames_list[0]}  and name {self.working_lip[resnames_list[0]]['head']}) "

        string +=  " ) "
        return string

    @staticmethod
    def build_name(resnames_list):

        if isinstance(resnames_list, str):
            return f"(name {resnames_list})"

        string = " (name " + resnames_list[0]
        for resname in resnames_list[1:]:
            string = string + " or name " + resname

        string = string + ") "
        return string






    def surface(self,
                start = None,
                final = None,
                step = None,
                lipids = "DSPC",
                layer = 'top',
                filename = None, include_charge = False):
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


        if start == None:
            start = self.start
        if final == None:
            final = self.final
        if step == None:
            step = self.step

        lipid_list = self.lipid_list

        if self.verbose:
            print("######### Running surface function ########## ")
            print(f"We will compute the surface files for a membrane with there lipids {lipid_list}")
            print(f"Currently working on: {lipids}")
            print(f"Layer: {layer}")
            print(f"Writing under the name of {filename}")


        if layer == "top":
            sign = " > "
        elif layer == "bot":
            sign = " < "


        ##### Select all the P atoms to find the middle of the membrane
        all_p = self.all_head

        if isinstance(lipids, list):

            names = f" {self.build_name(self.working_lip[lipids[0]]['head'])}"
            selection_string = f"(resname {lipids[0]} and {names}) "
            for lipid in lipids[1:]:
                selection_string += f" or (resname {lipid} and {self.build_name(self.working_lip[lipid]['head'])}) "
            if filename == None:
                filename = f"{lipids[0]}_{layer}_{start}_{final}.dat"
        else:
            selection_string = f"(resname {lipids} and {self.build_name(self.working_lip[lipids]['head'])}) "
            if filename == None:
                filename = f"{lipids}_{layer}_{start}_{final}.dat"




        #### Loop over trajectory to find the lipids in the requested membrane
        pos_data = []
        for ts in self.u.trajectory[start:final:step]:
            positions = all_p.positions[:,2]
            mean_z = positions.mean()

            # Selects the lipid head and of the working lipid
            if layer == "both":
                final_selection_string = selection_string
            else:
                #print(selection_string)
                final_selection_string = f"({selection_string}) and prop z {sign} {str(mean_z)}"

            # Find the positions of the P atoms
            atoms = self.u.select_atoms(final_selection_string)

            ### Get positions
            atom_pos = atoms.center_of_mass(compound="residues")
            #print(atom_pos.shape)



            ### Get resids
            atom_resid = atoms.residues.resids
            atom_resid = atom_resid[:,np.newaxis]
            #print(atom_resid.shape)

            atom_pos = np.concatenate((atom_pos, atom_resid), axis = 1)
            atom_pos[:,2] = np.abs(atom_pos[:,2]-mean_z)

            pos_data.append(atom_pos)



        pos_data = np.concatenate(pos_data, axis = 0)
        df_data = pd.DataFrame(pos_data, columns = ["x", "y", "z", "id"])
        df_data["id"] = df_data["id"].astype(int)
        #if include_charge:
        #    df_data["charge"] = self.charge_li[lipid]
        #    df_data.to_csv(f"pd_{filename}", index = False)
        df_data.to_csv(f"pd_{filename}", index = False)

        return df_data   # Maybe have to change, it does not make sense to return thi


    def height_matrix(self, lipids, layer, edges = None, start = None, final = None, step = None, nbins = 50, clean = True):
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

        filename = f"{lipids[0]}_{layer}_{start}-{final}.dat"
        if not clean:
            try:
                df_data = pd.read_csv(f"pd_{filename}")
            except:
                df_data = self.surface(lipids = lipids, layer = layer, filename = filename, include_charge = True, start = start, final = final)
        else:
            df_data = self.surface(lipids = lipids, layer = layer, filename = filename, include_charge = True, start = start, final = final)


        data = df_data
        #print(data)

        if edges is not None:
            xmin = edges[0]
            xmax = edges[1]
            ymin = edges[2]
            ymax = edges[3]
        else:
            xmin = data["x"].min()
            xmax = data["x"].max()
            ymin = data["y"].min()
            ymax = data["y"].max()


        H_height, x_edges, y_edges = np.histogram2d(x = data["x"], y = data["y"], weights = data["z"], bins = nbins, range = [[xmin,xmax], [ymin,ymax]])
        H_count, x_edges, y_edges = np.histogram2d(x = data["x"], y = data["y"], bins = nbins, range = [[xmin, xmax], [ymin,ymax]])

        H_count[H_count == 0] = 1.

        H_avg = H_height/H_count

        H_avg[H_avg == 0] =  np.nan

        H_avg = np.rot90(H_avg)

        np.savetxt(f'Height_{layer}_{start}_{final}.dat', H_avg, fmt = '%.2f')
        np.savetxt(f"edges_{layer}_{start}_{final}.dat", x_edges, fmt = "%.2f")
        return H_avg, [x_edges[0],x_edges[-1],y_edges[0], y_edges[-1]]



    def thickness(self, nbins, edges = None,lipids = None, start = 0, final=-1, step = 1):
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
        if lipids is None:
            lipids = list(self.lipid_list)
            lipids.remove("CHL1")

        matrix_up, edges = self.height_matrix(lipids,
                        "top",
                        edges = edges,
                        start = start,
                        final = final,
                        step = step,
                        nbins = nbins)
        matrix_bot, edges = self.height_matrix(lipids,
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



    def guess_minmax_space(self, all = False):
        """Check the minimun and max position in x,y

        Parameters
        ----------
        all : bool,
            If True return the max and min in each axis(x,y), else returns inly a minimun and maximun, by default False

        Returns
        -------
        float, float
            minimun and max position in x,y
        """
        positions = self.memb.positions[:,:2]
        if all:
            xmin = np.min(positions[:,0])
            xmax = np.max(positions[:,0])
            ymin = np.min(positions[:,1])
            ymax = np.max(positions[:,1])
            return xmin, xmax, ymin, ymax
        vmin = np.min(positions)
        vmax = np.max(positions)

        return vmin,vmax


    ################## Pcaking defects related code ###############





    def packing_defects(self,
                        layer = 'top',
                        nbins = 180,
                        height = False,
                        periodic = False,
                        edges = None,
                        area = True,
                        count = True,
                        verbose = True,
                        ):
        """Compute packing defects based on packmem: https://doi.org/10.1016/j.bpj.2018.06.025

        Parameters
        ----------
        layer : str, optional
            working layer (top/bot). Defaults to 'top', by default 'top'
        nbins : int, optional
            Number of bins of the xy grid, by default 180
        height : bool, optional
            Store height matrix (To study deepness of th packing defects), by default False
        periodic : bool, optional
            Defines if using periodicity or not. When True, takes into acount periodicity and returns a 2D grid with of the size of the periodic box, by default False
        vmin : float, optional
            Store the min value for x and y
        vmax : float, optional
            Store the max value for x and y

        Returns
        -------
        ndarray
            If height == Flase: matrix with packing defects
        ndarray, ndarray
            If height === True: matrix with packing defects, amtrix with height information
        """

        if edges is not None:
            xmin = edges[0]
            xmax = edges[1]
            ymin = edges[2]
            ymax = edges[3]
            dx = xmax-xmin
            dy = ymax-ymin
            if int(dx) != int(dy):
                print(dx,dy)
                print("Distances in x and y must be equal")
                return
        else:
            xmin = self.v_min
            ymin = self.v_min
            xmax = self.v_max
            ymax = self.v_max




        grid_size = abs(xmin - xmax)/nbins
        n_aug = int(4/grid_size)


        # Extend the grid 5 Amstrong to azure all the packing defects are correctly taken

        xmin_ex = xmin - n_aug * grid_size
        xmax_ex = xmax + n_aug * grid_size
        ymin_ex = ymin - n_aug * grid_size
        ymax_ex = ymax + n_aug * grid_size
        #vmin_ex = vmin - n_aug * grid_size
        #vmax_ex = vmax + n_aug * grid_size
        nbins_aug = nbins + 2 * n_aug



        edges = [xmin,xmax,ymin,ymax]






        lipid_list = list(self.lipid_list)


        if verbose:
            print((f'######### Running packing defects function ########## \
                  \nWe will compute packing defects for a membrane with lipids {lipid_list}'),
                  end = "\r")





        if layer == "top":
            sign = " > "
        elif layer == "bot":
            sign = " < "


        ##### Select all the P atoms to find the middle of the membrane
        all_p = self.all_head




        pos_data = []

        #### Define radious to be used for the different atoms
        mat_radii_dict = {}
        for atom in self.radii_dict.keys():
            mat_radii_dict[atom] = self.create_circle_array(grid_size, self.radii_dict[atom])


        # Create matrix to be filled
        matrix = np.zeros((nbins_aug+2, nbins_aug+2))
        if height:
            matrix_height = np.zeros((nbins_aug, nbins_aug))
        positions = all_p.positions[:,2]
        mean_z = positions.mean()
        dims = self.u.trajectory.ts.dimensions[:2]
        for lipid in self.lipid_list:


            selection_string = f"byres ((resname {lipid} and name {self.working_lip[lipid]['head']}) and prop z {sign} {mean_z})"
            #print(selection_string)
            layer_at = self.memb.select_atoms(selection_string)
            pos_ats = layer_at.positions
            elements = layer_at.elements
            names = layer_at.names

            if periodic:
                #temp = pos_ats.copy()
                pos_ats, others = self.extend_data(pos_ats, dims, self.periodicity, [elements, names])
                elements = others[0]
                names = others[1]

            if not height:
                indexes = self.get_indexes(pos_ats[:,:2], nbins_aug, edges=[xmin_ex,xmax_ex,ymin_ex,ymax_ex])
            else:
                indexes, matrix_temp = self.get_indexes(pos_ats, nbins_aug, edges=[xmin_ex,xmax_ex,ymin_ex,ymax_ex], matrix_height = True)
                matrix_height = np.maximum(matrix_height.copy(), matrix_temp.copy())



            matrix = self.add_defects(matrix, indexes,elements, names, lipid, mat_radii_dict)
        #print(f"Shapeeeeeee::::: {matrix.shape}, {matrix_height.shape}, dims {dims}")
        core1 = 2*(slice(n_aug+1,-(n_aug+1)),)
        core = 2*(slice(n_aug,-n_aug),)
        matrix = matrix[core1]


        #print(f"Shapeeeeeee::::: {matrix.shape}, {matrix_height.shape}, dims {dims}")


        if periodic:
            n = round((dims[0]/grid_size))
            #print(n, dims[0]/grid_size, grid_size)
            xmax = xmin + n*grid_size
            ymax = ymin + n*grid_size
            matrix = matrix[:n, :n]

            if height:
                matrix_height = matrix_height[:n, :n]
            edges = [xmin,xmax,ymin,ymax]




        #print(f"Shapeeeeeee::::: {matrix.shape}, {matrix_height.shape}, dims {dims}")
        defects = matrix
        defects = np.where(matrix < 1, matrix, np.nan)
        #defects = np.where(defects > 0, defects, np.nan)


        return_dict = {
            "edges" : edges,
        }

        if area:
            non_nan = ~np.isnan(defects)
            count = np.sum(non_nan)
            area_v = count*grid_size*grid_size

            area_tot = (defects.shape[0])**2 * grid_size * grid_size

            return_dict["area"] = {"defects" : area_v,
                               "total": area_tot}

        if count:
            from scipy.ndimage import label
            binary_matrix = ~np.isnan(defects)
            structure = np.array([[1,1,1], [1,1,1], [1,1,1]])
            labeled_array, num_features = label(binary_matrix, structure = structure)
            cluster_sizes = np.bincount(labeled_array.ravel())[1:]
            return_dict["count"] = {"number":num_features, "sizes":cluster_sizes}
            return_dict["grid_size"] = grid_size

            #plt.hist(cluster_sizes, bins=40)
            #plt.show()




        #print(return_dict)
        if height:
            matrix_height = matrix_height[core]
            matrix_height[matrix_height == 0 ] = np.nan
            return defects, matrix_height, return_dict


        return defects, return_dict





    @staticmethod
    def create_circle_array(grid_size, radius_A, center=None):
        """Create a small matrix with a inner circle of the size of radius_A

        Parameters
        ----------
        grid_size : float
            define the grid size to create the optimun grid (Amstrongs)
        radius_A : float
            define the radius for the matrix (amstrongs)
        center : bool, optional
            Bool to set the center of the circle, by default None

        Returns
        -------
        ndarray
            array with a circle of ones
        """



        n = int(2*radius_A /grid_size) + 1
        if n % 2 == 0:
            n += 1

        # Create an n x n array initialized to 1
        array = np.zeros((n, n))

        # Default the center to the middle of the array if not provided
        if center is None:
            center = (n // 2, n // 2)

        # Generate a grid of indices
        y, x = np.ogrid[:n, :n]

        #   Calculate the distance from each grid point to the center
        distance_from_center = (x - center[1])**2 + (y - center[0])**2


        # Set values to 2 within the circle of the given radius
        array[distance_from_center <= (radius_A/grid_size)**2] = 1

        return array

    @staticmethod
    def add_small_matrix(big_matrix, small_matrix, center_i, center_j):
        """Add smmall matrix to a big matrix

        Parameters
        ----------
        big_matrix : ndarray(n,n)
            big matrix where a small  matrix will be added
        small_matrix : ndarray(m,m)
            small matrix to be added
        center_i : int
             i coordinate
        center_j : int
             j coordinate

        Returns
        -------
        ndarray(n,n)
            big matrix modified
        """
    # Calculate the top-left corner of the submatrix in big_matrix
        start_i = center_i - small_matrix.shape[0] // 2
        start_j = center_j - small_matrix.shape[1] // 2
        end_i = start_i + small_matrix.shape[0]
        end_j = start_j + small_matrix.shape[1]

    # Handle boundaries to ensure indices stay within big_matrix
        big_start_i = max(0, start_i)
        big_start_j = max(0, start_j)
        big_end_i = min(big_matrix.shape[0], end_i)
        big_end_j = min(big_matrix.shape[1], end_j)

    # Calculate the overlapping region for small_matrix
        small_start_i = big_start_i - start_i
        small_start_j = big_start_j - start_j
        small_end_i = small_start_i + (big_end_i - big_start_i)
        small_end_j = small_start_j + (big_end_j - big_start_j)

    # Add the overlapping region of small_matrix to the big_matrix
        big_matrix[big_start_i:big_end_i, big_start_j:big_end_j] += small_matrix[small_start_i:small_end_i, small_start_j:small_end_j]

        return big_matrix

    def add_defects(self,
                    matrix,
                    indexes,
                    elements,
                    names,
                    lipid,
                    mat_radii_dict):
        """Code to easily add defects in the 2d matrix

        Parameters
        ----------
        matrix : ndarray(n,n)
            Matrix where the defects are going to be added
        indexes : ndarray(i,j)
            List of indexes i,j in the matrix where the defects should be added
        elements : list
            type of element (needed to put the right radious)
        names : list
            names of the atoms (needed to map hydrophobic and not hydrophobic atoms)
        lipid : str
            lipid name
        mat_radii_dict : dict
            dictionary with the radii

        Returns
        -------
        ndarray
            matrix matrix filled with the defects
        """

        matrix = matrix

        for i in range(len(indexes[0])):

            small_matrix = mat_radii_dict[elements[i]]

            if names[i] in self.non_polar_dict[lipid]:
                small_matrix = small_matrix * 0.0001
            self.add_small_matrix(matrix, small_matrix, indexes[0][i], indexes[1][i])
        return matrix

    @staticmethod
    def extend_data(data, dimensions, percentage, others = None):
        """ Function to extent data for periodicity purposes

        Parameters
        ----------
        data : ndarray(:,2)
            data in 2D fashion
        dimensions : ndarray(2)
            periodicity dimension (usually gottem from u.ts.dimensions[:2])
        percentage : float
            Percentage of data to replicate with periodicity
        others : list(ndarray(:), ndarray(:),...,ndarray(:)), optional
            List of other data to replicate with periodicity (ndarrays initially match with data axis 0 dimension), by default None

        Returns
        -------
        ndarray
            Data extended
        ndarray, list
            Data extended, others extended
        """
        xmin = np.min(data[:,0])
        xmax = np.max(data[:,0])
        ymin = np.min(data[:,1])
        ymax = np.max(data[:,1])



        dist_x = xmax - xmin
        dist_y = ymax - ymin
        if len(data[0]) == 2:
            left_add = data[data[:,0] <= xmin + percentage*dist_x] + [dimensions[0],0]
            right_add = data[data[:,0] >= xmax - percentage*dist_x] - [dimensions[0],0]
        else:
            left_add = data[data[:,0] <= xmin + percentage*dist_x] + [dimensions[0],0,0]
            right_add = data[data[:,0] >= xmax - percentage*dist_x] - [dimensions[0],0,0]


        if others is not None:
            temp_right = []
            temp_left = []
            for i in range(len(others)):
                temp_left.append(others[i][data[:,0] <= xmin + percentage*dist_x])
                temp_right.append(others[i][data[:,0] >= xmax - percentage*dist_x])




        data = np.concatenate([data, left_add, right_add], axis = 0)

        if others is not None:
            for i in range(len(others)):
                others[i] = np.concatenate([others[i], temp_left[i], temp_right[i]], axis = 0)

        # Extent in y
        if len(data[0]) == 2:
            up_add = data[data[:,1] <= ymin + percentage*dist_y] + [0,dimensions[1]]
            low_add = data[data[:,1] >= ymax - percentage*dist_y] - [0,dimensions[1]]
        else:
            up_add = data[data[:,1] <= ymin + percentage*dist_y] + [0,dimensions[1],0]
            low_add = data[data[:,1] >= ymax - percentage*dist_y] - [0,dimensions[1],0]


        if others is not None:
            temp_up = []
            temp_low = []
            for i in range(len(others)):
                temp_up.append(others[i][data[:,1] <= ymin + percentage*dist_y])
                temp_low.append(others[i][data[:,1] >= ymax - percentage*dist_y])



        data = np.concatenate([data, up_add, low_add], axis = 0)
        if others is not None:
            for i in range(len(others)):
                others[i] = np.concatenate([others[i], temp_up[i], temp_low[i]], axis = 0)

        if others is not None:
            return data, others
        return data






    @staticmethod
    def get_highest(data, min_lenght):
        """ Code to get the highest value given two columns that are to be ordered in a 2D grid
        Parameters
        ----------
        data : (ndarray(:,2))
            Array with two columns (column1: map to a 2D grid, column2: values)
        min_lenght : (int)
            Size of squares in the 2D grid

        Returns
        -------
        ndarray(:,2)
            With the maximun of each grid square
        """


        columns = ["index", "weight"] # Data expected is an np array with columns ["index", "x", "y", "z"]
        df = pd.DataFrame(data, columns = columns)
        result = df.groupby('index', as_index=False)['weight'].max()
        result_dict = dict(zip(result['index'], result['weight']))
        hist = []
        for i in range(min_lenght):
            try:
                hist.append(result_dict[i])
            except:
                hist.append(np.nan)
        return np.array(hist)


    # Computes  and return the indexes of the data if where arranegd in a 2d histogram
    def get_indexes(self,
                    data,
                    bins = 10,
                    edges = None,
                    matrix_height = False):
        """given data in 2D, the code returns the indexes to locate the data in the 2D grid defined by edges and bins

        Parameters
        ----------
        data : ndarray(n,2 o 3)
            Data in 2D fashion, (3 columns if height = True)
        bins : int, optional
            number of bins of the grid, by default 10
        edges : list, optional
            Edges in the following way [xmin,xmax,ymin,ymax], by default None
        matrix_height : bool, optional
            returns the height matrix (matrix with only the lipids closer to water), by default False

        Returns
        -------
        tuple
            tuple with the indexes

        tuple, ndarray(bins,bins)
            tuple with indexes and 2D data of the highest point
        """



        if edges is not None:
            xmin = edges[0]
            xmax = edges[1]
            ymin = edges[2]
            ymax = edges[3]
        else:
            xmin = self.v_min
            ymin = self.v_min
            xmax = self.v_max
            ymax = self.v_max


        ran = [[xmin,xmax],[ymin,ymax]]


        nbin = np.empty(2,np.intp)
        edges = 2*[None]

        for i in range(2):
            edges[i] = np.linspace(ran[i][0], ran[i][1], bins +1)
            nbin[i] = len(edges[i]) + 1


        if not matrix_height:

            indexes = (tuple(np.searchsorted(edges[i], data[:,i], side = "right") for i in range(2)))
        else:
            indexes = (tuple(np.searchsorted(edges[i], data[:,i], side = "right") for i in range(2)))

            xy = np.ravel_multi_index(indexes, nbin)
            xy_test = xy.reshape(-1,1)

            #print("last shape", data[:,2].reshape(-1,1).shape, xy_test.shape)
            xy_test = np.concatenate((xy_test, data[:,2].reshape(-1,1)), axis = 1)
            hist = self.get_highest(xy_test, nbin.prod())


            hist = hist.reshape(nbin)
            hist = hist.astype(float, casting = "safe")
            hist[np.isnan(hist)] = 0
            core = 2*(slice(1,-1),)
            hist = hist[core]
            #print("here", hist[20,:])
            return indexes, hist

        #for i in range(2):
        #    on_edge = (data[:,i] == edges[i][-1])
        #    Ncount[i][on_edge] -= 1
        #print(np.min(Ncount[0]), np.max(Ncount[0]), np.min(Ncount[1]), np.max(Ncount[1]))
        #print(len(edges[0]), "edges len")




        return indexes

    def packing_defects_stats(self,
                                nbins = 180,
                                layer = "top",
                                edges = None,
                                periodic = False,
                                start = 0,
                                final = -1,
                                step = 1,
                                area_size = True,
                                verbose = True,
                                ):
        """ Run packing defects from `start` to `final` and stores data
          about area of defects, total area, number of defects,
          and size of defects

        Parameters
        ----------
        nbins : int, optional
            number of bins to consider, by default 180
        layer : str, optional
            layer to consider, can be top, bot, by default "top"
        edges : list(float), optional
            Edges in the shape [xmin,xmax,ymin,ymax], by default None
        periodic : bool, optional
            Consider or not periodicity, by default False
        start : int, optional
            stat frame, by default 0
        final : int, optional
            final frame, by default -1
        step : int, optional
            frames to skip, by default 1
        area_size : bool, optional
            If true return the areas of the different defects, by default True

        Returns
        -------
        pd.DataFrame, np.array
            pandas dataframe with the area information of packing defects over time and informationabout the size of the packing defects
        """

        results = list()
        sizes = list()
        for ts in self.u.trajectory[start:final:step]:
            _, packing_dict = self.packing_defects(layer = layer,
                             nbins = nbins,
                             height = False,
                             periodic = periodic,
                             edges = edges,
                             area =True,
                             count = True,
                             verbose = False)
            if verbose:
                print(f"Runing packing defects on frame {ts.frame}", end="\r")
            results.append([packing_dict["area"]["defects"],
                             packing_dict["area"]["total"],
                             packing_dict["count"]["number"],
                             ])
            if area_size:
                sizes.append(packing_dict["count"]["sizes"]*(packing_dict["grid_size"]*packing_dict["grid_size"]))
            else:
                sizes.append(packing_dict["count"]["sizes"])

        results = pd.DataFrame(results, columns = ["defects_area", "total_area", "n_defects"])

        sizes = np.concatenate(sizes, axis = 0)

        return results, sizes










    ############# End order parameters related code ###############



    ####### Code related to APL using voronoi tesselation and delunay triangulations

    # This part needs scypy to process data

    def voronoi_apl(self,
                        layer = 'top',
                        working_lip = None,
                        lipid_list = None,
                        ):
        """Computes the APL for membranes with different lipids

        Parameters
        ----------
        layer : str, optional
            layer to compute the apl, by default 'top'
        working_lip : dict, optional
            dictionary mapping lipid and atoms to work for APL, by default None
        lipid_list : list, optional
            list of lipids to be considered for APL, by default None

        Returns
        -------
        dict
            dictionary with vertices, areas, and APL
        """

        if layer == "top":
            sign = " > "
        elif layer == "bot":
            sign = " < "

        if lipid_list is None:
            lipid_list = list(self.lipid_list)

        if working_lip is None:
            working_lip = self.working_lip
        print(self.lipid_list)


        all_p = self.all_head
        positions = all_p.positions

        xmin = np.min(positions[:,0])
        xmax = np.max(positions[:,0])
        ymin = np.min(positions[:,1])
        ymax = np.max(positions[:,1])

        dist_x = xmax - xmin
        dist_y = ymax - ymin


        mean_z = positions[:,2].mean()

        selection_string = f"(((resname {lipid_list[0]} and name {working_lip[lipid_list[0]]['head']}) and prop z {sign} {mean_z}))"
        for lipid in lipid_list[1:]:
            selection_string += f" or (((resname {lipid} and name {working_lip[lipid]['head']}) and prop z {sign} {mean_z}))"

        print(selection_string)


        heads = self.memb.select_atoms(selection_string)
        heads_pos = heads.positions[:,:2]
        resnames_pos = heads.resnames
        orig_len = len(heads_pos)
        print("Here first shapes", heads_pos.shape, resnames_pos.shape)

        ## Extent data
        dimensions = self.u.trajectory.ts.dimensions[:2]
        cons = 0.1

        heads_pos, resnames_pos = self.extend_data(heads_pos, dimensions, cons, others = [resnames_pos])
        resnames_pos = resnames_pos[0]
        from scipy.spatial import Voronoi, voronoi_plot_2d
        from scipy.spatial import ConvexHull

        voronoi_dict = {"vertices":list(),
                        "points":heads_pos,
                        "areas":list(),
                        "resnames":resnames_pos,
                        "orig_len":orig_len
                         }


        voronoi = Voronoi(heads_pos)
        vertices = voronoi.vertices

        resnames = list(set(voronoi_dict["resnames"]))
        result_dict = {}
        for lipid in resnames:
            result_dict[lipid] = list()

        #for i, region in enumerate(voronoi.point_region[:orig_len]):
        update_points = list()
        for i, region in enumerate(voronoi.point_region):
            #print(i, region, len(voronoi.point_region), len(voronoi_dict["resnames"]))
            if -1 in voronoi.regions[region]:
                #print("here")
                continue

            vertex = vertices[voronoi.regions[region]]
            hull = ConvexHull(vertex)
            area = hull.volume
            voronoi_dict["areas"].append(area)
            update_points.append(voronoi_dict["points"][i])
            if i < orig_len:
                result_dict[voronoi_dict["resnames"][i]].append(area)

            voronoi_dict["vertices"].append(vertex)
        voronoi_dict["points"] = np.array(update_points)
        print(len(voronoi_dict["areas"]))

        for lipid in resnames:
            result_dict[lipid] = np.mean(np.array(result_dict[lipid]))

        voronoi_dict["apl"] = result_dict

        #voronoi_plot_2d(voronoi)
        #for region in voronoi.regions:
        #    if not -1 in region:
        #        polygon = [voronoi.vertices[i] for i in region]
        #        plt.fill(*zip(*polygon))

        #plt.show()



        return voronoi_dict


    def map_voronoitest(self,voronoi_vertices, voronoi_areas, nbins, dimensions):

        from shapely.geometry import Polygon, Point
        from shapely.strtree import STRtree

        voronois = [Polygon(vertices) for vertices in voronoi_vertices]
        tree = STRtree(voronois)

        xmin =dimensions[0]
        xmax =dimensions[1]
        ymin =dimensions[2]
        ymax =dimensions[3]



        xcoords = np.linspace(xmin, xmax, nbins)
        ycoords = np.linspace(ymin, ymax, nbins)

        xx, yy = np.meshgrid(xcoords, ycoords)
        points = np.vstack([xx.ravel(), yy.ravel()])
        grid_values = np.array([self.process_point(point, tree, voronois, voronoi_areas) for point in points.T])
        grid = grid_values.reshape(nbins, nbins)
        plt.imshow(grid, cmap = "Spectral")
        plt.show()

    def map_voronoi(self, voronoi_points, voronoi_areas, nbins, edges):
        """ Function to map voronoi APL to a 2D plane

        Parameters
        ----------
        voronoi_points : ndarray(:,2)
            [x,y] positions of the points to be considered in the voronoi plot
        voronoi_areas : ndarray(:)
            Areas correspondng to the points
        nbins : int
            number of bins for the grid
        edges : list
            A list with the lipids of the grid [xmin,xmax,ymin,ymax]

        Returns
        -------
        ndarray, edges
            numpy array (nbins,nbins), adn edges of this array
        """
        xmin =edges[0]
        xmax =edges[1]
        ymin =edges[2]
        ymax =edges[3]

        xcoords = np.linspace(xmin, xmax, nbins)
        ycoords = np.linspace(ymin, ymax, nbins)

        xx, yy = np.meshgrid(xcoords, ycoords)
        grid_points = np.vstack([xx.ravel(), yy.ravel()]).T
        points = voronoi_points


        distances = np.linalg.norm(grid_points[:,None, :]- points[None,:,:], axis = 2)


        closest_seed_indices = np.argmin(distances, axis=1).astype(int)

        voronoi_areas = np.array(voronoi_areas)
        grid = voronoi_areas[closest_seed_indices].reshape(nbins, nbins)

        return grid, edges



    def grid_apl(self,layer = "top", start = 0, final = -1, step = 1, lipid_list = None, nbins = 180, edges = None):
        """Function to compute and map the grid APL for several frames, map them to a 2D grid and average them

        Parameters
        ----------
        layer : str, optional
            working lipid layer, by default "top"
        start : int, optional
            Frame to start, by default 0
        final : int, optional
            final frame, by default -1
        step : int, optional
            Frames to skip, by default 1
        lipid_list : list, optional
            lipids involved in the computation, by default None
        nbins : int, optional
            number of bins for the grid, by default 180
        edges : list, optional
            A list with the limits of the grid [xmin,xmax,ymin,ymax]

        Returns
        -------
        ndarray
            Array with the averaged 2D APL, edges
        """
        if lipid_list == None:
            lipid_list = list(self.lipid_list)
        if layer == "top":
            sign = " > "
        elif layer == "bot":
            sign = " < "


        all_p = self.all_head
        positions = all_p.positions
        mean_z = positions[:,2].mean()

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



        #xmin1 = np.min(positions[:,0])
        #xmax1 = np.max(positions[:,0])
        #ymin1 = np.min(positions[:,1])
        #ymax1 = np.max(positions[:,1])

        #dist_x = xmax1 - xmin1
        #dist_y = ymax1 - ymin1

        matrices = []
        for ts in self.u.trajectory[start:final:step]:

            selection_string = f"(((resname {lipid_list[0]} and name {self.working_lip[lipid_list[0]]['head']}) and prop z {sign} {mean_z}))"
            for lipid in lipid_list[1:]:
                selection_string += f" or (((resname {lipid} and name {self.working_lip[lipid]['head']}) and prop z {sign} {mean_z}))"

            heads = self.memb.select_atoms(selection_string)

            head_pos = heads.positions[:2]
            voronoi_dict = self.voronoi_apl(layer = layer)
            matrix,_ = self.map_voronoi(voronoi_dict["points"], voronoi_dict["areas"], nbins, [xmin, xmax, ymin, ymax])
            matrices.append(matrix)
            print(matrix.shape)




        final_mat = np.mean(np.array(matrices), axis = 0)

        #print(f"final_mat {np.array(matrices).shape}")
        #fig , ax = plt.subplots(1,len(matrices)+1)
        #count = 0
        #indices = []
        #for mat in matrices:
            #ax[count].imshow(mat, cmap = "Spectral")
        #    x_ref = matrices[0].flatten()
        #    y = mat.flatten()
            #ax[count].scatter(x_ref,y)

        #    indices.append(np.corrcoef(x_ref, y)[0,1])
        #    count +=1
        #ax[count].imshow(final_mat, cmap = "Spectral")
        #plt.show()


        #plt.plot(indices)
        #plt.show()
        return final_mat, edges

    def windows_apl(self, layer = "top", start = 0, final = -1, step = 1, w_size = 10, lipid_list = None, nbins = 180, edges = None):
        """Function to compute APL and map it to a 2D grid in windows of time defined by the user.

        Parameters
        ----------
        layer : str, optional
            Wroking layer, by default "top"
        start : int, optional
            Start Frame, by default 0
        final : int, optional
            Final frame, by default -1
        step : int, optional
            Frames to skip, by default 1
        w_size : int, optional
            windows size (number of frames), by default 10
        lipid_list : list, optional
            lipids to be included in the analysis, by default None
        nbins : int, optional
            nimber of bins in the grid, by default 180
        edges : list, optional
            list with the edges of the grid [xmin,xmax,ymin,ymax], by default None

        Returns
        -------
        list(ndarrays(nbins,bins)), list
            list with the windows averages between the start time and the final time, and edges of these matrices
        """


        if final == -1:
            final = len(self.u.trajectory)
        n_windows = (final-start)//w_size
        matrices = []
        for i in range(n_windows):
            matrix = self.grid_apl(layer = layer, start = i*w_size + start, final =(i+1)*w_size + start , step = step, lipid_list = None, nbins = nbins)
            matrices.append(matrix)
        return matrices






















    def process_point(self, point, tree, polygons, values):
        from shapely.geometry import Polygon, Point
        point_geom = Point(point[0], point[1])
        polygon_value_map = dict(zip(polygons, values))
        candidates = tree.query(point_geom)
        for candidate in candidates:
            if polygons[candidate].contains(point_geom):
                return values[candidate]
        return 0














        #fig = voronoi_plot_2d(voronoi)
        #plt.xlim(xmin,xmax)
        #plt.ylim(ymin,ymax)
        #plt.show()

        #return heads_pos














