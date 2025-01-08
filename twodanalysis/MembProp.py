"""
MembProp
=============
Class created mainly to analyze lipid membranes in different ways

Classes
-------

.. autoclass:: MembProp
    :members:
    :undoc-members:
    :show-inheritance:


"""

import MDAnalysis as mda
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


class MembProp:

    def __init__(
                self,
                obj,
                lipid_list = None,
                verbose = False,
            ):
        """Class set up common properties and to automatically guess useful membrane properties
         useful for the analyses that can ber performed. It sets up lipid heads (used to identify membrane middle),
         periodicity (used to introduce periodicity), default edges for the 2D maps. On the other hand it helps with the guess
        of properties such as chain length for order parameter, last carbon for splay angle, polarity information and radii
        for packing defects.


        Parameters
        ----------
        obj : Universe or AtomGroup
            MDanalysis universe or AtomGroup with or without trajectory
        lipid_list : List of lipids in the membrane, optional
            If not provided, the code will guess them based on eliminating protein and RNA
            from the atoms. If you upload a trajectory with water and ions and does not specify this
            information it can lead to errors, by default None
        edges : list(str), optional
            edges in the for [xmin,xmax,ymin,ymax]
        verbose : bool, optional
            Print information of the membrane, by default False
        """



        # Read trajectory depending if tpr is provided or not
        if isinstance(obj, mda.Universe):
            self.u = obj
        elif isinstance(obj,mda.core.groups.AtomGroup):

            self.u = mda.core.universe.Merge(obj)
        else:
            raise TypeError("Input must be an MDAnalysis Universe or AtomGroup")




        # Select elements in the membrane (in principle only lipids)
        if lipid_list is None: # Select only elements of the membrane
            self.memb = self.u.select_atoms("all and not protein and not\
                                             (resname URA or resname GUA\
                                             or resname ADE or resname CYT\
                                             or resname THY)")
            self.lipid_list = list(set(self.memb.residues.resnames))
        else:
            self.memb = self.u.select_atoms(f"{self.build_resname(list(lipid_list))}")


        # Set percentage for periodicity
        self.periodicity = 0.1


        # Set radius sizes of different elements
        self.radii_dict = None
        self.chain_info = {}

        self.working_lip = {
                                "CHL1" : {"head" :"O3", "charge" : 0},
                                "DODMA" : {"head" :"N1", "charge" : -0.21},
                                "DSPC" : {"head" :"P", "charge" : 1.1},
                            }

        self.connection_chains = {

            "CHL1" : [("O3", "C3")],
            "DODMA" : [("C21", "C22"), ("C31", "C32")],
            "DSPC" : [("C21", "C22"), ("C31", "C32")],

        }


        for lipid in self.lipid_list:
            test_lips = self.working_lip.keys()
            if lipid not in test_lips:
                self.working_lip[lipid] = {"head" : "P", "charge" : 0}
                self.connection_chains[lipid] = [("C21", "C22"), ("C31", "C32")]


        self.all_head = self.memb.select_atoms(self.build_resname(self.lipid_list) + " and name P")





        self.verbose = verbose
        if verbose:
            print(f"This system contains the following lipids : {self.lipid_list}\n\n")
            print(f"The chain length is : \n{self.print_dict(self.chain_info)}\n")
            print(f"We will use the following heads and charges for the following lipids.\
                   If the lipid is not here we will use P as head as default \n{self.print_dict(self.working_lip)}\n")
            print("Note: To compute the middle of the membrane we use only P heads\n\n")
            #print(f"The default start frame is {self.start}, final {self.final}, step {self.step}\n\n")


    def add_radii(self, radii_dict = None):
        if radii_dict is None:
            radii_dict = {"H": 1.1,
                            "N": 1.55,
                            "C": 1.7,
                            "P": 1.8,
                            "O": 1.52,
                            "S" : 1.8,
                            }
        # Add radii as a topology attribute for Mdanalysis
        try:
            string_array = self.memb.elements
        except:
            guessed_elements = mda.topology.guessers.guess_types(self.u.atoms.names)
            self.u.add_TopologyAttr("elements", guessed_elements)
            string_array = self.memb.elements
        radii_array = np.array([self.radii_dict[element] for element in string_array])
        self.u.add_TopologyAttr("radii")
        self.memb.radii = radii_array



    def guess_chain_lenght(self):
        # Guess the chain length of lipids. Chain sn2 start with C2 and chain sn1 start with C3
        for lipid in self.lipid_list:
            first_lipid = self.memb.select_atoms(f"resname {lipid}").resids[0]
            actual_sn1 = self.memb.select_atoms(f"resid {first_lipid} and name C3*")
            actual_sn2 = self.memb.select_atoms(f"resid {first_lipid} and name C2*")
            actual_sn1 = actual_sn1.names
            actual_sn2 = actual_sn2.names
            self.chain_info[lipid] = [len(actual_sn1) - 2, len(actual_sn2) - 2]

    def guess_last_cs(self):
        for lipid in self.lipid_list:
            first_lipid = self.memb.select_atoms(f"resname {lipid}").resids[0]
            actual_sn1 = self.memb.select_atoms(f"resid {first_lipid} and name C3*")
            actual_sn2 = self.memb.select_atoms(f"resid {first_lipid} and name C2*")
            actual_sn1 = actual_sn1.names
            actual_sn2 = actual_sn2.names
            self.working_lip[lipid]["last_c"] = [actual_sn1[-1], actual_sn2[-1]]
        return self.working_lip


    def guess_polarity(self):
        self.non_polar_dict = {}
        self.first_lipids = {}
        self.non_polar_visualize = {}
        for lipid in self.lipid_list:
            first_lipid = self.memb.select_atoms(f"resname {lipid}").resids[0]
            self.first_lipids[lipid] = first_lipid
            if lipid == "CHL1":
                non_polar = self.memb.select_atoms(f"resid {first_lipid} and not (name O3 or name H3')")
            else:
                non_polar = self.memb.select_atoms(f"resid {first_lipid} and (name *C3* or name H*Y or name H*X or name H*Z  or name *C2* or name H*R or name H*S or name H*T) and not (name C3 or name C31 or name HY or name HX or name HZ  or name C2 or name C21 or name HR or name HS or name HT)")
            self.non_polar_dict[lipid] = list(non_polar.names)
            self.non_polar_visualize[lipid] = [all_lip, non_polar]


    @staticmethod
    def build_resname(resnames_list):
        resnames_list = list(resnames_list)
        string = " (resname " + resnames_list[0]
        for resname in resnames_list[1:]:
            string = string + " or resname " + resname
        string = string + ") "
        return string

    @staticmethod
    def build_resname_head(resnames_list):
        resnames_list = list(resnames_list)
        string = f"( (resname {resnames_list[0]}  and name {self.working_lip[resnames_list[0]]['head']}) "

        for resname in resnames_list[1:]:
            string += f" or (resname {resname}  and name {self.working_lip[resname]['head']}) "

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





    def nx_polarity(self):
        import matplotlib.pyplot as plt


        self.non_polar_dict = {}
        self.polar_dict = {}
        self.first_lipids = {}
        self.non_polar_visualize = {}

        for lipid in self.lipid_list:
            first_lipid = self.memb.select_atoms(f"resname {lipid}")
            first_id = first_lipid.resids[0]
            first_lipid = first_lipid.select_atoms(f"resid {first_id}")
            self.first_lipids[lipid] = first_id


            try:
                bonds = first_lipid.bonds.indices
                #print(bonds)
            except:
                bonds = mda.topology.guessers.guess_bonds(first_lipid, first_lipid.positions)
            #print(bonds.shape)
            nodes = set()
            for edge in bonds:
                nodes.update(edge)
            nodes = sorted(nodes)

            names = first_lipid.names
            map_node_name = {nodes[i]:names[i] for i in range(len(nodes))}

            edges_labeled = [(map_node_name[u], map_node_name[v]) for u,v in bonds]

            G = nx.Graph()
            G.add_edges_from(edges_labeled)

            non_polar = []

            for edge in self.connection_chains[lipid]:
                G.remove_edge(edge[0], edge[1])

                components = list(nx.connected_components(G))

                non_polar.append(next(comp for comp in components if edge[1] in comp))
                if lipid == "CHL1":
                    print(components, next(comp for comp in components if edge[1] in comp))

            self.non_polar_dict[lipid] = [elem for s in non_polar for elem in s]

            polar = next(comp for comp in components if self.connection_chains[lipid][0][0] in comp)

            self.polar_dict[lipid] = polar

            print(lipid , ":", self.polar_dict[lipid])


    def nx_chain_names(self):
        import matplotlib.pyplot as plt



        self.first_lipids = {}

        self.chain_conections = {}
        for lipid in self.lipid_list:
            first_lipid = self.memb.select_atoms(f"resname {lipid}")
            first_id = first_lipid.resids[0]
            first_lipid = first_lipid.select_atoms(f"resid {first_id}")
            self.first_lipids[lipid] = first_id


            if lipid != "CHL1":
                bonds = mda.topology.guessers.guess_bonds(first_lipid, first_lipid.positions)
                nodes = set()
                for edge in bonds:
                    nodes.update(edge)
                nodes = sorted(nodes)

                names = first_lipid.names
                map_node_name = {nodes[i]:names[i] for i in range(len(nodes))}

                edges_labeled = [(map_node_name[u], map_node_name[v]) for u,v in bonds]

                G = nx.Graph()
                G.add_edges_from(edges_labeled)

                chains = []
                for edge in self.connection_chains[lipid]:
                    G.remove_edge(edge[0], edge[1])
                    components = list(nx.connected_components(G))
                    chains.append(next(comp for comp in components if edge[1] in comp))

                chain_graphs = []
                carbons = []
                for chain in chains:
                    chain_graphs.append(G.subgraph(chain))
                    filtered_sorted_carbs = [item for item in chain if "C" in item]
                    carbons.append(filtered_sorted_carbs)
                nx.draw(G, with_labels = True)
                plt.show()

                carbon_conections = {}

                for chain, carbs in zip(chain_graphs,carbons):
                    for carbon in carbs:
                        connected_atoms = list(chain.neighbors(carbon))
                        connected_atoms = [atom for atom in connected_atoms if "C" not in atom]
                        carbon_conections[carbon] = connected_atoms




            self.chain_conections[lipid] = carbon_conections

        print(self.chain_conections)





    @staticmethod
    def print_dict(dictio):
        string = ""
        for key in dictio.keys():
            string += f"{key} : {dictio[key]}\n"
        return string

