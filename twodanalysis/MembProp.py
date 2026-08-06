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
import re

def natural_key(node_name):
    """
    Key function to sort strings like C1, C2, C10 numerically.
    Extracts the prefix and numeric part.
    """
    match = re.match(r"([A-Za-z]+)(\d+)", str(node_name))
    if match:
        prefix, number = match.groups()
        return (prefix, int(number))
    else:
        return (node_name, 0)  # fallback


class MembProp:

    def __init__(
                self,
                obj,
                lipid_list = None,
                connection_chains = None,
                working_lip = None,
                verbose = False,
                forcefield = "charmm",
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
        connection_chains : dict, optional
            Dictionary with the information of the atoms that connect the lipid tails and the hidrophilic head.
            This dictionary should be as follows {"CHL1" : [("O3", "C3")], "DODMA" : [("C21", "C22"), ("C31", "C32")],}
            and should be provided for each lipid. Defaults to None (If None we use ([("C21", "C22"), ("C31", "C32")] for most lipids))
        working_lip : dict, optional
            Dictionary with the information of the lipid heads to use for the multiple analysis. The dictionary should be as follows
            {"CHL1" : {"head" : "O3"}, "DODMA" : {"head" : "N1"}, "DSPC" : {"head":"P"}} and should have entries for all the lipids
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
            self.lipid_list = list(set(self.memb.residues.resnames))


        # Set percentage for periodicity
        self.periodicity = 0.1


        # Set radius sizes of different elements
        self.radii_dict = None
        self.forcefield = forcefield
        self.chain_info = {}
        self.chain_conections = {}
        self.non_polar_dict = {}
        self.polar_dict = {}
        self.first_lipids = {}
        self.non_polar_visualize = {}

        self.amber_heads = ["PC", "PE", "PS", "PGR", "PGS", "PH", "SPM"]


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

        #self.connection_chains = self.connection_chains if connection_chains is None else connection_chains
        if connection_chains is not None:
            for lip in connection_chains:
                self.connection_chains[lip] = connection_chains[lip]

        for lipid in self.lipid_list:
            test_lips = self.working_lip.keys()
            if lipid not in test_lips:
                self.working_lip[lipid] = {"head" : "P", "charge" : 0}
            if lipid not in list(self.connection_chains.keys()):
                self.connection_chains[lipid] = [("C21", "C22"), ("C31", "C32")]


        self.working_lip = self.working_lip if working_lip is None else working_lip


        head_atoms = set([self.working_lip[lipid]["head"] for lipid in self.working_lip])
        if "P" in head_atoms:
            self.all_head = self.memb.select_atoms(self.build_resname(self.lipid_list) + " and name P")
        else:
            self.all_head = self.memb.select_atoms(self.build_resname_head(self.lipid_list))


        self.map_layers = {"top" : " > ",
                           "bot": " < ",
                           "up" : " > ",
                           "down" : " < "}



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
            self.radii_dict = {"H": 1.1,
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


    #def infer_chain_names(self, lipid, chain):
    #    lipid_w = self.memb.select_atoms(f"resname {lipid}")



    def guess_chain_lenght(self):
        # Guess the chain length of lipids. Chain sn2 start with C2 and chain sn1 start with C3
        for lipid in self.lipid_list:
            first_lipid = self.memb.select_atoms(f"resname {lipid}").resids[0]
            #print("###########",self.memb.select_atoms(f"resname {lipid}").resids)
            #print("first lipid", first_lipid)
            chain_l = []
            tail = self.memb.select_atoms(f"resid {first_lipid}")


            if self.forcefield == "charmm":
                for bond in self.connection_chains[lipid]:

                    _ , chaintemp = self.cut_structure(
                        tail,
                        bond,
                        )
                    c_in_node = [n for n in chaintemp.nodes if "C" in str(n)]
                    chain_l.append(len(c_in_node))
                #print(chain_l, lipid)
                try:
                    if len(chain_l) == 2:
                        chain_l = [chain_l[1]] + [chain_l[0]]
                    else:
                        chain_l = [chain_l[1]] +[chain_l[0]]  + [chain_l[2:]]
                except:
                    pass
                self.chain_info[lipid] = chain_l

            elif self.forcefield == "amber":
                if lipid != "CHL":

                    chain1 = tail.resids[0] - 1
                    chain2 = tail.resids[0] + 1

                    atoms_c1 = self.u.select_atoms(f"resid {chain1}")
                    atoms_c2 = self.u.select_atoms(f"resid {chain2}")
                    c1_in_node = [n for n in atoms_c1.names if "C" in str(n)]
                    c2_in_node = [n for n in atoms_c2.names if "C" in str(n)]
                    chain_l = [len(c1_in_node), len(c2_in_node)]
                else:
                    chain_l = [7]
                self.chain_info[lipid] = chain_l

            #print("################", self.chain_info)



    def guess_last_cs(self):

        if self.forcefield == "charmm":
            for lipid in self.lipid_list:
                first_lipid = self.memb.select_atoms(f"resname {lipid}").resids[0]
                actual_sn1 = self.memb.select_atoms(f"resid {first_lipid} and name C3*")
                actual_sn2 = self.memb.select_atoms(f"resid {first_lipid} and name C2*")
                actual_sn1 = actual_sn1.names
                actual_sn2 = actual_sn2.names
                self.working_lip[lipid]["last_c"] = [actual_sn1[-1], actual_sn2[-1]]
        elif self.forcefield == "amber":
            for lipid in self.lipid_list:
                structure = self.extract_chain_info(lipid)
                carbons_1 = [atom for conn in structure[0] for atom in conn if "C" in atom]
                carbons_2 = [atom for conn in structure[0] for atom in conn if "C" in atom]
                self.working_lip[lipid]["last_c"] = [carbons_1[-1], carbons_2[-1]]




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



    @staticmethod
    def build_resname(resnames_list):
        resnames_list = list(resnames_list)
        string = " (resname " + resnames_list[0]
        for resname in resnames_list[1:]:
            string = string + " or resname " + resname
        string = string + ") "
        return string


    def build_resname_head(self, resnames_list):
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
        if len(resnames_list) > 1:
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



        self.non_polar_dict = {}
        self.polar_dict = {}
        self.first_lipids = {}
        self.non_polar_visualize = {}

        for lipid in self.lipid_list:
            first_lipid = self.memb.select_atoms(f"resname {lipid}")
            first_id = first_lipid.resids[0]
            first_lipid = first_lipid.select_atoms(f"resid {first_id}")
            self.first_lipids[lipid] = first_id



            if self.forcefield == "charmm":
                G = self.atgroup2nx(first_lipid)
                non_polar = []
                for edge in self.connection_chains[lipid]:
                    G.remove_edge(edge[0], edge[1])
                    components = list(nx.connected_components(G))
                    non_polar.append(next(comp for comp in components if edge[1] in comp)) # Stores the tails

                self.non_polar_dict[lipid] = [elem for s in non_polar for elem in s]
                polar = next(comp for comp in components if self.connection_chains[lipid][0][0] in comp)
                self.polar_dict[lipid] = polar

            elif self.forcefield == "amber":


                if lipid != "CHL":
                    heads_ats = first_lipid
                    self.polar_dict[lipid] = heads_ats.names

                    tail1 = self.u.select_atoms(f"resid {first_id -1}")
                    tail2 = self.u.select_atoms(f"resid {first_id +1}")
                    #print("Here: ########", tail1.names, tail2.names, list(tail1.names) + list(tail2.names))
                    non_polar = list(set(list(tail1.names) + list(tail2.names)))
                    self.non_polar_dict[lipid] = non_polar
                else:
                    self.polar_dict[lipid] = ["O1", "HO1"]
                    non_polar = first_lipid.names
                    non_polar_final = [name for name in non_polar if name not in self.polar_dict[lipid]]
                    self.non_polar_dict[lipid] = non_polar_final




            #print(lipid , ":", self.polar_dict[lipid])

    def _store_fist_lipids(self):
        self.first_lipids = {}
        for lipid in self.lipid_list:
            first_lipid = self.memb.select_atoms(f"resname {lipid}")
            first_id = first_lipid.resids[0]
            first_lipid = first_lipid.select_atoms(f"resid {first_id}")
            self.first_lipids[lipid] = first_id

    def nx_chain_names(self):

        self.first_lipids = {}
        carbon_conections = {}

        for lipid in self.lipid_list:
            first_lipid = self.memb.select_atoms(f"resname {lipid}")
            first_id = first_lipid.resids[0]
            first_lipid = first_lipid.select_atoms(f"resid {first_id}")
            self.first_lipids[lipid] = first_id


            if lipid != "CHL1":


                G = self.atgroup2nx(first_lipid)


                chains = []
                for edge in self.connection_chains[lipid]:
                    G.remove_edge(edge[0], edge[1])
                    components = list(nx.connected_components(G))
                    chains.append(next(comp for comp in components if edge[1] in comp)) # Store the tail of the lipids

                chain_graphs = []
                carbons = []
                for chain in chains:
                    chain_graphs.append(G.subgraph(chain))
                    filtered_sorted_carbs = [item for item in chain if "C" in item]
                    carbons.append(filtered_sorted_carbs)


                carbon_conections = {}

                for chain, carbs in zip(chain_graphs,carbons):
                    for carbon in carbs:
                        connected_atoms = list(chain.neighbors(carbon))
                        connected_atoms = [atom for atom in connected_atoms if "C" not in atom]
                        carbon_conections[carbon] = connected_atoms


                #print("here", carbon_conections)

            self.chain_conections[lipid] = carbon_conections

        #print(self.chain_conections)

    def cut_structure(self, structure, bond):
        """Tool needed to cut structures in two parts based on a bond.
        This function is used to cut the structure in two parts based on a bond.
        It uses NetworkX to create a graph of the structure and then removes the bond.
        It then finds the connected components of the graph and returns the two parts.
        If the bond does not cut the structure into two parts, it raises a ValueError.


        Parameters
        ----------
        structure : MDAnalysis Universe or AtomGroup
            Any structure that we aim to cut
        bond : tuple
            tuple with the names of the atoms of the bond we want to cut
            e.g. ("C21", "C22") or ("C31", "C32")

        Returns
        -------
        part1, part2 : NetworkX Graph
            Two parts of the structure after cutting it based on the bond.
        Raises
        ------
        ValueError
            If the bond does not cut the structure into two parts, it raises a ValueError.

        """

        G = self.atgroup2nx(structure)
        G.remove_edge(bond[0], bond[1])
        components = list(nx.connected_components(G))
        subgraphs = [G.subgraph(c).copy() for c in components]
        if len(subgraphs) != 2:
            raise ValueError("The bond does not cut the structure into two parts")

        part1 = subgraphs[0] if bond[0] in subgraphs[0].nodes else subgraphs[1]
        part2 = subgraphs[1] if bond[1] in subgraphs[1].nodes else subgraphs[0]

        return part1, part2

    def atgroup2nx(self, atgroup):

        """Function to convert an MDAnalysis AtomGroup to a NetworkX Graph.

        Parameters
        ----------
        atgroup : MDAnalysis AtomGroup
            AtomGroup to convert

        Returns
        -------
        NetworkX Graph
            Graph representing the structure of the AtomGroup
        """
        try:
            bonds = atgroup.bonds.indices
            #print(bonds)
        except:
            bonds = mda.topology.guessers.guess_bonds(atgroup, atgroup.positions)

        # Map each atom's global (universe) index directly to its name. Both
        # ``atgroup.indices`` and ``atgroup.names`` are ordered the same way, so
        # they refer to the same atom position by position. Building the mapping
        # this way (instead of positionally zipping the sorted set of bonded-atom
        # indices with the names) keeps the labels correct even when some atoms
        # have no bonds, when the atomgroup is not sorted by index, or when a bond
        # reaches an atom that is not part of the atomgroup.
        map_node_name = dict(zip(atgroup.indices, atgroup.names))

        edges_labeled = [
            (map_node_name[a], map_node_name[b])
            for a, b in bonds
            if a in map_node_name and b in map_node_name
        ]

        G = nx.Graph()
        G.add_edges_from(edges_labeled)

        return G


    def chain_structure(self, chain):
        """Function to get the structure of the chain in the shape [C,Hi,Hj] for any nomemclature.

        Parameters
        ----------
        chain : NetworkX Graph
            Graph representing a lipid tail

        Returns
        -------
        list
            List containing various list of the form [C, Hi, Hj] where C is the carbon and Hi and Hj are the hydrogens connected to it.
        """
        c_in_node = [n for n in chain.nodes if "C" in str(n)]
        c_in_node = sorted(c_in_node, key=natural_key)
        neighbors = [list(chain.neighbors(n)) for n in c_in_node]
        neighbors = [[c_in_node[i]] + [item for item in neighbors[i] if "C" not in item] for i in range(len(c_in_node))]
        return neighbors

    def extract_chain_info(self, lipid):
        self._store_fist_lipids()
        lipids_ids = self.first_lipids


        chains = []

        if self.forcefield == "charmm":
            for chain in self.connection_chains[lipid]:

                _, chaintemp = self.cut_structure(
                        self.memb.select_atoms(f"resid {lipids_ids[lipid]}"),
                        chain,
                        )
                chains.append(chaintemp)

        elif self.forcefield == "amber":
            if lipid != "CHL":
                chain1 = lipids_ids[lipid] - 1
                chain2 = lipids_ids[lipid] + 1

                atoms_c1 = self.u.select_atoms(f"resid {chain1}")
                atoms_c2 = self.u.select_atoms(f"resid {chain2}")

                Ch1 = self.atgroup2nx(atoms_c1)
                Ch2 = self.atgroup2nx(atoms_c2)

                chains = [Ch1, Ch2]
            if lipid == "CHL" and "CHL" in self.connection_chains:
                for chain in self.connection_chains[lipid]:

                    _, chaintemp = self.cut_structure(
                            self.memb.select_atoms(f"resid {lipids_ids[lipid]}"),
                            chain,
                            )

                    chains.append(chaintemp)


        return [self.chain_structure(tail) for tail in chains]



    @staticmethod
    def print_dict(dictio):
        string = ""
        for key in dictio.keys():
            string += f"{key} : {dictio[key]}\n"
        return string

