"""
Voronoi2D
=============
Class created mainly to analyze lipid membranes in different ways

Classes
-------

.. autoclass:: Voronoi2D
    :members:
    :undoc-members:
    :show-inheritance:


"""


import MDAnalysis as mda
import numpy as np
import pandas as pd
from twodanalysis import MembProp
from twodanalysis.analysis import OrderParameters



class Voronoi2D(MembProp):

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




    ####### Code related to APL using voronoi tesselation and delunay triangulations
    # This part needs scypy to process data
    def voronoi_properties(self,
                    layer = 'top',
                        working_lip = None,
                        lipid_list = None,
                        splay = False,
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
        #print(self.lipid_list)


        all_p = self.all_head
        positions = all_p.positions






        mean_z = positions[:,2].mean()

        selection_string = f"(((resname {lipid_list[0]} and name {working_lip[lipid_list[0]]['head']}) and prop z {sign} {mean_z}))"
        for lipid in lipid_list[1:]:
            selection_string += f" or (((resname {lipid} and name {working_lip[lipid]['head']}) and prop z {sign} {mean_z}))"






        #print(selection_string)


        heads = self.memb.select_atoms(selection_string)
        heads_pos = heads.positions[:,:2]
        height_pos =np.abs(heads.positions[:,2] - mean_z)
        resnames_pos = heads.resnames
        orig_len = len(heads_pos)
        #print("Here first shapes", heads_pos.shape, resnames_pos.shape)

        ## Extent data
        dimensions = self.u.trajectory.ts.dimensions[:2]

        others = [resnames_pos, height_pos]


        if splay:
            carbons1 = [self.working_lip[lipid]["last_c"][0] for lipid in lipid_list]
            carbons2 = [self.working_lip[lipid]["last_c"][1] for lipid in lipid_list]
            heads = [self.working_lip[lipid]["head"] for lipid in lipid_list]
            heads = self.build_resname_atom(lipid_list, heads)
            carbons1 = self.build_resname_atom(lipid_list, carbons1)
            carbons2 = self.build_resname_atom(lipid_list, carbons2)
            selection_string_byres = f"byres {selection_string}"
            lipid_ats = self.memb.select_atoms(selection_string_byres)


            head_p = lipid_ats.select_atoms(heads)
            c1 = lipid_ats.select_atoms(carbons1)
            c2 = lipid_ats.select_atoms(carbons2)
            v1 = c1.positions - head_p.positions
            v2 = c2.positions - head_p.positions
            #print(heads,carbons2,carbons1, head_p.n_atoms, c1.n_atoms)
            costheta = np.sum(v1 * v2, axis=1)/(np.linalg.norm(v1, axis = 1)* np.linalg.norm(v2, axis = 1))# Compute the cos of splay angle, must bet lenmght nlipids
            costheta = np.rad2deg(costheta)
            others.append(costheta)


        heads_pos, others = self.extend_data(heads_pos, dimensions, self.periodicity, others = others)
        resnames_pos = others[0]
        height_pos = others[1]
        from scipy.spatial import Voronoi, voronoi_plot_2d
        from scipy.spatial import ConvexHull

        voronoi_dict = {"vertices":list(),
                        "points":heads_pos,
                        "heights":height_pos,
                        "areas":list(),
                        "resnames":resnames_pos,
                        "orig_len":orig_len
                         }
        if splay:
            voronoi_dict["splay"] = others[2]


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
                voronoi_dict["areas"].append(np.nan)
                continue

            vertex = vertices[voronoi.regions[region]]
            hull = ConvexHull(vertex)
            area = hull.volume
            voronoi_dict["areas"].append(area)
            #update_points.append(voronoi_dict["points"][i])
            if i < orig_len:
                result_dict[voronoi_dict["resnames"][i]].append(area)

            voronoi_dict["vertices"].append(vertex)
        #voronoi_dict["points"] = np.array(update_points)
        #print(len(voronoi_dict["areas"]))

        for lipid in resnames:
            result_dict[lipid] = np.mean(np.array(result_dict[lipid]))

        voronoi_dict["apl"] = result_dict


        return voronoi_dict


    def map_voronoi(self, voronoi_points, voronoi_property, nbins, edges):
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

        voronoi_property = np.array(voronoi_property)
        grid = voronoi_property[closest_seed_indices].reshape(nbins, nbins)

        return grid, edges



    def voronoi_apl(self,layer = "top", start = 0, final = -1, step = 1, lipid_list = None, nbins = 180, edges = None):
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

        if edges is None:
            xmin = self.v_min
            xmax = self.v_max
            ymin = self.v_min
            ymax = self.v_max
            edges = [xmin,xmax,ymin,ymax]
        else:
            xmin = edges[0]
            xmax = edges[1]
            ymin = edges[2]
            ymax = edges[3]

        grid_size = (xmax-xmin)/nbins
        ng = int(nbins*0.05)
        self.guess_last_cs()
        no_present = [lipid for lipid in list(self.lipid_list) if lipid not in lipid_list]
        matrices = []
        for _ in self.u.trajectory[start:final:step]:

            voronoi_dict = self.voronoi_properties(layer = layer)

            areas = np.array(voronoi_dict["areas"])
            #print(voronoi_dict["resnames"])
            #print(type(areas),
            #       voronoi_dict["points"].shape,
            #       areas.shape,
            #       voronoi_dict["resnames"].shape,
            #       voronoi_dict["splay"].shape)
            for lipid in no_present:
                areas[voronoi_dict["resnames"] == lipid] = np.nan

            points = voronoi_dict["points"]


            selection_x = (points[:,0] > xmin - ng*grid_size) & (points[:,0] < xmax + ng*grid_size)
            selection_y = (points[:,1] > ymin - ng*grid_size) & (points[:,1] < ymax + ng*grid_size)
            selection = selection_x & selection_y

            points = points[selection]
            areas = areas[selection]


            matrix,_ = self.map_voronoi(points, areas, nbins, [xmin, xmax, ymin, ymax])
            matrices.append(matrix)
            #print(matrix.shape)




        final_mat = np.nanmean(np.array(matrices), axis = 0)
        final_mat = np.flipud(final_mat)

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
            matrix = self.voronoi_apl(layer = layer, start = i*w_size + start, final =(i+1)*w_size + start , step = step, lipid_list = None, nbins = nbins)
            matrices.append(matrix)
        return matrices


    @staticmethod
    def create_sel_string(lipid_list, sign, mean_z):
        selection_string = f"(((resname {lipid_list[0]} and name {self.working_lip[lipid_list[0]]['head']}) and prop z {sign} {mean_z}))"
        for lipid in lipid_list[1:]:
            selection_string += f" or (((resname {lipid} and name {self.working_lip[lipid]['head']}) and prop z {sign} {mean_z}))"
        return selection_string




    def voronoi_thickness(self, start = 0, final = -1, step = 1, lipid_list = None, nbins = 180, edges = None):
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

        if edges is None:
            xmin = self.v_min
            xmax = self.v_max
            ymin = self.v_min
            ymax = self.v_max
            edges = self.edges
        else:
            xmin = edges[0]
            xmax = edges[1]
            ymin = edges[2]
            ymax = edges[3]

        no_present = [lipid for lipid in list(self.lipid_list) if lipid not in lipid_list]
        matrices = []
        for _ in self.u.trajectory[start:final:step]:





            voronoi_dict = self.voronoi_properties(layer = "top",
                                            lipid_list=self.lipid_list)
            heights = voronoi_dict["heights"]
            for lipid in no_present:
                heights[voronoi_dict["resnames"] == lipid] = np.nan




            matrix_top,_ = self.map_voronoi(voronoi_dict["points"],
                                         heights,
                                         nbins,
                                         [xmin, xmax, ymin, ymax],

                                         )

            voronoi_dict = self.voronoi_properties(layer = "bot",
                                            lipid_list=self.lipid_list)

            heights = voronoi_dict["heights"]
            for lipid in no_present:
                heights[voronoi_dict["resnames"] == lipid] = np.nan
            matrix_bot,_ = self.map_voronoi(voronoi_dict["points"],
                                         voronoi_dict["heights"],
                                         nbins,
                                         [xmin, xmax, ymin, ymax],

                                         )

            matrix_thickness = matrix_top + matrix_bot


            matrices.append(matrix_thickness)

        final_mat = np.nanmean(np.array(matrices), axis = 0)
        final_mat = np.flipud(final_mat)
        return final_mat, edges


    def voronoi_height(self, layer = "top", start = 0, final = -1, step = 1, lipid_list = None, nbins = 180, edges = None):
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


        # Check which lipids are not requested
        no_present = [lipid for lipid in list(self.lipid_list) if lipid not in lipid_list]

        #print(no_present, lipid_list, self.lipid_list)
        if edges is None:
            xmin = self.v_min
            xmax = self.v_max
            ymin = self.v_min
            ymax = self.v_max
            edges = [xmin,xmax,ymin,ymax]
        else:
            xmin = edges[0]
            xmax = edges[1]
            ymin = edges[2]
            ymax = edges[3]


        matrices = []



        for _ in self.u.trajectory[start:final:step]:

            voronoi_dict = self.voronoi_properties(layer = layer)

            heights = voronoi_dict["heights"]
            for lipid in no_present:
                heights[voronoi_dict["resnames"] == lipid] = np.nan


            matrix_height,_ = self.map_voronoi(voronoi_dict["points"],
                                         heights,
                                         nbins,
                                         [xmin, xmax, ymin, ymax],

                                         )




            matrices.append(matrix_height)

        final_mat = np.nanmean(np.array(matrices), axis = 0)
        final_mat = np.flipud(final_mat)
        return final_mat, edges


    def voronoi_splay(self, layer = "top", start = 0, final = -1, step = 1, lipid_list = None, nbins = 180, edges = None):
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

        if edges is None:
            xmin = self.v_min
            xmax = self.v_max
            ymin = self.v_min
            ymax = self.v_max
            edges = self.edges
        else:
            xmin = edges[0]
            xmax = edges[1]
            ymin = edges[2]
            ymax = edges[3]

        self.guess_last_cs()

        no_present = [lipid for lipid in list(self.lipid_list) if lipid not in lipid_list]
        matrices = []



        for _ in self.u.trajectory[start:final:step]:

            voronoi_dict = self.voronoi_properties(layer = layer, splay=True)


            splay_vect = voronoi_dict["splay"]
            for lipid in no_present:
                splay_vect[voronoi_dict["resnames"] == lipid] = np.nan


            matrix_height,_ = self.map_voronoi(voronoi_dict["points"],
                                         splay_vect,
                                         nbins,
                                         [xmin, xmax, ymin, ymax],

                                         )




            matrices.append(matrix_height)

        final_mat = np.nanmean(np.array(matrices), axis = 0)
        final_mat = np.flipud(final_mat)
        return final_mat, edges

    @staticmethod
    def build_resname_atom(resnames, atomsnames):
        resnames = list(resnames)
        string = f"( (resname {resnames[0]}  and name {atomsnames[0]} ) "

        for i, resname in enumerate(resnames[1:]):
            string += f" or (resname {resname}  and name {atomsnames[ i + 1]}) "

        string +=  " ) "
        return string









