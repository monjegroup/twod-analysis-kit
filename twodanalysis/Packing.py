import MDAnalysis as mda
import numpy as np
import pandas as pd
from twodanalysis import MembProp

class PackingDefects(MembProp):
    def __init__(self,
                 universe,
                lipid_list = None,
                verbose = False,
                edges = None,
                nbins = 100,
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
        self.nbins = nbins

        self.nx_polarity()

        self.build_polarity_dict()

        self.radii_dict = {"H": 1.1,
                            "N": 1.55,
                            "C": 1.7,
                            "P": 1.8,
                            "O": 1.52,
                            "S" : 1.8,
                            "K" : 1.2,
                            "CL" : 1.81,
                            "CS" : 1.7
                            }

        self.add_radii(self.radii_dict)


    def packing_defects(self,
                        layer = 'top',
                        nbins = None,
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
            Number of bins of the xy grid, by default None
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

        if nbins is None:
            nbins = self.nbins

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


        # Extend the grid 5 Amstrong to make sure all the packing defects are correctly taken

        xmin_ex = xmin - n_aug * grid_size
        xmax_ex = xmax + n_aug * grid_size
        ymin_ex = ymin - n_aug * grid_size
        ymax_ex = ymax + n_aug * grid_size

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
        defects = np.rot90(defects)

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

            #if names[i] in self.non_polar_dict[lipid]:
            #    small_matrix = small_matrix * 0.0001
            small_matrix = small_matrix * self.polarity_dict[lipid][names[i]]
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

                connections = []


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


    def build_polarity_dict(self):

        polarity_dict = {}
        for key in self.polar_dict.keys():
            polarity_dict[key] = {atom : 1 for atom in self.polar_dict[key]} | {atom : 0.001 for atom in self.non_polar_dict[key]}

        self.polarity_dict = polarity_dict
