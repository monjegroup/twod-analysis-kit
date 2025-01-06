from twodanalysis import Voronoi2D
import numpy as np
import matplotlib.pyplot as plt
class Curvature2D(Voronoi2D):
    def __init__(self,
                universe,
            lipid_list = None,
            verbose = False,
            edges = None,
            nbins = None
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
        self.nbins = 180


    def GaussianCurvature(self, layer = "top", start = 0, final = -1, step = 1, lipid_list = None, nbins = 180, edges = None):
        heights, edges= self.voronoi_height(layer = layer,
                                             start = start,
                                             final = final,
                                             step = step,
                                             lipid_list = lipid_list,
                                             nbins = nbins,
                                             edges = edges)





        dx = edges[1]- edges[0]
        dy = edges[3]- edges[2]


        f_x = np.gradient(heights, dx, axis=1)
        f_y = np.gradient(heights, dy, axis=0)
        f_xx = np.gradient(f_x, dx, axis=1)
        f_yy = np.gradient(f_y, dy, axis=0)
        f_xy = np.gradient(f_x, dy, axis=0)

        K = (f_xx * f_yy - f_xy**2) / (1 + f_x**2 + f_y**2)**2

        return K, edges

    def MeanCurvature(self, layer = "top", start = 0, final = -1, step = 1, lipid_list = None, nbins = 180, edges = None):
        heights, edges= self.voronoi_height(layer = layer,
                                             start = start,
                                             final = final,
                                             step = step,
                                             lipid_list = lipid_list,
                                             nbins = nbins,
                                             edges = edges)
        plt.imshow(heights, extent = edges,cmap = "Spectral")
        plt.xlabel("$x [\AA$]", fontsize = 20)
        plt.ylabel("$y [\AA$]", fontsize = 20)
        plt.colorbar()
        plt.title("height")
        plt.savefig("cumulative_curv.png", dpi = 1000)
        plt.show()

        dx = edges[1]- edges[0]
        dy = edges[3]- edges[2]


        f_x = np.gradient(heights, dx, axis=1)
        f_y = np.gradient(heights, dy, axis=0)
        f_xx = np.gradient(f_x, dx, axis=1)
        f_yy = np.gradient(f_y, dy, axis=0)
        f_xy = np.gradient(f_x, dy, axis=0)

        H = (f_xx * f_yy - f_xy**2) / (1 + f_x**2 + f_y**2)**2

        return H, edges
