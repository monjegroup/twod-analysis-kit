"""
twod_analysis
=============


Classes
-------

.. autoclass:: BioPolymer2D
    :members:
"""

import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import matplotlib as mpl
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.exceptions import SelectionError
import sys
class BioPolymer2D:
    def __init__(self, obj):
        """Initializes the class with either an MDAnalysis Universe or AtomGroup.

        Parameters
        ----------
        obj : (Universe or AtomGroup)

        Raises
        ------
        TypeError
           Error if the is not being initialized with MDAnalysis.Universe or MDAnalysis.AtomGroup
        """
        if isinstance(obj, mda.Universe):
            self.universe = obj
            self.atom_group = obj.atoms  # Select all atoms
        elif isinstance(obj, mda.core.groups.AtomGroup):
            self.universe = obj.universe
            self.atom_group = obj
        else:
            raise TypeError("Input must be an MDAnalysis Universe or AtomGroup")
        
        self.startT=self.universe.trajectory[0].time*0.001
        self.endT=self.universe.trajectory[-1].time*0.001
        self.stepT=self.universe.trajectory.dt*0.001
        self.startF=int(self.startT/self.stepT) 
        self.endF=int(self.endT/self.stepT)
        self.stepF=int(self.stepT/self.stepT)
        self.times=np.arange(self.startT,self.endT,self.stepT)
        self.frames=np.arange(self.startF,self.endF)
        self.pos=None
        self.com=None
        self.system_name=None
        self.kdeanalysis = lambda : None  # Create an empty object-like container
        self.kdeanalysis.paths = None
        self.kdeanalysis.kde = None
        self.hbonds=None
    def __repr__(self):
        return f"<{self.__class__.__name__} with {len(self.atom_group)} atoms>"
    
    def INFO(self):
        """Prints general information of the Universe and/or Atomgroup 
        """
        # Universe-level information
        print("=== UNIVERSE INFO ===")
        print("  N atoms:", len(self.universe.atoms))
        print("  N residues:", len(self.universe.residues))
        print("  N segments:", len(self.universe.segments))
        print(f"  Time : {self.universe.trajectory[0].time/1000}-{self.universe.trajectory[-1].time/1000}ns dt={self.universe.trajectory.dt/1000}ns")
        print(f"  N frames : {len(self.universe.trajectory[0].frames)}")
        # AtomGroup-specific information (only if a subset is selected)
        if len(self.atom_group) < len(self.universe.atoms):
            print("=== SELECTION INFO ===")
            print("  N selected atoms:", len(self.atom_group))
            print("  N selected residues:", len(self.atom_group.residues))
            print("  N selected segments:", len(self.atom_group.segments))

    def getPositions(self,pos_type='COM', inplace=True, select=None):
        """Computes positions of selection from self.startT to self.endT with self.stepT steps of frames. 
        By default, these parameters are set to compute over the whole trajectory.

        Parameters
        ----------
        pos_type : str, optional
            Computes the positions of "all" atoms of the object or the "COM" (center of mass) of residues, by default 'COM'.
        inplace : bool, optional
           If True, position values are assigned to the self.pos attribute and None is returned. If False, positions are returned, by default True
        select : None or str, optional
             If None, all atoms in the Atom group are computed. Otherwise, it is a string selection analogue to MDAnalysis format. Selection must be a set of atoms of the Atom group.  Defaults to None., by default None

        Returns
        -------
        None or np.ndarray
            None if inplace=True, numpy array if inplace=False with the positions of the center of mass of residues (if pos_type="COM")or positions of all atoms (pos_type="all")
        """

        print('Getting positions from frame',self.startF, 'to', self.endF,'with steps of',self.stepF)

        prot=self.atom_group
        if select:
            prot=self.universe.select_atoms(select)
        pos=np.zeros((int((self.endF-self.startF)/self.stepF),len(prot.residues),4)) 
        if pos_type=='all':
            pos=np.zeros((int((self.endF-self.startF)/self.stepF),len(prot.atoms.positions),4)) 
   
        j=0
        for ts in self.universe.trajectory[self.startF:self.endF:self.stepF]:
            if pos_type=='COM':
                pos[j,:,0]=ts.time/1000
                pos[j,:,1:]=[r.atoms.center_of_mass() for r in prot.residues]
                # print('Getting COMs..')
            if pos_type=='all':
                pos[j,:,0]=ts.time/1000
                pos[j,:,1:]= prot.atoms.positions
            j+=1
        if inplace:
            self.pos=np.array(pos)
            return None
        else:
            return np.array(pos)

    def getCOMs(self, inplace=True, select=None):
        """Computes positions of selection from self.startT to self.endT with self.stepT steps of frames. 
        By default, these parameters are set to compute over the whole trajectory.

        Parameters
        ----------
        inplace : bool, optional
            If True, position values are assigned to the self.com attribute and None is returned. If False, center of mass at each frame are returned. By default True
        select : None or str, optional
            If None, all atoms in the Atom group are computed. Otherwise, it is a string selection analogue to MDAnalysis format. Selection must be a set of atoms of the Atom group. By default None

        Returns
        -------
        None or np.ndarray
           None if inplace=True, numpy array if inplace=False with the center of mass of the AtomGroup.
        """

        print('Getting center of masses from frame',self.startF, 'to', self.endF,'with steps of',self.stepF)

        prot=self.atom_group
        if select:
            prot=self.universe.select_atoms(select)

        com=np.zeros((int((self.endF-self.startF)/self.stepF),4)) 
   
        j=0
        for ts in self.universe.trajectory[self.startF:self.endF:self.stepF]:

            com[j,0]=ts.time/1000
            com[j,1:]= prot.atoms.center_of_mass()
            j+=1
        if inplace:
            self.com=np.array(com)
            return None
        else:
            return np.array(com)
        

    # def FilterMinFrames(self, zlim,Nframes,control_plots=False):
    #     """
    #     Selects a set of Nframes in which the AtomGroup is closer to the surface and bellow a zlim threshold distance to the surface. 

    #     Args:
    #         zlim (float): Distance (in angstroms) threshold limit in which the AtomGroup is considered adsorped to the surface. 
    #         Nframes (int): Nframes closest to the surface within the frames where the AtomGroup is < zlim. 
    #         control_plots (bool, optional): If control plots are to be shown. Defaults to False.

    #     Returns:
    #         np.ndarray (Nframes,Nresidues or Natoms,4 <t,x,y,z>): Numpy array with the AtomGroup positions in self.pos that are below zlim and closest to the surface.
    #     """

    #     pos=self.pos
    #     ##Take the mean of all selected residues
    #     mean_z_top=pos[:,:,3].mean(axis=1)
    #     print(mean_z_top.shape)
    #     # Select frames where the mean of selected residues is < zlim
    #     # to garanty that same frames are used for all residues. 
    #     zMask= mean_z_top < zlim
    #     pos_masked=pos[zMask]
    #     print('There are', len(pos_masked),' frames < %i A in Z'%zlim)
    #     print('Taking', Nframes, 'closest frames to surface...')
    #     pos_masked=pos_masked[np.argsort(pos_masked[:,:,3].mean(axis=1),axis=0)][:Nframes]
    #     # print(pos_masked[:5:,0,0],'pos_masked shuffled')

    #     if control_plots:
    #         ires=[0,1] ## If want to change defaut residue to make control plots, change here
    #         print(f"Doing control plot with residue {self.atom_group.residues[ires]}")
    #         plt.plot(pos[:,ires,0],pos[:,ires,3],)
    #         plt.plot(pos_masked[:,ires,0],pos_masked[:,ires,3],'.',ms=2)
    #         plt.show() 
    #     return pos_masked
    
    @staticmethod
    def FilterMinFrames(pos, zlim,Nframes,control_plots=False):
        """Selects a set of Nframes in which the AtomGroup is closer to the surface and bellow a zlim threshold distance to the surface. 


        Parameters
        ----------
        pos : list or np.ndarray (TotalFrames,Nresidues or Natoms,4 <t,x,y,z>)
            Positions over time to be filtered.
        zlim : float
            Distance (in angstroms) threshold limit in which the AtomGroup is considered adsorped to the surface. 
        Nframes : int
            Nframes closest to the surface within the frames where the AtomGroup is < zlim.
        control_plots : bool, optional
            If control plots are to be shown, by default False

        Returns
        -------
        np.ndarray
            Filtered positions (Nframes,Nresidues or Natoms,4 <t,x,y,z>)
        """
        # pos=self.pos
        ##Take the mean of all selected residues
        mean_z_top=pos[:,:,3].mean(axis=1)
        print(mean_z_top.shape)
        # Select frames where the mean of selected residues is < zlim
        # to garanty that same frames are used for all residues. 
        zMask= mean_z_top < zlim
        pos_masked=pos[zMask]
        print('There are', len(pos_masked),' frames < %i A in Z'%zlim)
        print('Taking', Nframes, 'closest frames to surface...')
        pos_masked=pos_masked[np.argsort(pos_masked[:,:,3].mean(axis=1),axis=0)][:Nframes]
        # print(pos_masked[:5:,0,0],'pos_masked shuffled')

        if control_plots:
            ires=[0,1] ## If want to change defaut residue to make control plots, change here
            # print(f"Doing control plot with residue {self.atom_group.residues[ires]}")
            plt.plot(pos[:,ires,0],pos[:,ires,3],)
            plt.plot(pos_masked[:,ires,0],pos_masked[:,ires,3],'.',ms=2)
            plt.show() 
        return pos_masked
    
    def PolarAnalysis(self,select_res,Nframes,max_hist=None,sort='max',plot=False,control_plots=False, zlim=14,Nbins=1000,resolution=5):
        """Makes a Polar Histogram of the positions of the center of mass of select_res residues considering Nframes closest to the surface within the < zlim threshold. self.pos attribute is used to compute the center of mass of the AtomGroup, which will be the referential center of the histograms. 
        The colors of the histogram are ordered according to sort parameter. 

        Parameters
        ----------
        select_res : str
            MDAnalysis string selection of the residues to which compute the histograms.
        Nframes : int
            Nframes closest to the surface within the frames where the AtomGroup is < zlim. 
        max_hist : None or float, optional
            Value to normalize the histograms. If None, highest histogram value is used to normalize the histograms., by default None
        sort : str, list, ndarray or None, optional
            How to sort the histograms in the legend and the coloring. If "max", histograms will be ploted descending from the residue with highest peak in its histogram to the flattest peak. If sort is a list or ndarray, sort is the index positions of the residues to which reorder the histograms. If None, they well be plotted by MDAnalysis default sort (ascending Resid values)., by default 'max'
        plot : bool, optional
            Show the polar plot (True) or only return the data (False), by default False
        control_plots : bool, optional
            Show control plots of the different steps of the polar analysis calculation., by default False
        zlim : float, optional
            Distance (in angstroms) threshold limit in which the AtomGroup is considered adsorped to the surface, by default 14
        Nbins : int, optional
             How many bins to use in the histograms., by default 1000
        resolution : float, optional
            One the position vectors of each residue are normalized, resolution is the value too which the normalized positions are multiplied. An increase in this value will make higher peaks in the histograms since position vector are further away. Reducing this value represents an increase of resolution of the histogram. By default 5
        
        Returns
        -------
        list,np.ndarray
            Histogram data (Nresidues, 2 <X_bin_values,Y_bin_values>, Nbins),
            Positions in order in which they were plotted (Nframes,Nresidues or Natoms,4 <t,x,y,z>)
        """

        # colors=['C%s'%(c+7) for c in range(10)]
        prot_center=self.atom_group
        to_center=prot_center.atoms.center_of_mass()

        fromF=self.startF
        endF=self.endF
        print(f"Computing Polar Analysis from frame {fromF} (t={self.startT}ns) to {endF} (t={self.endT}ns) ")
        pos=self.getPositions(select=select_res,inplace=False)
        prot=self.universe.select_atoms(select_res)

        print(prot.residues)
        print(pos.shape)
        # sys.exit()

        pos_centered=pos-np.array([0,to_center[0],to_center[1],0])
        pos_selected=BioPolymer2D.FilterMinFrames(pos_centered,zlim,Nframes,control_plots=control_plots)

        print(pos_selected.shape)

        
        # center_arr=np.array(center_arr)
        if control_plots==True:
            for r in range(len(pos_selected[0])):
                plt.plot(pos_selected[:,r,1],pos_selected[:,r,2],'o')

            plt.plot(0,0,'o',c='k')

            plt.show()   

        # sys.exit()

        MeanPos=np.mean(pos_selected, axis=0) 
        print(MeanPos.shape, 'MeanPos')
        sign_arr=np.where(MeanPos[:,2]<0)#[0]
        angles = np.arccos(MeanPos[:,1]/np.linalg.norm(MeanPos[:,[1,2]], axis=1))

        # angles[sign_arr]=-angles[sign_arr]
        # angles=2*np.pi-angles#[sign_arr]
        angles[sign_arr]=2*np.pi-angles[sign_arr]

        # print(angles)
        # pos_selected[:,:,0]= pos_selected[:,:,0]-MeanPos[:,0]
        # pos_selected[:,:,1]= pos_selected[:,:,1]-MeanPos[:,1]
        pos_selected[:,:,1]+=resolution*np.cos(angles)
        pos_selected[:,:,2]+=resolution*np.sin(angles)
        norms=np.linalg.norm(pos_selected[:,:,[1,2]], axis=2)
        # print(norms.shape)
        if control_plots==True:
            for r in range(len(pos_selected[0])):
                plt.plot(pos_selected[:,r,1],pos_selected[:,r,2],'o')
            plt.plot(0,0,'o',c='k')
            plt.show()
            
        pos_selected[:,:,1]=pos_selected[:,:,1]/norms
        pos_selected[:,:,2]=pos_selected[:,:,2]/norms
        
        if control_plots==True:
            for r in range(len(pos_selected[0])):
                plt.plot(pos_selected[:,r,1],pos_selected[:,r,2],'o')
            plt.plot(0,0,'o',c='k')
            plt.show()

        N = Nbins
        bottom = 1
        max_height = 1
        theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
        # time = np.linspace(0.0, 1501,1501)
        max_arr=[]
        hist_arr=[]
    #     t=np.arange(Univs[0].trajectory[fromF].time,Univs[0].trajectory[endF].time,Univs[0].trajectory[0].dt)/1000
        for r in range(len(pos_selected[0])):
            res_pos_selected=pos_selected[:,r]
            sign_arr=np.where(res_pos_selected[:,2]<0)#[0]
            radii = np.arccos(res_pos_selected[:,1]/np.linalg.norm(res_pos_selected[:,[1,2]], axis=1))
    #         # sign_arr=np.where(pos_selected[:,i,1]<0)#[0]
    #         # radii = np.arccos(pos_selected[:,i,0]/np.linalg.norm(pos_selected[:,i,[0,1]], axis=1))
            radii[sign_arr]=2*np.pi-radii[sign_arr]
            hist=np.histogram(radii, bins=theta)
            max_arr.append(np.max(hist[0]))
            hist_arr.append([hist[1][1:], hist[0]])
        width = (2*np.pi) / N
        hist_arr=np.array(hist_arr)

        if not max_hist:
            max_hist=np.max(max_arr)
        else:
            max_hist=max_hist

        if sort=='max':
            sort_i=np.argsort(max_arr)[::-1]
            ordered_selected_pos=pos_selected[:,sort_i]
        elif isinstance(sort,(list,np.ndarray)):
            sort_i=sort
            ordered_selected_pos=pos_selected[:,sort_i]
        else:
            ordered_selected_pos=pos_selected
            sort_i=np.arange(ordered_selected_pos.shape[1])
        print(pos_selected.shape, ordered_selected_pos.shape)
        print(pos_selected[:,:3].shape, ordered_selected_pos[:,:3].shape)
        if plot==True:
            norm_max_hist=1/max_hist
            ax = plt.subplot(111, polar=True)
            for i in range(ordered_selected_pos.shape[1]):
                fig_label='%i-%s'%(prot.residues.resids[sort_i[i]],prot.residues.resnames[sort_i[i]])

                with mpl.rc_context():
                    # mpl.style.use('classic')
                    bars = ax.bar(hist_arr[sort_i[i],0], hist_arr[sort_i[i],1]*max_height*norm_max_hist,
                            width=width,label='%s'%fig_label,
                            bottom=bottom,#color=colors[i],
                            alpha=0.8
                        )
                ax.set_yticklabels([])
            plt.legend(loc='center', bbox_to_anchor=(0.5, 0.5),framealpha=1)
            plt.ylim((0,2))
            # plt.title('Polar Histogram RBD-%s %s variant'%('PBL', plabels[ivar]))
            plt.tight_layout()
            plt.show()
        #     ax.set_theta_zero_location('N')
        return hist_arr,ordered_selected_pos
    
    ######## Radii of Gyration 2D Analysis ###################

    def computeRG2D(self, masses, total_mass=None):
        """Computes parallel, perpendicular and 3D radius of gyration in 1 frame. 

        .. math:: R_{\\textrm{g}\parallel}= \sqrt{ \\frac{1}{m_T}\sum_{i} m_{i}\left[ (x_i-x_{\\textrm{CM}})^2+(y_i-y_{\\text{CM}})^2\\right]}
        
        .. math:: R_{\\textrm{g}\perp} = \sqrt{\\frac{1}{m_T}\sum_{i} m_{i} (z_i-z_{\\text{CM}})^2,}

        where :math:`{\\bf R}_{\\textrm{CM}}=(x_{\\textrm{CM}}`, :math:`y_{\\textrm{CM}}`, :math:`z_{\\textrm{CM}})` is the position of the center of mass, :math:`m_{i}` the mass of each residue and :math:`m_T` the total mass of the residues.

        
        Parameters
        ----------
        masses : np.ndarray (Natoms)
            Mass vallues of each atom.
        total_mass : float, optional
            Sum of these masses. By default None

        Returns
        -------
        np.ndarray
            3D, perpendicular and parallel radius of gyration values (3, Natoms)
        """
        # coordinates change for each frame
        coordinates = self.atom_group.positions
        center_of_mass = self.atom_group.center_of_mass()

        # get squared distance from center
        ri_sq = (coordinates-center_of_mass)**2
        # sum the unweighted positions
        sq = np.sum(ri_sq, axis=1)
        sq_par = np.sum(ri_sq[:,[0,1]], axis=1) # sum over x and y, parallel Rg
        #print(sq_par.shape)
        sq_perp = ri_sq[:,2] # sum over z for perpendicular Rg

        # make into array
        sq_rs = np.array([sq, sq_perp, sq_par])

        # weight positions
        rog_sq = np.sum(masses*sq_rs, axis=1)/total_mass
        # square root and return
        return np.sqrt(rog_sq)
    
    def getRgs2D(self, plot=False):
        """Computes the radii of gyration over the trajectory using computeRG2D method. 

        Parameters
        ----------
        plot : bool, optional
            If True, it plots the raw Rg data values in time, by default False

        Returns
        -------
        np.ndarray
            3D, perpendicular and parallel radius of gyration values at each frame (self.endF-self.startF frames,3, Natoms)
        """
        colors=['tab:blue','tab:orange','tab:green']
        rg_arr=np.zeros((len(self.pos),4))
        i=0
        masses=self.atom_group.atoms.masses
        total_mass=np.sum(masses)

        for ts in self.universe.trajectory[self.startF:self.endF:self.stepF]:
            rg_arr[i,0]=ts.time/1000
            rgs=self.computeRG2D(masses, total_mass)
            rg_arr[i,1:]=rgs
            i+=1
        legend_names=['3D','Perp','Parallel']
        if plot:
            plt.plot(rg_arr[:,0], rg_arr[:,1], label=legend_names[0],color=colors[0])
            plt.plot(rg_arr[:,0], rg_arr[:,2], label=legend_names[1],color=colors[1])
            plt.plot(rg_arr[:,0], rg_arr[:,3], label=legend_names[2],color=colors[2])
            plt.legend()
            plt.xlabel('Time (ns)')    
            plt.ylabel('Radius of gyration (Angs)')    
            plt.show()
        return rg_arr    

    def RgPerpvsRgsPar(self,rgs,color, marker='s',plot=True,show=False):
        """Generates :math:`R_{g\perp}` vs. :math:`R_{g\parallel}` plots. Also, returns the :math:`\langle R_{g\perp}^2 \\rangle /\langle R_{g\parallel}^2 \\rangle` ratio

        Parameters
        ----------
        rgs : np.ndarray 
            Radii of gyration data. It must have Rg data as [Rg 3D, Rg parallel, Rg perpendicular], analogue to getRgs2D functions
        color : str
            Color used to plot points. Color names use the same of those of Matplotlib package.
        marker : str, optional
            Marker used to plot. Marker names use the same of those of Matplotlib package. , by default 's'
        plot : bool, optional
            If False, only the :math:`\langle R_{g\perp}^2 \\rangle /\langle R_{g\parallel}^2 \\rangle` ir returned with out make plot, by default True
        show : bool, optional
            If True, matplotlib.pyplot.show() is executed. This is left optional in case further characterization of plot is desired. Also, this enables showing multple data in the same figure. By default False.

        Returns
        -------
        float
            :math:`\langle R_{g\perp}^2 \\rangle /\langle R_{g\parallel}^2 \\rangle`
        """
        data=rgs[:,2:]
        print(data.shape)
        rg_ratio=(data[:,0]**2).mean()/(data[:,1]**2).mean()
        if plot:
            plt.plot(data[:,1],data[:,0],'o',markersize=1,c=color)
            label='%s (%.3f)'%(self.system_name,rg_ratio)
            plt.plot(data[:,1].mean(),data[:,0].mean(),marker,markersize=10, label=label,color='k')
            plt.legend(title=r'Syst ($\langle Rg_\perp^2\rangle /\langle Rg_\parallel^2 \rangle$)')
            plt.ylabel(r'$Rg_\parallel$ (angs)')
            plt.xlabel(r'$Rg_\perp$ (angs)')
            if show:
                plt.show()
        return rg_ratio

    ############# Compute Contour Area #################
    @staticmethod
    def ListPathsInLevel(kde_plot,contour_level,plot_paths=False):
        """Lists vertices of a path in a contour level of the KDE plot. 

        Parameters
        ----------
        kde_plot : kde_plot
            Results of seaborn.kde_plot function. Once computed getKDEAnalysis, self.kde stores this input.
        contour_level : int
            Contour Level in which to enlist paths or vertices.
        plot_paths : bool, optional
            If True,all the paths of a given contour level are plotted, by default False
        """
        
        # Extract the contour levels
        contour_collections = kde_plot.collections
        outermost_contour = contour_collections[-1]
        
        # Extract the contour paths
        contour_paths = outermost_contour.get_paths()
        # inner_contour_paths = innermost_contour.get_paths()
        # print(contour_paths)
        
        # Extract the vertices of the outermost contour
        vertices = contour_paths[contour_level].vertices
        codes = contour_paths[contour_level].codes
        
        # inner_vertices = contour_paths[-1].vertices
    #    print(vertices.shape)
    #    print(codes.shape)
        
        # Find how many borders
        paths=[]
        i=1
        # print(codes)
        while i < len(vertices):
            one_path=[]
            first_i_of_path=i
            while codes[i] == 2:
                # plt.scatter(vertices[i,0],vertices[i,1],c='k')
                one_path.append([vertices[i,0],vertices[i,1]])
                i+=1
                if i == len(vertices)-1:
                    break
            one_path.append([vertices[first_i_of_path,0],vertices[first_i_of_path,1]])
            one_path=np.array(one_path)
            paths.append(one_path)
            if plot_paths:
                plt.plot(one_path[:,0],one_path[:,1],color='orange')    
            i+=2
        return paths
        
    def getKDEAnalysis(self,zlim,Nframes,inplace=True,control_plots=False):
        """Computes KDE Contours using seaborn.kde_plot() function and extracts the paths of each contour level. The output of the seaborn.kde_plot() is stored in self.kdeanalysis.kde, and the paths of each contour level is stored in self.kdeanalysis.paths if inplace=True.

        Parameters
        ----------
        zlim : float
            zlim of BioPolymer2D.FilterMinFrames(). Only use frames under a zlim threshold, to avoid using frames with desorbed molecule.
        Nframes : int
            Nframes of BioPolymer2D.FilterMinFrames(). To ensure to have a controled number of frames under zlim threshold.
        inplace : bool, optional
            If True, stores the paths of al contour levels in self.kdeanalysis.paths. Otherwise, it only returns it. By default True
        control_plots : bool, optional
            Make control plots, by default False

        Returns
        -------
        list
            List of all paths in all the contour levels.
        """
        pos_selected=BioPolymer2D.FilterMinFrames(self.pos,zlim,Nframes,control_plots=control_plots)
        ## Concatenate positions of all residues
        print(pos_selected.shape)
        pos_selected_reshape=np.reshape(pos_selected,(pos_selected.shape[0]*pos_selected.shape[1],pos_selected.shape[2]))
        print(pos_selected_reshape.shape)
        fig,ax=plt.subplots()  
        df0=pd.DataFrame(pos_selected_reshape, columns=['t','x','y','z'])
        # kde_plot=sns.kdeplot(df0, x='x',y='y', color = 'black',alpha=1,fill=True)
        kde_plot = sns.kdeplot(df0, x='x', y='y', fill=True, cmap="Purples", color='black')
        kde_plot=sns.kdeplot(df0, x='x',y='y', color = 'black',alpha=0.5)
        self.kdeanalysis.kde=kde_plot

        paths_arr=[]
        Nlvls=len(kde_plot.collections[-1].get_paths())
        print(f"There are {Nlvls} levels in the KDE.")
        for lvl in range(Nlvls):
            paths=BioPolymer2D.ListPathsInLevel(kde_plot,lvl,plot_paths=control_plots)
            if not paths:
                continue
            # print(np.shape(paths[0]))
            # print(np.shape(paths[1]))
            paths_arr.append(paths)
            # plt.pause(2)
            if control_plots:
                plt.pause(2)
        plt.close()
        if inplace:
            self.kdeanalysis.paths=paths_arr
        return paths_arr

    # def plotPathsInLevel(self, contour_lvl,color='k',alpha=0.3,show=False):
    #     paths_in_lvl=self.kdeanalysis.paths[contour_lvl]
    #     for p in range(len(paths_in_lvl)):
    #         x_val,y_val=paths_in_lvl[p].T
    #         plt.plot(x_val,y_val,color=color, alpha=alpha)
    #     if show:
    #         plt.show()


    def getAreas(self,contour_lvl,getTotal=False):
        """Computes the area of each path a given contour level. Negative values of area are holes in the contour level. If getTotal=True, computes the area of the whole contour level.

        Parameters
        ----------
        contour_lvl : int
            Contour level to compute the area.
        getTotal : bool, optional
            If False, gives the area of each path of the contour level. If getTotal=True, sums up the contribution of each path in the contour level returning the area of the whole contour level. By default False

        Returns
        -------
        list or float
            A list with the area of each path in the contour level, or the total area of the contour level (if getTotal=True)
        """
        if not self.kdeanalysis.paths:
            print("Must compute contour paths first. Please compute getKDEAnalysis method first.")
            return 

        Areas=[]    
        paths_in_lvl=self.kdeanalysis.paths[contour_lvl]
        for i_paths in range(len(paths_in_lvl)):
            x_values,y_values=paths_in_lvl[i_paths].T
            area_outline = simpson(y_values[::-1], x=x_values[::-1])
            print("Area of the outline %s:"%i_paths,np.abs(area_outline))
            Areas.append(np.abs(area_outline))        
        if getTotal:
            return np.sum(Areas)
        else:
            return Areas
    def getHbonds(self,selection1,selection2, update_selections=True,trj_plot=False, inplace=True ):
        """Computes H-bonds between to selection1 and selection2 of the trajectory using MDAnalysis.analysis.hydrogenbonds.HydrogenBondAnalysis.

        Parameters
        ----------
        selection1 : str
            First selection to which compute H-bonds
        selection2 : str
            Second selection to which compute H-bonds
        update_selections : bool, optional
            Fills parameter update_selection from MDAnalysis.analysis.hydrogenbonds.HydrogenBondAnalysis. If True, it updates donor and acceptor at each frame. , by default True
        trj_plot : bool, optional
            True to have a look at the Hbonds over the trajectory. By default False
        inplace : bool, optional
            If True, it stores results of H-bond run in self.hbonds. by default True
        
        
        Returns
        -------
        MDAnalysis.analysis.base.Results
            A list with the area of each path in the contour level, or the total area of the contour level (if getTotal=True)
        """

        u=self.universe
        hbonds = HydrogenBondAnalysis(
            universe=u,
#             donors_sel=None,
            between=[selection2, selection1],
            d_a_cutoff=3.0,
            d_h_a_angle_cutoff=150,
            update_selections=update_selections
        )


        selection1_hydrogens_sel = hbonds.guess_hydrogens(selection1)
        selection1_acceptors_sel = hbonds.guess_acceptors(selection1)

        selection2_hydrogens_sel = hbonds.guess_hydrogens(selection2)
        selection2_acceptors_sel = hbonds.guess_acceptors(selection2)
        
        if not selection1_hydrogens_sel:
            raise SelectionError(f"Selection1 has no hydrogen donors.")
        if not selection1_acceptors_sel:
            raise SelectionError(f"Selection1 has no acceptors.")
        if not selection2_hydrogens_sel:
            raise SelectionError(f"Selection2 has no hydrogen donors.")
        if not selection2_acceptors_sel:
            raise SelectionError(f"Selection2 has no acceptors.")


        hbonds.hydrogens_sel = f"({selection1_hydrogens_sel}) or ({selection2_hydrogens_sel})"
        hbonds.acceptors_sel = f"({selection1_acceptors_sel}) or ({selection2_acceptors_sel})"
        hbonds.run(verbose=True)

        if inplace:
            self.hbonds=hbonds.results

        if trj_plot:
            plt.plot(hbonds.times/1000, hbonds.count_by_time(), lw=2, label='%s'%self.system_name)
            plt.xlabel("Time (ns)")
            plt.ylabel(r"$N_{HB}$")
            plt.show()
        return hbonds.results
    
    def HbondsPerResidues(self,sorted=True):
        """Computes the number of H-bonds of each residue during the simulation. self.getHbonds(inplace=True) must be computed prior to the use this function.

        Parameters
        ----------
        sorted : bool, optional
            If True, returns data sorted by number of H-bonds. By default True

        Returns
        -------
        pandas.DataFrame
            DataFrame showing all the residues with H-bonds.
        """
        result=np.array(self.hbonds.hbonds[:,[2,3]], dtype=int)
        resids=np.array([self.universe.atoms[result[:,0]].resids,self.universe.atoms[result[:,1]].resids])
        df=pd.DataFrame(data=resids.T, columns=['Hydrogen', 'Acceptors'])
        Acc_resids = df[np.isin(df['Acceptors'],self.atom_group.residues.resids)].pivot_table(columns=['Acceptors'], aggfunc='size')
        H_resids = df[np.isin(df['Acceptors'],self.atom_group.residues.resids)].pivot_table(columns=['Hydrogen'], aggfunc='size')
        # print(H_resids,Acc_resids)
        final_count=Acc_resids.add(H_resids, fill_value=0)
        df_final=pd.DataFrame(np.array([self.universe.residues[final_count.index-1].resids,self.universe.residues[final_count.index-1].resnames,final_count]).T,
                            columns=['ResIDs','ResNames','Count'])

        if sorted:
            return df_final.sort_values('Count', ascending=False)
        else:
            return df_final
    def plotHbondsPerResidues(self, paths_for_contour,top=-1,contour_lvls_to_plot=None, print_table=True): ### Add a Residue Filter option
        """Makes a figure showing the center of mass of the residues with H-bonds. Figure shows a contour plot as a reference of position of the whole molecule. Legend of the Figure shows the percentage of time in which there were Hbonds during the simulation of the plotted residues. 

        Parameters
        ----------
        paths_for_contour : list
            List of paths of all the contour levels.
        top : int, optional
            Residues are plotted ranked by residues with most contact to least. This parameters indicates how many residues to plot of these ranked residues, e.g. top=5 wil plot the 5 residues with most Hbonds during the simulations. By default -1, plots all the residues with H-bonds.
        contour_lvls_to_plot : list, optional
            Contour Levels to show in plot, by default None
        print_table : bool, optional
            Whether or not to print the pandas.DataFrame with the data shown in figure, by default True

        Returns
        -------
        pandas.DataFrame
            Sorted DataFrame from the residues with most contact to the residue with least contacts.
        """

        df=self.HbondsPerResidues(sorted=False)
        max_val=df['Count'].max()
        str_resids=' '.join(np.array(df['ResIDs'],dtype=str))
        print(str_resids)

        pos=self.getPositions(select=f'resid {str_resids}', inplace=False)
        df[['X','Y', 'Z']]=pos[:,:,1:].mean(axis=0)
        sorted_df=df.sort_values('Count', ascending=False)
        if print_table:
            print(sorted_df.iloc[:top])

        if not contour_lvls_to_plot:
            contour_lvls_to_plot=range(len(paths_for_contour))
        for lvl in contour_lvls_to_plot:
            BioPolymer2D.plotPathsInLevel(paths_for_contour,lvl)

        colors = ['C%s' % i for i in range(10)]  # Define color palette
        num_colors = len(colors)
        for i in range(pos[:,:top].shape[1]):
            # Use modular indexing to cycle through colors
            color = colors[i % num_colors]
            norm_val=sorted_df['Count'].iloc[i]/len(self.universe.trajectory) #max_val
            norm_val_plot=sorted_df['Count'].iloc[i]/max_val
            pos=sorted_df[['X','Y','Z']].iloc[i]
            plt.plot(pos['X'],pos['Y'], 'o',color=color,
                    label='%s-%s (%.2f)'%(sorted_df['ResIDs'].iloc[i],
                                        sorted_df['ResNames'].iloc[i],
                                        norm_val*100),)
            plt.scatter(pos['X'],pos['Y'], s=(8*20*norm_val_plot)**2, alpha=.5, color=color)
        plt.xlabel(r'X-axis($\AA$)',)#fontsize=20)
        plt.ylabel(r'Y-axis($\AA$)',)#fontsize=20)
        plt.gca().set_aspect('equal')
        plt.tight_layout()    
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),#prop={'size':22}, 
                    title="ResID-ResName(Hbond %)",)#title_fontsize=20)
        plt.show()
        return sorted_df
    
    @staticmethod   
    def plotPathsInLevel(paths, contour_lvl,color='k',alpha=0.3,show=False):
        """Plots the paths of a given contour level. 

        Parameters
        ----------
        paths : list
            List with all the paths of all the contour levels.
        contour_lvl : int
            Contour level to plot
        color : str, optional
            Color to plot the contour level , by default 'k', corresponding to a black color. 
        alpha : float, optional
            Sets trasparency. Fills parameter alpha of matplolib.pyplot.plot method by default 0.3
        show : bool, optional
            Where to show or not the plot yet with matplolib.pyplot.plot, by default False
        """
        paths_in_lvl=paths[contour_lvl]
        for p in range(len(paths_in_lvl)):
            x_val,y_val=paths_in_lvl[p].T
            plt.plot(x_val,y_val,color=color, alpha=alpha)
        if show:
            plt.show()