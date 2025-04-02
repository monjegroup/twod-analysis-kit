"""
BioPolymer2D
=============


Classes
-------

.. autoclass:: BioPolymer2D
    :members:
    :undoc-members:
    :show-inheritance:
"""
import matplotlib.pyplot as plt

import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import matplotlib as mpl
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.exceptions import SelectionError
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colors
from matplotlib.lines import Line2D
import sys
import logging





class BioPolymer2D:
    """Class Initialization.

    Parameters
    ----------
    obj : AtomGroup or Universe
        Selection to initialize the class. If initialized with a Universe it will take the whole Universe as AtomGroup.
        The MDAnalysis Universe and AtomGroup will be accesible to with the ``self.universe`` and ``self.atom_group`` class attributes.
    surf_selection : str, optional
        String selection in MDAnalysis format to define the surface (recomended). If ``surf_selection`` is given, ``surf_pos`` attribute will save the mean position 
        of the surface over the time. This will particularly important for setting the ``zlim`` variable some of the analysis require. By default None
    biopol_selection : str, optional
        String selection in MDAnalysis format to define the biopolymer as a selection of a Universe or AtomGroup.
        If None, all the AtomGroup or Universe from ``obj`` is used. Default in None.
    start : float, optional
        Starting time or frame from the trajectory to initialize the object. If ``by_frames=True``, it will be considered as the frame number of the trajectory to start 
        the object. If ``by_frames=False``, it will set the time in ns. If None, `start` is frame 0 of the trajectory. By default None
    step : float, optional
        Time or frame steps from the trajectory to initialize the object. If ``by_frames=True``, it will be considered as the step of frames of the trajectory to
        consider. If ``by_frames=False``, it will set the time steps in ns. If None, all the frames of the trajectory are considered. By default None
    end : float, optional
        Final time or frame from the trajectory to initialize the object. If ``by_frames=True``, it will be considered as the final frame number of the trajectory to consider.
        If ``by_frames=False``, it will set the time in ns. If None, `end` is frame 0 of the trajectory. By default None
    by_frame : bool, optional
        Whether or not to consider 'start', 'step' and 'end' inputs as frame values. If False, they are consider as time values in ns.
    surf_axis : str, optional
        Set in which distance is the surface. By default 'z'

    Raises
    ------
    TypeError
        Error raised if class is not initialized with a Universe or AtomGroup
    """
    def __init__(self, obj,surf_selection=None,biopol_selection=None , start=None, step=None,end=None,by_frames=True,surf_axis='z'):

        if isinstance(obj, mda.Universe) and biopol_selection is None:
            self.universe = obj
            self.atom_group = obj.atoms  # Select all atoms
        elif isinstance(obj, mda.core.groups.AtomGroup) and biopol_selection is None:
            self.universe = obj.universe
            self.atom_group = obj
        elif not biopol_selection is None:
            if isinstance(biopol_selection, str):
                self.universe = obj.universe
                self.atom_group = obj.select_atoms(biopol_selection)
            else:
                raise TypeError("`biopol_selection` must be a string and a valid MDAnalysis selection")
        else:
            raise TypeError("Input must be an MDAnalysis Universe or AtomGroup")
        
        # if not biopol_selection is None:
        #     if isinstance(biopol_selection, str): 



        
        #Error or default value for surface direction? 
        # if surf_selection and not surf_axis:
        #     raise TypeError("If a surface selections if assigned, you must indicate the axis of the surface with surf_axis parameter. Options are 'x','y' or 'z')"

        # Initialize trajectory attributes
        self._startT = self.universe.trajectory[0].time * 0.001
        self._endT = self.universe.trajectory[-1].time * 0.001
        self._stepT = self.universe.trajectory.dt * 0.001

        self._startF = self.universe.trajectory[0].frame
        self._endF = self.universe.trajectory[-1].frame
        self._stepF = 1
        # Initialize trajectory attributes
        
        if not start is None:
            if isinstance(start,(float,int)):
                if by_frames is True:
                    if isinstance(start,int):
                        self._startF=start
                        # self._recalculate_frames(triggered_by='frame',to_update='start')
                    else:
                        raise TypeError("`start` must be an integer.")

                elif by_frames is False:
                    self._startT=start
                    # self._recalculate_frames(triggered_by='time',to_update='start')
                else:
                    raise TypeError("`by_frames` must be boolean. Can only be True or False.")
            else:
                raise TypeError("`start` must be a float, or integer if `by_frames` is True.")
        

        if not step is None:
            if isinstance(step,(float,int)):
                if by_frames is True:
                    if isinstance(step,int):
                        self._stepF=step
                        # self._recalculate_frames(triggered_by='frame',to_update='step')
                    else:
                        raise TypeError("`step` must be an integer.")

                elif by_frames is False:
                    self._stepT=step
                    # self._recalculate_frames(triggered_by='time',to_update='step')
                else:
                    raise TypeError("`by_frames` must be boolean. Can only be True or False.")
            else:
                raise TypeError("`step` must be a float, or integer if `by_frames` is True.")
        

        if not end is None:
            if isinstance(end,(float,int)):
                if by_frames is True:
                    if isinstance(end,int):
                        self._endF=end
                    else:
                        raise TypeError("`end` must be an integer.")

                elif by_frames is False:
                    self._endT=end
                else:
                    raise TypeError("`by_frames` must be boolean. Can only be True or False.")
            else:
                raise TypeError("`end` must be a float, or integer if `by_frames` is True.")
        
        self._sim_startT = self.universe.trajectory[0].time * 0.001
        self._sim_dt = self.universe.trajectory.dt * 0.001

        # Call _recalculate_frames only after setting the value
        self._recalculate_frames(triggered_by="frame" if by_frames else "time")
        self.surf_axis=surf_axis
        self.surf_pos = None

        if not surf_selection is None:
            if isinstance(surf_selection, str):
                self.surf_pos=self.getPositions(select=surf_selection,inplace=False,surf_is_zero=False).mean(axis=(0,1)) #If want mean only over atoms set axis=1
                print(' Mean position of the surface is', self.surf_pos)
            else:
                raise TypeError("`surf_selection` must be a string and a valid MDAnalysis selection")



        self.pos = None
        self.com = None
        self.system_name = None
        self.kdeanalysis = lambda: None  # Create an empty object-like container
        self.kdeanalysis.paths = None
        self.kdeanalysis.kde = None
        self.hbonds = None


    def _recalculate_frames(self, triggered_by="time",):
        """Recalculate frame-related attributes based on ``startT``, ``endT``, and ``stepT`` 
        or startF, endF, and stepF."""
        if triggered_by == "time":
            # Update frame-related attributes based on time-related attributes
            self._startF = round(self._startT / self._sim_dt)
            self._endF = round(self._endT / self._sim_dt)
            self._stepF = round(self._stepT / self._sim_dt) 
            # self._stepF = int(self._stepT / self._stepT)  # Should always be 1

        if triggered_by == "frame":
            self._stepT = self._stepF * self._sim_dt
            self._startT = self._startF * self._sim_dt
            self._endT = self._endF * self._sim_dt

        self.times = np.arange(self._startT, self._endT, self._stepT)
        self.frames = np.arange(self._startF, self._endF,self._stepF)

    # Properties for time attributes
    @property
    def startT(self):
        return self._startT

    @startT.setter
    def startT(self, value):
        self._startT = value
        self._recalculate_frames(triggered_by="time")

    @property
    def endT(self):
        return self._endT

    @endT.setter
    def endT(self, value):
        self._endT = value
        self._recalculate_frames(triggered_by="time")

    @property
    def stepT(self):
        return self._stepT

    @stepT.setter
    def stepT(self, value):
        self._stepT = value
        self._recalculate_frames(triggered_by="time")

    # Properties for frame attributes
    @property
    def startF(self):
        return self._startF

    @startF.setter
    def startF(self, value):
        self._startF = value
        self._recalculate_frames(triggered_by="frame")

    @property
    def endF(self):
        return self._endF

    @endF.setter
    def endF(self, value):
        self._endF = value
        self._recalculate_frames(triggered_by="frame")

    @property
    def stepF(self):
        return self._stepF

    @stepF.setter
    def stepF(self, value):
        self._stepF = value
        self._recalculate_frames(triggered_by="frame")

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
        print(f"  Time : {self.startT}-{self.endT}ns dt={self.stepT}ns")
        print(f"  N frames : {self.frames.shape[0]}")
        # AtomGroup-specific information (only if a subset is selected)
        # if len(self.atom_group) < len(self.universe.atoms):
        print("=== SELECTION INFO ===")
        print("  N selected atoms:", len(self.atom_group))
        print("  N selected residues:", len(self.atom_group.residues))
        print("  N selected segments:", len(self.atom_group.segments))

    def getPositions(self,pos_type='COM', surf_is_zero=True, inplace=True, select=None,getselection=False):
        r"""Computes positions of selection from `self.startT` to ``self.endT`` with ``self.stepT`` steps of times.
        By default, these parameters are set to compute over the whole trajectory.

        Parameters
        ----------
        pos_type : str, optional
            Computes the positions of "all" atoms of the object or the "COM" (center of mass) of residues, by default 'COM'.
        surf_is_zero : bool, optional
            Whether to set position coordinate of the BioPolymer respect to the surface, computing :math:`z_{\textrm{biopol}}^i-z_{\textrm{surf}}` for each residue or atom. By default is True.        
        inplace : bool, optional
           If True, position values are assigned to the self.pos attribute and None is returned. If False, positions are returned, by default True
        select : None or str, optional
            If None, all atoms in the Atom group are computed. Otherwise, it is a string selection analogue to MDAnalysis format. Selection must be a set of atoms of the Atom group.  Defaults to None.
            Whether or not to return the selected MDAnalysis AtomGroup as output.
        getselection:
            Whether or not to return the positions `and` the AtomGroup of residues.

        Returns
        -------
        None or np.ndarray
            None if ``inplace=True``, numpy array if inplace=False with the positions of the center of mass of residues (if pos_type="COM") or positions of all atoms (pos_type="all")
        """

        logging.debug('Getting positions from frame',self.startF, 'to', self.endF,'with steps of',self.stepF)
        if select is None and getselection:
            logging.warning("Ignoring `getselection=True` since `select` is None.")
        ag=self.atom_group
        if select:
            ag=self.universe.select_atoms(select)
        pos=np.zeros((int((self.endF-self.startF)/self.stepF),len(ag.residues),4))
        if pos_type=='all':
            pos=np.zeros((int((self.endF-self.startF)/self.stepF),len(ag.atoms.positions),4))

        j=0
        for ts in self.universe.trajectory[self.startF:self.endF:self.stepF]:
        # for ts in self.universe.trajectory[self.startF:self.endF:self.stepF]:
            if pos_type=='COM':
                pos[j,:,0]=ts.time/1000
                pos[j,:,1:]=[r.atoms.center_of_mass() for r in ag.residues]
                # print('Getting COMs..')
            if pos_type=='all':
                pos[j,:,0]=ts.time/1000
                pos[j,:,1:]= ag.atoms.positions
            j+=1
        if not self.surf_pos is None and surf_is_zero==True:
            dict_axis={'x':1,'y':2,'z':3}
            # pos_centered=pos-np.array([0,to_center[0],to_center[1],self.surf_pos[2]])
            pos[:,:,dict_axis[self.surf_axis]]-=self.surf_pos[dict_axis[self.surf_axis]]
        if inplace:
            self.pos=np.array(pos)
            return None
        elif getselection:
            return np.array(pos),ag
        else:
            return np.array(pos)

    def getCOMs(self, inplace=True, select=None):
        """Computes positions of selection from ``self.startT`` to ``self.endT`` with ``self.stepT`` steps of time.
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
           None if ``inplace=True``, numpy array if ``inplace=False`` with the center of mass of the AtomGroup.
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

    @staticmethod
    def FilterMinFrames(pos, zlim,Nframes,control_plots=False,getNframes=False):
        """Selects a set of ``Nframes`` in which the AtomGroup is closer to the surface and bellow a zlim threshold distance to the surface.

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
        getNframes : bool, optional
            If True, returns the number of frames selected, by default False

        Raises
        ------
        ValueError
            Error raised if ``Nframes`` is bigger than the number of frames within the ``zlim`` angs threshold from the surface.
        
        Returns
        -------
        np.ndarray
            Filtered positions (Nframes,Nresidues or Natoms,4 <t,x,y,z>)
        """
        # pos=self.pos
        ##Take the mean of all selected residues
        # mean_z_top=pos[:,:,3].mean(axis=1)
        mean_z_top=pos[:,:,3].min(axis=1)
        print(mean_z_top.shape)
        # Select frames where the mean of selected residues is < zlim
        # to garanty that same frames are used for all residues.
        zMask= mean_z_top < zlim
        pos_masked=pos[zMask]
        print('There are', len(pos_masked),' frames < %i A in Z'%zlim)

        if Nframes is None:
            pos_masked=pos_masked[np.argsort(pos_masked[:,:,3].mean(axis=1),axis=0)]
        elif len(pos_masked) >= Nframes:
            print('Taking', Nframes, 'closest frames to surface...')
            pos_masked=pos_masked[np.argsort(pos_masked[:,:,3].mean(axis=1),axis=0)][:Nframes]
        else:
            raise ValueError(f"There are less than {Nframes} within the {zlim} angs threshold from the surface. Set `Nframes=None` if you want to consider all the frames < zlim. Otherwise, educe `Nframes` or  increase `zlim`.")

        if control_plots:
            ires=[0,1] ## If want to change defaut residue to make control plots, change here
            # print(f"Doing control plot with residue {self.atom_group.residues[ires]}")
            plt.plot(pos[:,ires,0],pos[:,ires,3],)
            plt.plot(pos_masked[:,ires,0],pos_masked[:,ires,3],'.',ms=2)
            plt.show()
        if getNframes:
            return pos_masked, len(pos_masked)
        else:
            return pos_masked


    def PolarAnalysis(self,select_res,Nframes, zlim=14,max_hist=None,sort=None,ax=None,plot=False,control_plots=False,Nbins=1000,resolution=5):
        """Makes a Polar Histogram of the positions of the center of mass of ``select_res`` residues considering ``Nframes`` closest to the surface within the < zlim threshold. 
        ``self.pos`` attribute is used to compute the center of mass of the AtomGroup, which will be the referential center of the histograms.
        The colors of the histogram are ordered according to sort parameter.

        Parameters
        ----------
        select_res : str
            MDAnalysis string selection of the residues to which compute the histograms.
        Nframes : int
            Number of frames closest to the surface to consider that are at ``zlim`` distance from ``self.surf_pos``. ``self.surf_pos`` `must not` be `None`. 
            If ``Nframes`` > `the total number of frames` in tha ``zlim`` adsorption threshold,the total of frames in the thresh are considered. 
        zlim : float, optional
            Distance (in angstroms) threshold limit in which the AtomGroup is considered adsorped to the surface, by default 14
        max_hist : None or float, optional
            Value to normalize the histograms. If None, highest histogram value is used to normalize the histograms, by default None
        sort : str, list, ndarray or None, optional
            How to sort the histograms in the legend and the coloring. If "max", histograms will be ploted descending from the residue with highest peak in its histogram to the flattest peak. If sort is a list or ndarray, sort is the index positions of the residues to which reorder the histograms. If None, they well be plotted by MDAnalysis default sort (ascending Resid values)., by default None
        ax : matplotlib.axes.Axes, optional
            Set plot on a predefined Axis. If None, function defines if own axis internally. By default None
        plot : bool, optional
            Show the polar plot (True) or only return the data (False), by default False
        control_plots : bool, optional
            Show control plots of the different steps of the polar analysis calculation., by default False
        Nbins : int, optional
             How many bins to use in the histograms., by default 1000
        resolution : float, optional
            One the position vectors of each residue are normalized, resolution is the value too which the normalized positions are multiplied. An increase in this value will make higher peaks in the histograms since position vector are further away. Reducing this value represents an increase of resolution of the histogram. By default 5

        Raises
        ------
        KeyError
            Error raised if ``self.surf_axis`` is not 'x', 'y' or 'z'.
        
        Returns
        -------
        list,np.ndarray
            Histogram data (Nresidues, 2 <X_bin_values,Y_bin_values>, Nbins),
            Positions in order in which they were plotted (Nframes,Nresidues or Natoms,4 <t,x,y,z>)
        """
        if self.surf_axis in ['x','y']:
            raise KeyError("`PolarAnalysis` method is currently only working with surfaces with normal axis in 'z' direction. Future versions will include computing this method for 'x', 'y' and 'z' directions.")
        elif self.surf_axis is None:
            raise KeyError("Must indicate the normal axis of the surfaces with `surf_axis` attribute. `PolarAnalysis` method is currently only working with surface with normal axis in 'z' direction. Future versions will include computing this method for 'x', 'y' and 'z' directions.")

        # colors=['C%s'%(c+7) for c in range(10)]
        prot_center=self.atom_group
        to_center=prot_center.atoms.center_of_mass()

        fromF=self.startF
        endF=self.endF
        print(f"Computing Polar Analysis from frame {fromF} (t={self.startT}ns) to {endF} (t={self.endT}ns) ")
        pos=self.getPositions(select=select_res,inplace=False,surf_is_zero=True)
        prot=self.universe.select_atoms(select_res)

        print(prot.residues)
        print(pos.shape)
        # sys.exit()
        print(pos.mean(axis=(0,1)), 'pos mean')
        # if self.surf_axis=='x':
        #     pos_centered=pos-np.array([0,0,to_center[1],to_center[2]])
        # if self.surf_axis=='y':
        #     pos_centered=pos-np.array([0,to_center[0],0,to_center[2]])
        if self.surf_axis=='z':                                        # If surf_axis if not 'z' by default, this must change slightly
            pos_centered=pos-np.array([0,to_center[0],to_center[1],0])
        else: 
            raise KeyError("Must indicate the normal axis of the surfaces with `surf_axis` attribute. `PolarAnalysis` method is currently only working with surface with normal axis in 'z' direction. Future versions will include computing this method for 'x', 'y' and 'z' directions.")
            # raise KeyError("Must indicate the normal axis of the surface with surf_axis parameter. Options are 'x','y' or 'z'")
        # if not self.surf_pos is None:
        #     dict_axis={'x':1,'y':2,'z':3}
        #     # pos_centered=pos-np.array([0,to_center[0],to_center[1],self.surf_pos[2]])
        #     pos_centered[:,:,dict_axis[self.surf_axis]]-=self.surf_pos[dict_axis[self.surf_axis]]
        print(pos_centered.mean(axis=(0,1)),'pos centered')
        
        pos_selected=self.FilterMinFrames(pos_centered,zlim,Nframes,control_plots=control_plots)
        print(pos_selected.shape)

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
        # print(pos_selected.shape, ordered_selected_pos.shape)
        # print(pos_selected[:,:3].shape, ordered_selected_pos[:,:3].shape)
        if plot==True:
            norm_max_hist=1/max_hist
            if ax is None:
                ax = plt.subplot(111, polar=True)  # Create a new figure and axes if not provided
            for i in range(ordered_selected_pos.shape[1]):
                fig_label='%i-%s'%(prot.residues.resids[sort_i[i]],prot.residues.resnames[sort_i[i]])

                bars = ax.bar(hist_arr[sort_i[i],0], hist_arr[sort_i[i],1]*max_height*norm_max_hist,
                        width=width,label='%s'%fig_label,
                        bottom=bottom,#color=colors[i],
                        alpha=0.8
                        )
            ax.set_yticklabels([])
            ax.legend(loc='center', bbox_to_anchor=(0.5, 0.5),framealpha=1)
            ax.set_ylim((0,2))
            # plt.title('Polar Histogram RBD-%s %s variant'%('PBL', plabels[ivar]))
            plt.tight_layout()
            # plt.show()
        #     ax.set_theta_zero_location('N')
        return hist_arr,ordered_selected_pos

    ######## Radii of Gyration 2D Analysis ###################

    def computeRG2D(self, masses, total_mass):
        r"""Computes parallel, perpendicular and 3D radius of gyration in 1 frame.

        .. math:: R_{\textrm{g}\parallel}= \sqrt{ \frac{1}{m_T}\sum_{i} m_{i}\left[ (x_i-x_{\textrm{CM}})^2+(y_i-y_{\text{CM}})^2\right]}

        .. math:: R_{\textrm{g}\perp} = \sqrt{\frac{1}{m_T}\sum_{i} m_{i} (z_i-z_{\text{CM}})^2,}

        where :math:`{\bf R}_{\textrm{CM}}=(x_{\textrm{CM}}`, :math:`y_{\textrm{CM}}`, :math:`z_{\textrm{CM}})` is the position of the center of mass, :math:`m_{i}` the mass of each residue and :math:`m_T` the total mass of the residues.


        Parameters
        ----------
        masses : np.ndarray (Natoms)
            Mass vallues of each atom.
        total_mass : float, optional
            Sum of these masses

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
            3D, perpendicular and parallel radius of gyration values at each frame (``self.endF``-``self.startF`` frames,3, ``Natoms``)
        """

        # Define the classic color cycle
        colors = ['b', 'm', 'c', 'r', 'g', 'y', 'k']
        rg_arr=np.zeros((len(self.universe.trajectory[self.startF:self.endF:self.stepF]),4))
        i=0
        masses=self.atom_group.atoms.masses
        total_mass=np.sum(masses)
        # This section could be optimized to use previously computed positions and coms to avoid recomputing them 
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
            plt.ylabel(r'Radius of gyration ($\mathrm{\AA}$)')
            # plt.show()
        return rg_arr

    def RgPerpvsRgsPar(self,rgs,color, marker='s',plot='both',ax=None,legend=True,show=False,mfc='k',markersize=10,system_label=None):
        r"""Generates :math:`R_{g\perp}` vs. :math:`R_{g\parallel}` plots. Also, returns the :math:`\langle R_{g\perp}^2 \rangle /\langle R_{g\parallel}^2 \rangle` ratio

        Parameters
        ----------
        rgs : np.ndarray
            Radii of gyration data. It must have Rg data as [Rg 3D, Rg parallel, Rg perpendicular], analogue to getRgs2D functions
        color : str
            Color used to plot points. Color names use the same of those of Matplotlib package.
        marker : str, optional
            Marker used to plot. Marker names use the same of those of Matplotlib package. , by default 's'
        plot : str, optional 
            Whether to plot only the raw data (``plot='data'``), only the mean value of the data (``plot='mean'``), or both (``plot='both.'``)
        ax : matplotlib.axes.Axes, optional
            Set plot on a predefined Axis. If None, function defines if own axis internally. By default None
        legend : bool, optional
            Whether or not to make the default legend, by default True
        show : bool, optional
            If True, matplotlib.pyplot.show() is executed. This is left optional in case further characterization of plot is desired. Also, this enables showing multple data in the same figure. By default False.
        mfc : str, optional
            Set marker face color. Color options are intrinsic colors used in Matplotlib library. 
        markersize : int, optional
            Size of colors. Uses the same values as intrinsic values of Matplotlib library. 
        system_label : str, optional
            Label to assign to the system in the legend. If None, the system name is used. By default None
        
        Raises
        ------
        KeyError
            Error raised if ``plot`` is different than 'data', 'mean' or 'both'.
            Error raised if ``legend=True``,``system_label`` is None and the system name is not defined.
        Returns
        -------
        float
            :math:`\langle R_{g\perp}^2 \rangle /\langle R_{g\parallel}^2 \rangle`
        """
        data=rgs[:,2:]
        print(data.shape)
        rg_ratio=(data[:,0]**2).mean()/(data[:,1]**2).mean()
        if system_label is None:
            system_label=self.system_name
            
        if ax is None:
            fig, ax = plt.subplots()  # Create a new figure and axes if not provided
        if plot=='both':
            label='%s (%.3f)'%(system_label,rg_ratio)
            ax.plot(data[:,1],data[:,0],'o',markersize=1,c=color,)
            ax.plot(data[:,1].mean(),data[:,0].mean(),marker,markersize=markersize, label=label,color='k',mfc=mfc)
        elif plot=='data':
            label='%s (%.3f)'%(system_label,rg_ratio)
            ax.plot(data[:,1],data[:,0],'o',markersize=1,c=color)
            ax.plot([], [], 'o', c=color,label=label)  # Add a single label for all data points
        elif plot=='mean':
            label='%s (%.3f)'%(system_label,rg_ratio)
            # ax.plot(data[:,1].mean(),data[:,0].mean(),label=label,**kwargs)
            ax.plot(data[:,1].mean(),data[:,0].mean(),marker,markersize=markersize, label=label,color='k',mfc=mfc)
        else:
            raise KeyError("`plot` must be 'data', 'mean' or 'both' ")
        ax.set_xlabel(r'$Rg_\parallel$ ($\mathrm{\AA}$)')
        ax.set_ylabel(r'$Rg_\perp$ ($\mathrm{\AA}$)')
        if legend:
            if not system_label is None:
                ax.legend(title=r'Syst ($\langle Rg_\perp^2\rangle /\langle Rg_\parallel^2 \rangle$)', fontsize=10)
            else:
                raise KeyError("Must indicate the system name with `system_name` attribute to make a legend, or assign a label to the system with `system_label` parameter.")
        if show:
            plt.show()
        return rg_ratio

    ############# Compute Contour Area #################
    @staticmethod
    def makeCmapColor(def_color):
        RGB_color=colors.to_rgba(def_color)
        cdict1 = {
            'red': (
                (0.0,RGB_color[0]*0.6, RGB_color[0]*0.6),
                # (0.5, 0.0, 0.1),
                (1.0, RGB_color[0], RGB_color[0]),
            ),
            'green': (
                (0.0, RGB_color[1]*0.6, RGB_color[1]*0.6),
                (1.0,RGB_color[1], RGB_color[1]),
            ),
            'blue': (
                (0.0,  RGB_color[2]*0.6, RGB_color[2]*0.6),
                # (0.5, 0.1, 0.0),
                (1.0, RGB_color[2], RGB_color[2]),
            ),
            'alpha': (
                (0.0, 0.5, 0.5),
                # (0.5, 0.1, 0.0),
                (1.0, 1.0, 1.0),
            )
        }
        cmap_color = LinearSegmentedColormap('Browns', cdict1)
        return cmap_color
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

    def getKDEAnalysis(self,zlim,Nframes,axis=['x','y'],inplace=True,ax=None,show=False,control_plots=False):
        """Computes KDE Contours using seaborn.kde_plot() function and extracts the paths of each contour level. The output of the seaborn.kde_plot() is stored in self.kdeanalysis.kde, and the paths of each contour level is stored in self.kdeanalysis.paths if inplace=True.

        Parameters
        ----------
        zlim : float
            zlim of BioPolymer2D.FilterMinFrames(). Only use frames under a zlim threshold, to avoid using frames with desorbed molecule.
        Nframes : int
            Nframes of BioPolymer2D.FilterMinFrames(). To ensure to have a controled number of frames under zlim threshold.
        axis : list, optional
            Axes for which to make the KDE plot. Default is x and y axis.
        inplace : bool, optional
            If True, stores the paths of al contour levels in self.kdeanalysis.paths. Otherwise, it only returns it. By default True
        ax : matplotlib.axes.Axes, optional
            Set plot on a predefined Axis. If None, function defines if own axis internally. By default None
        show : bool, optional
            If True, matplotlib.pyplot.show() is executed. This is left optional in case further characterization of plot is desired. Also, this enables showing multple data in the same figure. By default False.
        control_plots : bool, optional
            Make control plots, by default False

        Returns
        -------
        list
            List of all paths in all the contour levels.
        """
        pos=self.pos
            
        if not self.surf_pos is None:
            dict_axis={'x':1,'y':2,'z':3}
            # pos_centered=pos-np.array([0,to_center[0],to_center[1],self.surf_pos[2]])
            pos[:,:,dict_axis[self.surf_axis]]-=self.surf_pos[dict_axis[self.surf_axis]]
        pos_selected=self.FilterMinFrames(pos,zlim,Nframes,control_plots=control_plots)
        ## Concatenate positions of all residues
        print(pos_selected.shape)
        pos_selected_reshape=np.reshape(pos_selected,(pos_selected.shape[0]*pos_selected.shape[1],pos_selected.shape[2]))
        print(pos_selected_reshape.shape)
        if ax is None:
            fig, ax = plt.subplots()  # Create a new figure and axes if not provided
        df0=pd.DataFrame(pos_selected_reshape, columns=['t','x','y','z'])
        kde_plot = sns.kdeplot(df0, x=axis[0], y=axis[1], fill=True, cmap="Purples", color='black',ax=ax)
        kde_plot=sns.kdeplot(df0, x=axis[0],y=axis[1], color = 'black',alpha=0.5,ax=ax)
        self.kdeanalysis.kde=kde_plot

        paths_arr=[]
        Nlvls=len(kde_plot.collections[-1].get_paths())
        print(f"There are {Nlvls} levels in the KDE.")
        for lvl in range(Nlvls):
            paths=self.ListPathsInLevel(kde_plot,lvl,plot_paths=control_plots)
            if not paths:
                continue
            # print(np.shape(paths[0]))
            # print(np.shape(paths[1]))
            paths_arr.append(paths)
            # plt.pause(2)
            if control_plots:
                plt.pause(2)
        plt.gca().set_aspect('equal', 'box')
        if not show:
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

    @staticmethod
    def getAreas(paths,contour_lvl,getTotal=False):
        """Computes the area of each path a given contour level. Negative values of area are holes in the contour level. If ``getTotal=True``, computes the area of the whole contour level.

        Parameters
        ----------
        paths : list
            Paths of all contour levels.         
        contour_lvl : int
            Contour level to compute the area.
        getTotal : bool, optional
            If False, gives the area of each path of the contour level. If ``getTotal=True``, sums up the contribution of each path in the contour level returning the area of the whole contour level. By default False

        Returns
        -------
        list or float
            A list with the area of each path in the contour level, or the total area of the contour level (if getTotal=True)
        """

        Areas=[]
        paths_in_lvl=paths[contour_lvl]
        for i_paths in range(len(paths_in_lvl)):
            x_values,y_values=paths_in_lvl[i_paths].T
            area_outline = simpson(y_values[::-1], x=x_values[::-1])
            print("Area of the outline %s:"%i_paths,area_outline)
            Areas.append(area_outline)
        if getTotal:
            return np.sum(Areas)
        else:
            return Areas
        
        
    def KDEAnalysisSelection(self,select_res,Nframes,zlim=15,ax=None,show=False,legend=False,plot_COM=True,getNframes=False):
        r"""KDE Contours for a set of selected residues. This computes the paths of all the contour levels of each residue.

        Parameters
        ----------
        select_res : str
            MDAnalaysis-like selection of residues to compute their KDE Contour paths
        Nframes : int, optional
            Number of frames to use within :math:`mean(z_{select\_res})` < ``zlim``. This fills value of :code:`self.FilterMinFrames`.
        zlim : float, optional
            Height limit to consider as adsorted frames. This fills value of :code:`self.FilterMinFrames`, by default 15
        ax : matplotlib.axes.Axes, optional
            Set plot on a predefined Axis. If None, function defines if own axis internally. By default None
        show : bool, optional
            Where or not to show plot. Convinient to use False if you want modify the default plot, by default False
        legend : bool, optional
            Whether or not to make the default legend, by default False
        plot_COM : bool, optional
            Whether or not to plot the center of mass of ``select_res``. Particularly relevant if comparing results with ``PolarAnalysis`` method, 
            since this will be the center of the polar histograms if ``select_res`` are the same in both analysis. By default True
        getNframes : bool, optional
            If True, returns the number of frames selected, by default False

        Returns
        -------
        (list(list)), AtomGroup
            Returns a list with the contour levels paths of each selected residue, and selected AtomGroup. The list of contour levels paths and the reisudes in the AtomGroup are ordered consistently.
        """

        def_colors=["C%i"%(i) for i in range(10)]

        COM=self.atom_group.center_of_mass()

        pos_res_contour, res=self.getPositions(select=select_res,inplace=False,getselection=True,surf_is_zero=True,)
        ### Maybe this section can be taken of if surf_is_zero=True?
        # if not self.surf_pos is None:
        #     dict_axis={'x':1,'y':2,'z':3}
        #     # pos_centered=pos-np.array([0,to_center[0],to_center[1],self.surf_pos[2]])
        #     pos_res_contour[:,:,dict_axis[self.surf_axis]]-=self.surf_pos[dict_axis[self.surf_axis]]
        print(pos_res_contour.mean(axis=(0,1)))
        # fig,ax=plt.subplots()   
        # ListPaths(vertices,codes,plot_paths=True)
        # print(RBM.residues)
        # print(pos.shape)
        if getNframes:
            all_pos_selected,n_used_frames=self.FilterMinFrames(pos_res_contour, zlim,Nframes,control_plots=False,getNframes=True)
        else:
            all_pos_selected=self.FilterMinFrames(pos_res_contour, zlim,Nframes,control_plots=False)
        print(all_pos_selected.shape)
        # all_pos_selected_reshaped=np.reshape(all_pos_selected,(all_pos_selected.shape[0]*all_pos_selected.shape[1],4))
        # print(all_pos_selected_reshaped.shape)
        Nres=len(all_pos_selected[0])
        if plot_COM:
            plt.plot(COM[0],COM[1],c='k',marker='o')
        resnames=[]
        handles=[]
        paths_arr_arr=[]
        if ax is None:
            fig, ax = plt.subplots()  # Create a new figure and axes if not provided
        for ires in range(Nres):
            res_pos=all_pos_selected[:,ires]
            df0=pd.DataFrame(res_pos, columns=['t','x','y','z'])
            # kde_plot=sns.kdeplot(df0, x='x',y='y', color = 'black',alpha=1,fill=True)
            cmap_color=self.makeCmapColor(def_colors[ires])
            kde_plot=sns.kdeplot(df0, x='x',y='y', fill=True,cmap=cmap_color,ax=ax)
            # Create a legend handle manually
            handles.append(Line2D([0], [0], color=def_colors[ires], lw=4))
            resnames.append(f"{res.residues[ires].resid}-{res.residues[ires].resname}")
            paths_arr=[]
            Nlvls=len(kde_plot.collections[-1].get_paths())
            # print(f"There are {Nlvls} levels in the KDE.")
            for lvl in range(Nlvls):
                paths=self.ListPathsInLevel(kde_plot,lvl,plot_paths=False)
                if not paths:
                    continue
                # print(np.shape(paths[0]))
                # print(np.shape(paths[1]))
                paths_arr.append(paths)
            paths_arr_arr.append(paths_arr)

        ax.set_title(f'Contour Positions {self.system_name}')
        ax.set_xlabel(r'X ($\AA$)')
        ax.set_ylabel(r'Y ($\AA$)')
        plt.gca().set_aspect('equal', 'box')

        # Add the custom legend
        if legend:
            plt.legend(handles=handles, labels=resnames, loc='upper right')
        if show:
            plt.show()
        if getNframes:
            return paths_arr_arr, res, n_used_frames
        else:
            return paths_arr_arr, res

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
        
        Raises
        ------
        SelectionError
            Error raised if any of the two selections lack of hydrogens or acceptors.

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
        hbonds.run(verbose=True,start=self.startF,stop=self.endF,step=self.stepF)

        if inplace:
            self.hbonds=hbonds.results

        if trj_plot:
            plt.plot(hbonds.times/1000, hbonds.count_by_time(), lw=2, label='%s'%self.system_name)
            plt.xlabel("Time (ns)")
            plt.ylabel(r"$N_{HB}$")
            plt.show()
        return hbonds.results

    def HbondsPerResidues(self,sorted=True, unique_res=False):
        """Computes the number of H-bonds of each residue during the simulation. ``self.getHbonds(inplace=True)`` must be computed prior to the use this function.

        Parameters
        ----------
        sorted : bool, optional
            If True, returns data sorted by number of H-bonds. By default True
        unique_res : bool, optional
            Whether or not to consider only one hydrogen bond per residue. Important to consider percentage in which a residue participates in a H-bond. By default False
        Returns
        -------
        pandas.DataFrame
            DataFrame showing all the residues with H-bonds and H-bond count
        """
        result=np.array(self.hbonds.hbonds[:,[2,3]], dtype=int) #Indexes of Hydrogen and acceptors starting from 0.
        # print(result[0,0],self.universe.atoms[result[0,0]-1], 'result check') ## [result[0,0]-1 is an O, not a hydrogen, not correct to substract 1
        # resids=np.array([self.universe.atoms[result[:,0]].resids,self.universe.atoms[result[:,1]].resids])
        # df=pd.DataFrame(data=resids.T, columns=['Hydrogens', 'Acceptors'])
        resids=np.array([self.hbonds.hbonds[:,0],self.universe.atoms[result[:,0]].resids,self.universe.atoms[result[:,1]].resids], dtype=int)
        df=pd.DataFrame(data=resids.T, columns=['Frame','Hydrogens', 'Acceptors'])
        # print(df[df['Hydrogens']==199].iloc[:10])
        print(resids.shape, result.shape, 'resid', 'result')
        print(df)
        if unique_res:
            df.drop_duplicates(subset=['Frame','Hydrogens'],inplace=True,ignore_index=True)
            df.drop_duplicates(subset=['Frame','Acceptors'],inplace=True,ignore_index=True)
        print(resids.shape, result.shape, 'resid', 'result')
        print(df)
        Acc_resids = df[np.isin(df['Acceptors'],self.atom_group.residues.resids)].pivot_table(columns=['Acceptors'], aggfunc='size')
        H_resids = df[np.isin(df['Hydrogens'],self.atom_group.residues.resids)].pivot_table(columns=['Hydrogens'], aggfunc='size')
        # print(H_resids,Acc_resids)
        final_count=Acc_resids.add(H_resids, fill_value=0) # dividing by 2 is not necesarry since a a atom will be a donot and acceptor at the same time. 
        print(final_count.index-1, 'final_count.index')
        df_final=pd.DataFrame(np.array([self.universe.residues[final_count.index-1].resids,self.universe.residues[final_count.index-1].resnames,final_count]).T,
                            columns=['ResIDs','ResNames','Count'])

        if sorted:
            return df_final.sort_values('Count', ascending=False)
        else:
            return df_final
    def plotHbondsPerResidues(self,
                            paths_for_contour,
                            max_circle_size=160,
                            top=-1,contour_lvls_to_plot=None,
                            contour_colors=None,
                            contour_ls=None,
                            contour_alphas=None,
                            filter=None,
                            print_table=True,
                            show=False,
                            ax=None):
        """Makes a figure showing the center of mass of the residues with H-bonds. Figure shows a contour plot as a reference of position of the whole molecule. Legend of the Figure shows the percentage of time in which there were Hbonds during the simulation of the plotted residues.

        Parameters
        ----------
        paths_for_contour : list
            List of paths of all the contour levels.
        max_circle_size : int, optional
            Maximum size of the circle representing the residue with most H-bonds. By default 160
        top : int, optional
            Residues are plotted ranked by residues with most contact to least. This parameters indicates how many residues to plot of these ranked residues, e.g. top=5 wil plot the 5 residues with most Hbonds during the simulations. By default -1, plots all the residues with H-bonds.
        contour_lvls_to_plot : list, optional
            Contour Levels to show in plot, by default None
        contour_colors : list, optional
            Colors to use to show the reference contour levels.This list must be the same size than `contour_lvls_to_plot` parameter. Default None, which all reference contour plots are shown in black.
        contour_ls : list, optional
            Line style of the contour levels. Default is a solid line.
        contour_alphas : list, optional
            Transparency of the contour levels. Default is 0.3.
        filter : str, optional
            Filter out a residue from the plot. By default None
        print_table : bool, optional
            Whether or not to print the pandas.DataFrame table with the data shown in figure, by default True
        ax : matplotlib.axes.Axes, optional
            Set plot on a predefined Axis. If None, function defines if own axis internally. By default None
        show : bool, optional
            Whether or not to show the figure. Convinient to use False if want to do further tunning the plot. By default False

        Returns
        -------
        pandas.DataFrame
            Sorted DataFrame from the residues with most contact to the residue with least contacts.
        """

        df=self.HbondsPerResidues(sorted=False, unique_res=True)
        df_resname=df['ResNames']
        if filter is not None:
            if isinstance(filter, list):
                df = df[~df_resname.isin(filter)]
            else:
                df = df[df_resname != filter]
        # max_val=df['Count'].max()
        str_resids=' '.join(np.array(df['ResIDs'],dtype=str))
        print(str_resids)

        pos=self.getPositions(select=f'resid {str_resids}', inplace=False,surf_is_zero=False)
        df[['X','Y', 'Z']]=pos[:,:,1:].mean(axis=0)
        sorted_df=df.sort_values('Count', ascending=False)
        if print_table:
            print(sorted_df.iloc[:top])

        if not contour_lvls_to_plot:
            contour_lvls_to_plot=range(len(paths_for_contour))
        if not contour_colors:
            contour_colors=['k' for _ in paths_for_contour]
        if not contour_ls:
            contour_ls=['-' for _ in  paths_for_contour]
        if not contour_alphas:
            contour_alphas=[0.3 for _ in  paths_for_contour]
        if ax is None:
            fig,ax = plt.subplots()  # Create a new figure and axes if not provided

        for lvl,c,i in zip(contour_lvls_to_plot,contour_colors, range(len(contour_lvls_to_plot))):
            # self.plotPathsInLevel(paths_for_contour,lvl,color=c,ax=ax)
            self.plotPathsInLevel(paths_for_contour,lvl,ax=ax,ls=contour_ls[i],alpha=contour_alphas[i])

        colors = ['C%s' % i for i in range(10)]  # Define color palette
        num_colors = len(colors)
        for i in range(pos[:,:top].shape[1]):
            # Use modular indexing to cycle through colors
            color = colors[i % num_colors]
            norm_val=sorted_df['Count'].iloc[i]/len(self.universe.trajectory) #max_val
            # norm_val_plot=sorted_df['Count'].iloc[i]/max_val
            pos=sorted_df[['X','Y','Z']].iloc[i]
            ax.plot(pos['X'],pos['Y'], 'o',color=color,
                    label='%s-%s (%.2f)'%(sorted_df['ResIDs'].iloc[i],
                                        sorted_df['ResNames'].iloc[i],
                                        norm_val*100),)
            # ax.scatter(pos['X'],pos['Y'], s=(max_circle_size*norm_val_plot)**2, alpha=.5, color=color)
            ax.scatter(pos['X'],pos['Y'], s=(max_circle_size*norm_val)**2, alpha=.5, color=color)
        ax.set_xlabel(r'X-axis($\AA$)',)#fontsize=20)
        ax.set_ylabel(r'Y-axis($\AA$)',)#fontsize=20)
        plt.gca().set_aspect('equal')
        plt.tight_layout()
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),#prop={'size':22},
                    title="ResID-ResName(Hbond %)",)#title_fontsize=20)
        if show:
            plt.show()
        return sorted_df

    @staticmethod
    def plotPathsInLevel(paths, contour_lvl,color='k',alpha=0.3,ls='-',lw=None,show=False,ax=None,):
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
        ax : matplotlib.axes.Axes, optional
            Set plot on a predefined Axis. If None, function defines if own axis internally. By default None
        show : bool, optional
            Where to show or not the plot yet with matplolib.pyplot.plot, by default False
        """
        if ax is None:
            fig,ax = plt.subplots()  # Create a new figure and axes if not provided
        paths_in_lvl=paths[contour_lvl]
        for p in range(len(paths_in_lvl)):
            x_val,y_val=paths_in_lvl[p].T
            ax.plot(x_val,y_val,color=color, alpha=alpha,ls=ls,linewidth=lw)
            # ax.plot(x_val,y_val,**kwargs)
        if show:
            plt.show()

