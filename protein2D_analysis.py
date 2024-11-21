import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import matplotlib as mpl
import sys
class protein2D_analysis:
    def __init__(self, obj):
        """
        Initializes the class with either an MDAnalysis Universe or AtomGroup.
        
        Parameters:
        obj (Universe or AtomGroup): The object to initialize with.
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

    def __repr__(self):
        return f"<{self.__class__.__name__} with {len(self.atom_group)} atoms>"
    
    def INFO(self):
        """
        Print information about the Universe and/or AtomGroup.
        """
        # Universe-level information
        print("=== UNIVERSE INFO ===")
        print("  N atoms:", len(self.universe.atoms))
        print("  N residues:", len(self.universe.residues))
        print("  N segments:", len(self.universe.segments))
        
        # AtomGroup-specific information (only if a subset is selected)
        if len(self.atom_group) < len(self.universe.atoms):
            print("=== SELECTION INFO ===")
            print("  N selected atoms:", len(self.atom_group))
            print("  N selected residues:", len(self.atom_group.residues))
            print("  N selected segments:", len(self.atom_group.segments))

    def getPositions(self,pos_type='COM', inplace=True, select=None):

        print('Getting positions from frame',self.startF, 'to', self.endF,'with steps of',self.stepF)

        prot=self.atom_group
        if select:
            prot=self.universe.select_atoms(select)
        pos=np.zeros((int((self.endF-self.startF)/self.stepF),len(prot.residues),4)) 
        if pos_type=='all':
            pos=np.zeros((int((self.endF-self.startF)/self.stepF),len(prot.atoms.positions),4)) 
   
        j=0
        for ts in self.universe.trajectory[self.startF:self.endF:self.stepF]:
            if pos_type=='CA':
                pos[j,:,0]=ts.time/1000
                pos[j,:,1:]=prot.atoms.positions
                # print('Getting C-alphas..')
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
        

    def FilterMinFrames(self, zlim,Nframes,control_plots=False):
        '''
        pos: 
        array of shape (TotalFrames,Nresidues,4 <t,x,y,z>)
        
        returns:
        array of shape (Nframes,Nresidues,4 <t,x,y,z>) of frames 
        
        '''
        pos=self.pos
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
            print(f"Doing control plot with residue {self.atom_group.residues[ires]}")
            plt.plot(pos[:,ires,0],pos[:,ires,3],)
            plt.plot(pos_masked[:,ires,0],pos_masked[:,ires,3],'.',ms=2)
            plt.show() 
        return pos_masked
    
    def PolarAnalysis(self,select_res,Nframes,max_hist=None,sort='max',plot=False,control_plots=False, zlim=15,Nbins=1000,resolution=5):
        ## ------------------####
        # Makes the histogram based on the angles respect to the center of mass of the protein
        ## -------------------####
        # RBM_resids=self.atom_group.residues.resids
        # RBM_str=' '.join(RBM_resids)


        colors=['C%s'%(c+7) for c in range(10)]
        prot_center=self.atom_group
        to_center=prot_center.atoms.center_of_mass()

        fromF=self.startF
        endF=self.endF
        print(f"Computing Polar Analysis from frame {fromF} (t={self.startT}ns) to {endF} (t={self.endT}ns) ")
        select_res_to_compute=self.universe.select_atoms(select_res)
        res_to_analyse=protein2D_analysis(select_res_to_compute)
        res_to_analyse.getPositions()
        prot=res_to_analyse.atom_group
        pos_prot=res_to_analyse.pos
        print(prot.residues)
        print(pos_prot.shape)
        # sys.exit()

        res_to_analyse.pos=pos_prot-np.array([0,to_center[0],to_center[1],0])
        pos_selected=res_to_analyse.FilterMinFrames(zlim,Nframes,control_plots=control_plots)

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
                # if res.residues.resids[i]>=193 and res.residues.resids[i]>=200:
                #     fig_label='%ith Glyc. Carb.\n%s'%(res.residues.resids[i]-192,res.residues.resnames[i])
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
    def RG2D(self, masses, total_mass=None):
        # coordinates change for each frame
        coordinates = self.atom_group.positions
        center_of_mass = self.atom_group.center_of_mass()

        # get squared distance from center
        ri_sq = (coordinates-center_of_mass)**2
        # sum the unweighted positions
        sq = np.sum(ri_sq, axis=1)
        sq_par = np.sum(ri_sq[:,[0,1]], axis=1) # sum over x and y, parallel Rg
        #print(sq_par.shape)
        sq_perp = ri_sq[:,2] # sum over zfor perpendicular Rg

        # make into array
        sq_rs = np.array([sq, sq_perp, sq_par])

        # weight positions
        rog_sq = np.sum(masses*sq_rs, axis=1)/total_mass
        # square root and return
        return np.sqrt(rog_sq)
    
    def getRgs2D(self, plot=True):
        colors=['tab:blue','tab:orange','tab:green']
        rg_arr=np.zeros((len(self.universe.trajectory),4))
        i=0
        masses=self.atom_group.atoms.masses
        total_mass=np.sum(masses)

        for ts in self.universe.trajectory:
            rg_arr[i,0]=ts.time/1000
            rgs=self.RG2D(masses, total_mass)
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
    
###  POTENTIAL UPGADE: RG2D could compute directly all frames to avoid reruning trajectory. A "getCOM" Method would be required.