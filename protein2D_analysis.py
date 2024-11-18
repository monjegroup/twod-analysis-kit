import MDAnalysis as mda
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.integrate import simpson

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

    def getPositions(self,pos_type='COM'):

        print('Getting positions from frame',self.startF, 'to', self.endF,'with steps of',self.stepF)

        prot=self.atom_group
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
        return np.array(pos)
        