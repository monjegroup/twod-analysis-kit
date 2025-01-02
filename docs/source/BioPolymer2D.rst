BioPolymer2D
----------------------------------------


This quick user guide includes instructions on how compute the following analysis:

 - Polar histogram analysis
 - Kernel Density Estimation contours
 - Parallel and perpendicular radii of gyrations
 - Hydrogen bonds per residues

Before any analysis, we must import ``MDAnalysis`` (since our class is initialized with a 
MDAnalysis AtomGroup or Universe), and our class.

.. code-block:: python

    import MDAnalysis as mda
    import matplotlib.pyplot as plt
    from twodanalysis import BioPolymer2D



And we initialize our object using a MDAnalysis Universe, or a AtomGroup. 

.. code-block:: python

    tpr = "md.tpr" # Replace this with you own topology file
    xtc = "md.xtc" # Replace this with your trajectory file

    universe = mda.Universe(tpr,xtc) # Define a universe with the trajectories
    sel=universe.select_atoms("protein") # Use any convinient selection.

    biopol = Biopolymer2D(sel)   # Initialize object by loading selection.

If `BioPolymer2D` is initialized with a Universe, the whole universe is considered to be the selection. Note 
that the biopolymer is not necesarilly a protein but can also be a nucleic acid. With ``INFO`` method we can 
retreive general information on our universe and selection, and we can set a system name with ``system_name`` 
attribute. Especially convinient if working with more than object at the same time, since the names will appear in the legends. 

.. code-block:: python

    ag_analysis.system_name='Omicron PBL1'
    ag_analysis.INFO()

In general, we would also like to compute the positions of the residues in our object. This will store position values of each frame on 
the ``pos`` attribute. Compute:

.. code-block:: python

    ag_analysis.getPositions()

If you want to consider only a time section of your whole trajectory set the attributes ``startT``, ``endT``, and ``stepT`` before 
computing ``getPositions``,e.g.

.. code-block:: python
    
    ag_analysis.startT=100
    ag_analysis.endT=200
    ag_analysis.stepT=0.4
    ag_analysis.INFO()
    ag_analysis.getPositions()


``INFO`` to confirm that  ``startT``, ``endT``, and ``stepT`` have been overwriten.

Polar histogram analysis
^^^^^^^^^^^^^^^^^^^^^^^^

Since we are interested in only sampling the adsorption, ``PolarAnalysis`` method filters the frames in which the object is not 
adsorbed using a ``zlim`` and ``Nframes`` parameters. Frames in which the center of mass of the object  is lower than ``zlim`` to the surface,
are considered adsorbed.

Now, in general, the number of adsorbed frames will vary for different trajectories, and we would like to compare between trajectories. To this matter,
``Nframes`` paramater will set the number of frames we want to take from the total adsorbed frames, take the ``Nframes`` most adsorbed.

Then, we compute the  ``PolarAnalysis``, setting these parameters,

.. code-block:: python

    select_res='resid 198 200 12 8 40 45 111 115 173'
    zlim=15
    Nframes=900
    hist_arr,pos_hist=ag_analysis.PolarAnalysis(select_res,Nframes, 
                                                zlim=zlim,control_plots=False,plot=True)
    plt.show()

If we only want to compute the histogram, set ``plot=False``. ``control_plots`` is to visualize the diferent steps of the PolarAnalysis calculations.
Titles and further figure costumization can be added to the plot using standard ``matplotlib.pyplot`` methods before ``plt.show()``.

.. note::
    Typically, the surface in the trajectory will no be set in ``z=0``. We suggest overwriting the ``surf_pos`` attribute with the position of the 
    surface (<x,y,z>) before computing the ``PolarAnalysis`` method. Only the z value will be used. 
    
    .. code-block:: python

        surface_selection='resname DOL and name O1 and prop z > 16'
        surface_pos=ag_analysis.getPositions(select=surface_selection, inplace=False)
        ag_analysis.surf_pos=surface_pos

    With the ``inplace=False`` it will not overwrite the ``pos`` attribute of the object, but only return it.
    
Kernel Density Estimation (KDE) contours
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In general, we would like to have a reference of the position of the whole biopolymer to have insight ont the flexible regions. Therefore, 
we first compute the KDE of whole molecule, and then compute the KDE of selected residues:

.. code-block:: python

    paths=biopol.getKDEAnalysis(zlim,Nframes,)
    biopol.plotPathsInLevel(paths,0,show=False)
    all_residues_paths,residues_in_contour=biopol.KDEAnalysisSelection('resid 198 200 12 8 40 45 111 115 173',Nframes,zlim,show=False,legend=True)
    plt.show()

.. note:: Setting the same ``zlim`` and ``Nframes`` paramater values for ``PolarAnalysis`` , ``getKDEAnalysis`` and ``KDEAnalysisSelection`` is suggested.

We now can compute the Areas of the paths computed by ``KDEAnalysisSelection`` with the ``getAreas`` attribute as follows:


.. code-block:: python

    data=[]
    for p in range(len(all_residues_paths)):
        areas=BioPolymer2D.getAreas(all_residues_paths[p],0,getTotal=True)
        data.append([residues_in_contour.residues[p].resid,residues_in_contour.residues[p].resname,areas])
    df=pd.DataFrame(data=data, columns=["ResIDs", "Resnames", "Area (angs^2)"])
    df

``df`` will show a table with the areas of the outer contour levels (level 0 in ``getAreas`` , is outer).


Parallel and perpendicular radii of gyrations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parallel and perpendicular radii of gyration gives structural information during the adsorption,

* :math:`R_{g\parallel}`: Gives information on how the biopolymer is expanded by the sides (parallel to the surface). 

* :math:`R_{g\perp}` : Gives information on how the biopolymer is streched or flattened.

.. figure:: TOC12_Final.png
   :alt: Example of radii of gyration correlation
   :width: 100%
   :align: center

   **Figure 1:** Example of radii of gyration correlation figures that can be made with method on the left and a schematic 
   representacion of the parallel and perpendicular radii of gyrations on the right. Figure taken from the TOC figure of `Bosch et\.al`_ (2024).


To notice significant results, we need to select a region that is in contact with the surface as our object, e.g.

.. code-block:: python

    sel_in_Contact=u.select_atoms('resid 4-15 or resid 34-45 or resid 104-117 or resid 170-176') # Select region in contact with surface
    Contact_region = BioPolymer2D(sel_in_Contact) # Initialize object
    Contact_region.system_name='Contact Omicron PBL1'# Set system name
    Contact_region.getPositions() # Compute positions
    ratio=Contact_region.RgPerpvsRgsPar(rgs, 'tab:green',show=False) # Make RgPerp vs Rg parallel plot
    
The output will be similar to Figure 1 (left), with only one system instead of six. The ``ratio`` will give
the :math:`\langle R_{g\perp}^2 \rangle /\langle R_{g\parallel}^2 \rangle` ratio, which is relevant on charactertizing the adsorption of polymers (CITE 2 PAPERS).


Hydrogen bonds per residues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. _Bosch et\.al: https://pubs.acs.org/doi/10.1021/acs.jcim.4c00460
