Membrane Simulations
--------

To study biophysical properties of membranes in 2D we developed 3 classes: :code:`Cumulative2D`, :code:`Voronoi2D`, and :code:`PackingDefects`.
The first two classes project membrane properties and features to the membrane surface plane. Current structural properites the code computes and projects include:

 - Membrane thickness
 - Deuterium order parameters (Only for Cumulative2D)
 - Area-per-lipid (Only for Voronoi2D)
 - Splay angle

The third class :code:`PackingDefects` identifies regions where the hydrophobic membrane core is exposed. The implementation in Python is built and optimized from a method implemented in PackMem. 

Below are concise explanations and examples of the Cumulative2D, Voronoi2D, and PackingDefects workflows. For a detailed tutorial notebook visit https://github.com/pyF4all/2DanalysisTutorials/tree/main

Cumulative2D
^^^^^^^^^^^^^^^^^^^^^^^^^^

`Cumulative2D` projects membranes properties to the two dimensional plane of the membrane surface, perpendicular to the z axis. 
The image below illustrates the protocol this class uses for the projection. The process begins by dividing the space into a :math:`m\times m` grid. The xy positions of lipid phosphorus atoms and the respective analysis are collected over a user-set number of frames. These values are averaged within each grid square and stored in a :math:`m\times m` matrix. Alongside this matrix, the grid edges are recorded in the format :math:`[x_\text{min},x_\text{max},y_\text{min},y_\text{max}]`.


.. image:: cumulative.png

To import :code:`Cumulative2D`, type:

.. code-block:: python

    import MDAnalysis as mda
    from twodanalysis import Cumulative2D



To use the class, call it using :code:`mda.AtomGroup` or a :code:`mda.Universe` as follows:

.. code-block:: python

    tpr = "membrane.tpr" # Replace with your tpr or gro file
    xtc = "membrane.xtc" # Replace with your xtc file

    universe = mda.Universe(tpr,xtc) # Define a universe

    membrane = Cumulative2D(universe,   # load the universe
                    verbose = False, # Does not print intial information
                    )


.. note::
    If your trajectory contains water and/or ions, pass the list of lipids in the membrane by specifying :code:`lipid_list`.


Membrane Thickness
++++++++++++++++++

This code requires the user to set the number of bins, the edges, and the time interval. Additional options
are listed in the documentation.

.. code-block:: python

    mat_thi, edges = membrane.thickness(50,           # nbins
                                        start = 61,   # Initial frame
                                        final = 110,  # Final Frame
                                        step = 1,     # Frames to skip
                                        )

The output is a matrix :math:`nbins\times nbins` and the edges in the form
 :math:`[x_\text{min},x_\text{max},y_\text{min},y_\text{max}]`.

To visualize with :code:`plt.imshow`:

 .. code-block:: python

    import matplotlib.pyplot as plt

    plt.imshow(mat_thi, extent=edges, cmap="Spectral")
    plt.xlabel("x $\AA$")
    plt.ylabel("y $\AA$")
    plt.title("Membrane thichness from frames 61-110")
    cbar = plt.colorbar()
    cbar.set_label('Thickness $\AA$')

 .. image:: thickness.png


Membrane order parameters
+++++++++++++++++++++++++

To compute the order parameters the user must select the leaflet for which to run the analysis (top, bottom, or both) as shown below.

.. code-block:: python

    scd_top, edges = membrane.all_lip_order("top",
                                                50,
                                                start = 61,
                                                final=110,
                                                step = 1)
    scd_bot, edges = membrane.all_lip_order("bot",
                                                50,
                                                start = 61,
                                                final=110,
                                                step = 1)


To plot the results:


 .. code-block:: python

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    # Plot
    fig, ax = plt.subplots(1,2, sharex = True, sharey = True)
    first = ax[0].imshow(scd_top, extent=edges, cmap="Spectral")
    ax[0].set_xlabel("x $\AA$")
    ax[0].set_ylabel("y $\AA$")
    ax[0].set_title("Top layer")
    divider1 = make_axes_locatable(ax[0])
    cax1 = divider1.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(first, cax = cax1)
    # Point to a low ordered region
    ax[0].add_patch(patches.Rectangle((48, 98), 20,20, linewidth = 1, edgecolor = "black", facecolor = "none"))
    # High ordered region
    ax[0].add_patch(patches.Rectangle((90, 120), 20,20, linewidth = 1, edgecolor = "black", facecolor = "none"))



    second = ax[1].imshow(scd_bot, extent=edges, cmap="Spectral")
    ax[1].set_xlabel("x $\AA$")
    ax[1].set_title("Bot layer")
    divider2 = make_axes_locatable(ax[1])
    cax2 = divider2.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(second, cax = cax2)
    cbar.set_label('|SCD| $\AA$')

 .. image:: scd.png

The image shows regions where the order parameters are low (in red) and high (in blue). Visual examination of those regions shows the lipids have the following configurations:

 .. image:: image1aa.png

path_arr_arr,res=obj.KDEAnalysisSelection(select_res,Nframes=1000,zlim=15,show=False,legend=False)



Splay Angle
+++++++++++

The splay angle between lipid tails can also be projected to a 2D grid using :code:`Cumulative2D`. To do so, the user defines two vectors from the lipid head (usually a P-atom) to the last carbons of the lipid tails, respectively. The angle between these vectors is mapped and averaged over the set number of frames to get the following plot.


.. code:: python

    splay, edges = membrane.splay_matrix(lipid_list = ["DSPC", "DODMA", "POPE"],
                                            layer = "top",
                                            nbins = 150,
                                            start = 61,
                                            final = 110,
                                            step = 1)





.. code:: python

    plt.imshow(splay, extent = edges, cmap = "Spectral")
    plt.xlabel("$x [\AA]$")
    plt.ylabel("$y [\AA]$")
    plt.title("Splay angle")
    cbar = plt.colorbar()
    cbar.set_label('Angle $[\AA^2]$')

.. image:: splay_cumu.png



Voronoi2D
^^^^^^^^^^

:code:`Voronoi2D` also projects properties to a 2D grid, but using a different method. 
:code:`Voronoi2D` first constructs a Voronoi diagram using the positions of lipid head groups (typically lipid P-atoms), and mapping them into a :math:`m\times m` grid. The mapping step is done on each frame as illustrated in the figure below, and averages computed across n frames. At each step, the value of the computed property is assigned to the grid squares that correspond to the xy position of each lipid.  The output, similar to :code:`Cumulative2D`, is a matrix :math:`m \times m`, along with the edges :math:`[x_{\text{min}}, x_{\text{max}}, y_{\text{min}}, y_{\text{max}}]`.

.. image:: voronoii.png



To import :code:`Voronoi2D` type:

.. code-block:: python

    import MDAnalysis as mda
    from twodanalysis import Voronoi2D



Call the class using an :code:`mda.AtomGroup` or :code:`mda.Universe` as follows:

.. code-block:: python

    tpr = "membrane.tpr" # Replace with you own tpr or gro file
    xtc = "membrane.xtc" # Replace with you xtc file

    universe = mda.Universe(tpr,xtc) # Define a universe with the trajectories

    membrane = Voronoi2D(universe,   # load the universe
                    verbose = False, # Does not print initial information
                    )


.. note::
    If your trajectory contains water and/or ions, pass the list of lipids in the membrane by specifying :code:`lipid_list`.


Membrane Thickness
++++++++++++++++++

The user must set the number of bins, the edges, and the time interval. Additional options are available in the documentation.

.. code-block:: python

    lipids = membrane.lipid_list.copy()
    lipids.remove("CHL1")
    mat_thi, edges = membrane.voronoi_thickness(lipid_list=lipids,
                                            nbins = 150,           # nbins
                                            start = 61,   # Initial frame
                                            final = 110,  # Final Frame
                                            step = 1,     # Frames to skip
                                            )

The output is a matrix :math:`nbins\times nbins` and the edges in the form :math:`[x_{\text{min}}, x_{\text{max}}, y_{\text{min}}, y_{\text{max}}]`.

Visualize the output with :code:`plt.imshow`:

 .. code-block:: python

    import matplotlib.pyplot as plt

    plt.imshow(mat_thi, extent = edges, cmap = "Spectral")

    plt.xlabel("x $[\AA]$")
    plt.ylabel("y $[\AA]$")

    plt.title("Membrane thickness from frames 61-110")
    cbar = plt.colorbar()
    cbar.set_label('Thickness $\AA$')
    plt.show()

 .. image:: voronoi_thickness.png

Area per lipid
++++++++++++++

The area per lipid (APL) is a metric of lipid lateral packingm typically used to determine thermal equilibrium of a lipid bilayer. This code plots the Voronoi APL for a single frame, output images can be merged into a giff or short movies.


To run this analysis type:

.. code:: python

    voronoi_dict = membrane.voronoi_properties(layer = "top")


This will return a dictionary that contains the APL per residue in the top bilayer, accesible as :code:`voronoi_dict["apl"]`.

To map the Voronoi APL and compute its mean over time use:

.. code:: python

    areas, edges = membrane.voronoi_apl(layer = "top",
                                        nbins = 150,
                                        start = 61,
                                        final = 110,
                                        step = 1)



To render the plot use:

.. code:: python

    plt.imshow(areas, extent = edges, cmap = "Spectral")
    plt.xlabel("$x [\AA]$")
    plt.ylabel("$y [\AA]$")
    plt.title("Area per lipid")
    cbar = plt.colorbar()
    cbar.set_label('Area per lipid $[\AA^2]$')

.. image:: multiple_apl.png


Splay Angle
+++++++++++

:code:`Voronoi2D` can also project the splay angle between lipid tails to a 2D grid. Similar to :code:`Cumulative2D`, the user must set the two vectors that define the lipid tails. Using the :code:`Voronoi2D` protocol, the splay angle is plot for a set number of frames as follows.


.. code:: python

    splay, edges = membrane.voronoi_splay(layer = "top",
                                            nbins = 150,
                                            start = 61,
                                            final = 110,
                                            step = 1)



.. code:: python

    plt.imshow(splay, extent = edges, cmap = "Spectral")
    plt.xlabel("$x [\AA]$")
    plt.ylabel("$y [\AA]$")
    plt.title("Splay angle")
    cbar = plt.colorbar()
    cbar.set_label('Angle $[\AA^2]$')

.. image:: splay.png



PackingDefects
^^^^^^^^^^^^^^^

The membrane surface topology is highly dynamic, different lipid species and interactions with other biomolecules result in local changes that can be identified using :code:`2Danalysis` methods. Lipid packing defects analysis is used to quantify the exposure of the hydrophobic membrane core. :code:`PackingDefects`code allows efficient and robust statistical analysis of lipid packing deffects on the membrane surface. The analysis can be done for a single frame as well as for the full trajectory.

Import this class as follows:


.. code-block:: python

    import MDAnalysis as mda
    from twodanalysis import PackingDefects


Call the class using an :code:`mda.AtomGroup` or :code:`mda.Universe` as follows:

.. code-block:: python

    tpr = "membrane.tpr" # Replace with you own tpr or gro file
    xtc = "membrane.xtc" # Replace with you xtc file

    universe = mda.Universe(tpr,xtc) # Define a universe with the trajectories

    membrane = PackingDefects(universe,   # load the universe
                    verbose = False, # Does not print intial information
                    )

Single Frame
++++++++++++

To run the analysis for a single frame, set the frame number of interest and run:

.. code-block:: python

    membrane.u.trajectory[100] # Compute deffects for the 80 frame
    defects, defects_dict = membrane.packing_defects(layer = "top",         # layer to compute packing defects
                                                    periodic = True,  # edges for output
                                                    nbins = 400,            # number of bins
                                                    )



To plot and visualize the output run:

.. code-block:: python

    plt.imshow(defects, cmap = "viridis", extent = defects_dict["edges"]) # Plot defects
    plt.xlabel("x  $[\AA]$")
    plt.ylabel("y  $[\AA]$")
    plt.show()

.. image:: packing_defects.png

The following figure shows: (A) the packing deffects plot on VMD, (B) the output from :code:`PackingDefects`, and (C) the overlay of both approaches for comparison and validation

.. image:: packing1.png





Multiple Frames
+++++++++++++++

For statistical analysis of packing deffects across several frames, :code:`PackingDefects` returns a pandas dataframe and an array with the size of individual packing defects along the trajectory.

To run the analysis over n frames type:
.. code-block:: python

    data_df, numpy_sizes = membrane.packing_defects_stats(nbins = 400,
                                                      layer = "top",
                                                      periodic = True,
                                                      start = 0,
                                                      final = -1,
                                                      step=1)


To plot the distribution of packing defects areas type: 

.. code-block:: python

    unique, counts = np.unique(numpy_sizes, return_counts = True)
    probabilities = counts/counts.sum()

    plt.figure(figsize=(8, 5))
    plt.scatter(unique*defects_dict["grid_size"]*defects_dict["grid_size"], probabilities)
    plt.xlabel('Area $\AA$')
    plt.yscale('log')
    plt.ylabel('Probability')
    plt.title('Probability Distribution of Area')
    plt.axvline(x = 5, color = "black")
    plt.show()

.. image:: sizedefetc.png


