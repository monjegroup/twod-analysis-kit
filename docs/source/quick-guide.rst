Quick guide
==========

Here we show a quick guide of our code and some of the plots and anslysis you can get from it.

Memb2D
------

This quick tutorial will include computation of

 - Membrane Thickness
 - 2D order paramaters
 - Packing defects
 - Voronoi APL

Import the class, and mdanalysis

.. code-block:: python

    import MDAnalysis as mda
    from twodanalysis import Memb2D



All this analysis have in common the definition of a class that must be done with a universe

.. code-block:: python
    tpr = "membrane.tpr" # Replace this with you own membrane tpr or gro file
    xtc = "membrane.xtc" # Replace this with you xtc fila

    universe = mda.Universe(tpr,xtc) # Define a universe with the trajectories

    membrane = Memb2D(universe,   # load universe
                    verbose = False, # Does not print intial information
                    add_radii = True) # Add radii (used in packing defects)


Membrane Thickness
^^^^^^^^^^^^^^^^^^

This code requires the number of bins, the edges, and the time interval that you want to use. Other options
are also available, check the documentation for mor information.

.. code-block:: python

    mat_thi, edges = membrane.thickness(50,           # nbins
                                        start = 61,   # Initial frame
                                        final = 110,  # Final Frame
                                        step = 1,     # Frames to skip
                                        )

The output is a matrix :math:`nbins\times nbins` and the edges in the form :math:`[xmin,xmax,ymin,ymax]`.

We can visualize with `plt.imshow`

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
^^^^^^^^^^^^^^^^^^^^^^

The computation of order parameters is as easy as the computation of thickness. In this case
you can also choose which layer the analysis will run (top, bot, both). Follows an example of running order parameters

.. code-block:: python

    scd_top, edges = membrane.thickness( "top",       # top layer
                                        50,           # nbins
                                        start = 61,   # Initial frame
                                        final = 110,  # Final Frame
                                        step = 1,     # Frames to skip
                                        )

    scd_bot, edges = membrane.thickness( "bot",       # bot layer
                                        50,           # nbins
                                        start = 61,   # Initial frame
                                        final = 110,  # Final Frame
                                        step = 1,     # Frames to skip
                                        )

Now we can plot the results


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

Here we highligted regions where the order parameters are low (red region) and high (blue region). From this region
the lipids looks as follows

 .. image:: image1aa.png


