Usage
=====

This software contains four main classes, `Cumulative2D`, `Voronoi2D`, `PackingDefects`, `analysis`, and `BioPolymer2D`, which should be imported as follows

.. code-block:: python

   from twodanalysis import Cumulative2D
   from twodanalysis import Voronoi2D
   from twodanalysis import PackingDefects
   from twodanalysis import BioPolymer2D

Cumulative and Voronoi approaches project membrane properties to the membrane surface plane using a slightly
different logic behind that is explained in the Quick Guide. Current membrane analysis include membrane thickness, deuterium order paramters, area-per-lipid, and lipid tail splay angle.

Biopolymer2D currently includes four routines to characterize the adsorption and confinement mechanisms of biopolymers onto surfaces:  (i) parallel and perpendicular radii of gyration; (ii) polar histogram; (iii) 2D-projected density; and (iv) H-bonds per residue/nucleotide/bead

2Danalysis allows the user to compute and project additional structural and dynamic properties to quantify molecular interactions, conpute statistical correlations from the output data, and readily integrate the results into other analysis pipelines.
