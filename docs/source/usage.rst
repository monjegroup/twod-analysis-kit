Usage
=====

This software contains four main classes, `Cumulative2D`, `Voronoi2D`, `PackingDefects`, `analysis`, and `BioPolymer2D` which should
be imported as follows

.. code-block:: python

   from twodanalysis import Cumulative2D
   from twodanalysis import Voronoi2D
   from twodanalysis import PackingDefects
   from twodanalysis import BioPolymer2D


We split our work in three approaches: Cumulative approaches (Cumulative2D), Voronoi approaches (Voronoi2D), and Packing defects (PackingDefects).
Cumulative and Voronoi approaches intend to project properties in the 2D space using a slightly
different logic behind that is explained in [Add article, or link to quick guide]. Some of the properties
that we can currently project are:

 - Thickness
 - Order parameters (Only for Cumulative)
 - APL (Only for Voronoi)
 - Splay angle

Furthermore, our code allows the user to project further properties in an easy way.

On the other hand, Packing defects allows us to identify regions where the Hydrophobic core is exposed,
using a method implemented in PackMem. This method has been optimized and implemented in python.





BioPolymer2D
------------

To Do
