Usage
=====

This software contains two main classes, `Memb2D` and `BioPolymer2D` which should
be imported as follows

.. code-block:: python

   from twodanalysis import Memb2D
   from twodanalysis import BioPolymer2D


Memb2D
------

Memb2D is a class focused on the analysis of membranes in a 2D fashion. We split
our work in three approaches: Cumulative approaches, Voronoi approaches, and Packing defects.
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
