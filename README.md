forward
=======

Pipelines for preparing accurate electromagnetic head models and generating EEG lead field matrices from T1 / DTI

The structural and DTI meshing is a rework of the SimNIBS tool (http://simnibs.org/). The appropriate citation is:

The shell scripts have been rewritten as Nipype pipelines for parallel processing and to make it easier to start/stop/restart the processing.

At present these pipelines depend on interfaces for Gmsh (http://geuz.org/gmsh/) and Meshfix (https://code.google.com/p/meshfix/) that can only be found in the enh/conductivity branch of my fork of Nipype (https://github.com/swederik/nipype).

Install
-------

    python setup.py install

Set environment variable FWD_DIR to the root path of the installation directory

    export FWD_DIR=/Code/forward

in your .bashrc or .zshrc file


References
==============
The citations for these tools are:

SimNIBS
-------
Windhoff, M., Opitz, A. and Thielscher, A. (2011), Electric field calculations in brain stimulation based on finite elements: An optimized processing pipeline for the generation and usage of accurate individual head models. Human Brain Mapping

Gmsh
----
C. Geuzaine and J.-F. Remacle. Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities. International Journal for Numerical Methods in Engineering 79(11), pp. 1309-1331, 2009


Meshfix
-------
M. Attene - A lightweight approach to repairing digitized polygon meshes.
The Visual Computer, 2010. (c) Springer.

Nipype
------
Gorgolewski K, Burns CD, Madison C, Clark D, Halchenko YO, Waskom ML, Ghosh SS. (2011). Nipype: a flexible, lightweight and extensible neuroimaging data processing framework in Python. Front. Neuroimform. 5:13.