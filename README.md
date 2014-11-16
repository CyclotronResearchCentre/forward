# Overview

This project aims to simplify the preparation of accurate electromagnetic head models for EEG forward modeling.

It builds off of the seminal [SimNIBS](http://simnibs.org/) tool the field modelling of transcranial magnetic stimulation (TMS) and transcranial direct current stimulation. Human skin, skull, cerebrospinal fluid, and brain meshing pipelines have been rewritten with Nipype to ease access parallel processing and to allow users to start/stop the workflows. Conductivity tensor mapping from diffusion-weighted imaging is also included.

At present these pipelines depend on interfaces for [Gmsh](http://geuz.org/gmsh/) and [Meshfix](https://code.google.com/p/meshfix/) that can only be found in the enh/conductivity branch of my fork of [Nipype](https://github.com/swederik/nipype). 

A paper describing the tool has been accepted for publication in NeuroImage:  
E. Ziegler, S.L. Chellappa, G. Gaggioni, J.Q.M. Ly, G. Vandewalle, E. André, C. Geuzaine, C. Phillips, A finite-element reciprocity solution for EEG forward modeling with realistic individual head models. NeuroImage, in press, 2014.  
Direct link: http://www.sciencedirect.com/science/article/pii/S1053811914007307  
Permalink : http://hdl.handle.net/2268/171896

# Install

To install the package, run:

    python setup.py install

Set environment variable FWD_DIR to the root path of the installation directory by placing

    export FWD_DIR=/path/to/forward

in your .bashrc or .zshrc file.
# References
The citations for these tools are:

**SimNIBS: Simulation of Non-invasive Brain Stimulation**

Windhoff, M., Opitz, A. and Thielscher, A. (2011), Electric field calculations in brain stimulation based on finite elements: An optimized processing pipeline for the generation and usage of accurate individual head models. Human Brain Mapping

**Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities**

C. Geuzaine and J.-F. Remacle. Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities. International Journal for Numerical Methods in Engineering 79(11), pp. 1309-1331, 2009

**GetDP: a General Environment for the Treatment of Discrete Problems**

C. Geuzaine and J.-F. Remacle. Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities. International Journal for Numerical Methods in Engineering 79(11), pp. 1309-1331, 2009

**Meshfix**

M. Attene - A lightweight approach to repairing digitized polygon meshes.
The Visual Computer, 2010. (c) Springer.

**Nipype: Neuroimaging in Python - Pipelines and Interfaces**

Gorgolewski K, Burns CD, Madison C, Clark D, Halchenko YO, Waskom ML, Ghosh SS. (2011). Nipype: a flexible, lightweight and extensible neuroimaging data processing framework in Python. Front. Neuroimform. 5:13.
