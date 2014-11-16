## Overview

This project aims to simplify the preparation of accurate electromagnetic head models for EEG forward modeling.

It builds off of the seminal [SimNIBS](http://simnibs.org/) tool the field modelling of transcranial magnetic stimulation (TMS) and transcranial direct current stimulation. Human skin, skull, cerebrospinal fluid, and brain meshing pipelines have been rewritten with Nipype to ease access parallel processing and to allow users to start/stop the workflows. Conductivity tensor mapping from diffusion-weighted imaging is also included.

At present these pipelines depend on interfaces for [Gmsh](http://geuz.org/gmsh/) and [Meshfix](https://code.google.com/p/meshfix/) that can only be found in the enh/conductivity branch of my fork of [Nipype](https://github.com/swederik/nipype). 

The reference for the tools found in this repository is:

> **A finite-element reciprocity solution for EEG forward modeling with realistic individual head models**.  
E. Ziegler, S.L. Chellappa, G. Gaggioni, J.Q.M. Ly, G. Vandewalle, E. André, C. Geuzaine<sup>1</sup>, C. Phillips<sup>1</sup>
_**NeuroImage**_. Volume 103, December 2014, Pages 542-551, ISSN 1053-8119, [http://dx.doi.org/10.1016/j.neuroimage.2014.08.056](doi:10.1016/j.neuroimage.2014.08.056).

<sup>1</sup> Contributed equally

## Install

To install the package, run:

    python setup.py install

Set environment variable FWD_DIR to the root path of the installation directory by placing

    export FWD_DIR=/path/to/forward

in your .bashrc or .zshrc file.

## Links

#### Pubmed

http://www.ncbi.nlm.nih.gov/pubmed/25204867

#### ScienceDirect

http://www.sciencedirect.com/science/article/pii/S1053811914007307

#### ORBi (University of Liège)

http://orbi.ulg.ac.be/handle/2268/171896

#### NITRC

http://www.nitrc.org/projects/forward/

## BibTex code

[Direct link to .bib file](https://github.com/CyclotronResearchCentre/forward/raw/master/citation.bib)

    @article{Ziegler2014b,
    title = "A finite-element reciprocity solution for \{EEG\} forward modeling with realistic individual head models ",
    journal = "NeuroImage ",
    volume = "103",
    number = "0",
    pages = "542 - 551",
    year = "2014",
    note = "",
    issn = "1053-8119",
    doi = "http://dx.doi.org/10.1016/j.neuroimage.2014.08.056",
    url = "http://www.sciencedirect.com/science/article/pii/S1053811914007307",
    author = "Erik Ziegler and Sarah L. Chellappa and Giulia Gaggioni and Julien Q.M. Ly and Gilles Vandewalle and Elodie André and Christophe Geuzaine and Christophe Phillips",
    keywords = "Electroencephalography",
    keywords = "\{EEG\}",
    keywords = "Forward model",
    keywords = "Diffusion ",
    abstract = "Abstract We present a finite element modeling (FEM) implementation for solving the forward problem in electroencephalography (EEG). The solution is based on Helmholtz's principle of reciprocity which allows for dramatically reduced computational time when constructing the leadfield matrix. The approach was validated using a 4-shell spherical model and shown to perform comparably with two current state-of-the-art alternatives (OpenMEEG for boundary element modeling and SimBio for finite element modeling). We applied the method to real human brain \{MRI\} data and created a model with five tissue types: white matter, gray matter, cerebrospinal fluid, skull, and scalp. By calculating conductivity tensors from diffusion-weighted \{MR\} images, we also demonstrate one of the main benefits of FEM: the ability to include anisotropic conductivities within the head model. Root-mean square deviation between the standard leadfield and the leadfield including white-matter anisotropy showed that ignoring the directional conductivity of white matter fiber tracts leads to orientation-specific errors in the forward model. Realistic head models are necessary for precise source localization in individuals. Our approach is fast, accurate, open-source and freely available online. "
    }


## References
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
