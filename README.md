# pFemView
An open-source C++ library to visualize p-hierarchical basis finite element (p-FEM) solutions on the scientific visualization application ParaView.

  Author: Janitha Gunatilake, developed in 2023.

ParaView:

As p-FEM solutions does not directly match the input file format of VTK files, this library brideges that gap.
Here, the input: hierarchical vtk file and the output is a .vtk file that can be read in ParaView.

Citing the package:
If you use pFemView in your research, please cite it! In BibTeX format, the following entry can be used:

  @manual{pfemview,
    key     = {pFemView},
    author  = {Janitha Gunatilake},
    title   = {pFemView: a C++ library for visualizing p-FEM solutions on ParaView (version 1.0.0)},
    note    = {{\tt [https://github.com/janithag/pFemView}},
    year    = {2023},
  }

This might render as:

    The mpmath development team. mpmath: a Python library for arbitrary-precision floating-point arithmetic (version 1.3.0), 2023. http://mpmath.org/.
