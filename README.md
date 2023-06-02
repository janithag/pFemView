# pFemView

**Introduction**

pFemView is an open-source C++ library to visualize p-hierarchical basis finite element (p-FEM) solutions on the scientific visualization application ParaView.

As p-FEM solutions does not directly support the VTK file format in ParaView, additional work is required to visualize p-FEM solutions.
This library is meant to bridege this gap. Specifically, this library reads the p-FEM solution in a pFemView Dataifle format, and 
generates a VTK file that could be input to ParaView.

**Documentation**

The documentation is found in the *docs* directory.

**Licence**

pFemView uses an LGPL licence.

**Citing pFemView**

If you use pFemView in your research, please cite it.  

Janitha Gunatilake. *pFemView*, https://github.com/janithag/pFemView, 2023. 

Following is the BibTeX format:
```
  @manual{pfemview,
    key     = {pFemView},
    author  = {Janitha Gunatilake},
    title   = {pFemView},
    note    = {https://github.com/janithag/pFemView},
    year    = {2023},
  }
```

