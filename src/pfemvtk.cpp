/**
The module pfemview reads a hierarchical vtk file, and writes a vtk file that can be directly visualized by ParaView.
A sample input file is given for testing purposes.

Author: Janitha Gunatilake
**/

#include "hVisualize.hpp"
int main(int argc, char **argv){
  hVisualize vtkFile("../input/Hybrid3Dhierarchicalp8.vtk");
  vtkFile.writeVtk();
}
