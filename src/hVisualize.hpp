/**
The module hVisualize reads a "hierarchical" vtk file, extracts the information of the mesh as well as the solution, 
and write a vtk file that can be visualized on ParaView via the projections: hierarchical-to-nodal projection, and p-to-h projection.
 

Author: Janitha Gunatilake
**/

#ifndef HVISUALIZE_HPP
#define HVISUALIZE_HPP

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include "basis.hpp"
#include "element.hpp"
#include "basish.hpp"

using std::cerr;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;
using std::cout;

class hVisualize
{

private:

/* Attributes read directly from the input hierarchical vtk file */

  int numberOfVertices;		// Number of vertices in the mesh
  vector < vector <double> > vertex;	// Coordinates of the vertices of each element
  int numberOfElements;		// Number of elements in the mesh
  vector <short> elemType;		// Type of each element in the mesh: 0 - hex, 1 - tet, 2 - prism, 3 - quad, 4 - tri, 5 - line
  vector <int> p;		// Polynomial order for each element in the mesh (hierarchical basis)
  int numberOfDofs;			// Total number of (global) degrees of freedom (size of the global stiffness matrix)
  vector <int> elemDof;		// Number of degrees of freedom on each element in the mesh
  vector < vector <int> > globalDof;// Global dofs of each element in the mesh
  vector <double> Sol; 		//hierarchical basis (global) FE solution 
  string SolName;


/* Attributes used to generate the ParaView compatible vtk file */

  int spatialDim;		// Spatial dimension of the domain
  int maxP;			// Highest "p" of the mesh

  vector < vector <int> > sign;	// orientationFlag of each element
  vector < vector <int> > type;	// type of the reference element for each element (for orientation flags)

  int p_l;			// Lagrange "order p" of each element
  vector < vector <double> > LNode;	// Lagrange nodes in the projected mesh.
  vector <double> LSol;		// Projected solution to Lagrange basis
  vector < vector <const element*> > type_elem;
  vector <vector <int> > Cell;	// h-refined cells
  vector < int> CellType;	// h-refined cell types
  

/* Methods */

  void readHVtk(const char hVtk_file[]);
  void setOrientationFlags();
  void setSignAndType();
  void hierarchicalToNodalProjection2D();
  void hierarchicalToNodalProjection3D();
  int pTohProjection2D();
  int pTohProjection3D();

public:
  ~hVisualize(){
};
hVisualize();
hVisualize(const char hVtk_file[]);
void writeVtk();

};

#endif
