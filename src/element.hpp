/**
This module generates the global projection matrix needed for p-hierarchical to nodal projection in visualization.
Thus, the Lagrange nodes of order p, and the entries of the projection matrix are generated, enabling calculation of the approximate nodal (Lagrange) solution for the problem.

Author: Janitha Gunatilake
**/

#ifndef ELEMENT_HPP
#define ELEMENT_HPP
#include "basis.hpp"
#include "basish.hpp"
#include <string>
#include <algorithm>
#include <stdexcept>

using std::string;
using std::vector;
using std::max;

class element{

public:

  element(const int solid_, const int p, const int orderOfMapping);
  ~element();

  void ProjectionMatrix2D(const int p_, const int* sign, vector <vector <double> > vertex, double** P, vector <vector <double> > & node, const int* type) const;
  void ProjectionMatrix3D(const int p_, const int* sign, vector <vector <double> > vertex, double** P, vector <vector <double> > & node, const int* type) const;
  void (element::*projectionMatrix)(const int p_, const int* sign, vector <vector <double> > vertex, double** P, vector <vector <double> > & node, const int* type) const;

  int** getIndices() const{
    return IND;
  }

  unsigned ElementDOF() const{
    return numOfShapeFuns;
  }


private:
  int numOfShapeFuns;
  int solid;
  int** IND;
  int** I; //lagrange basis indices.
  int numOfNodes;  
  int dim;
  basish* pt_h[6];
  basis* pt_lag;
  unsigned hierarchicP;

  void error(const string& msg) const;

};

#endif
