/**
This module performs the following tasks for the Lagrange basis reference hexahedron

Author: Janitha Gunatilake
**/


#include "basis.hpp"

// Evaluates the function value of the basis function indexed I at a point X
double hex::eval_phi(const int* I,const double* X, const int n) const{
  return evalPhi(I[0], X[0], n)*evalPhi(I[1], X[1], n)*evalPhi(I[2], X[2], n);
}

// Evaluates the partial derivatives of the basis function indexed I at a point X

double hex::eval_dphidx(const int* I,const double* X, const int n) const{
  return evalDphi(I[0], X[0], n)*evalPhi(I[1], X[1], n)*evalPhi(I[2], X[2], n);
}

double hex::eval_dphidy(const int* I,const double* X, const int n) const{
  return evalPhi(I[0], X[0], n)*evalDphi(I[1], X[1], n)*evalPhi(I[2], X[2], n);
}

double hex::eval_dphidz(const int* I,const double* X, const int n) const{
  return evalPhi(I[0], X[0], n)*evalPhi(I[1], X[1], n)*evalDphi(I[2], X[2], n);
}

// The number of degrees of freedom of the Lagrange element for polynomial order p
int hex::ElementDOF(const int p) const{ // p is 1-based
  return pow(p+1,3);
}

// h-refinement for the reference hexahedron
void hex::RefinedNodes(const int p, double** x) const{ // x first then y, then z, bottom to top
    double h = 2./p;

  for(int i=0;i<p+1;i++)
    for(int j=0;j<p+1;j++)
      for(int k=0;k<p+1;k++){
	x[i*(p+1)*(p+1)+j*(p+1)+k][0]=-1+k*h;
	x[i*(p+1)*(p+1)+j*(p+1)+k][1]=-1+j*h;
	x[i*(p+1)*(p+1)+j*(p+1)+k][2]=-1+i*h;
      }
return;
}



