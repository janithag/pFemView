/**
This module contains the functions for the Lagrange basis reference prism (wedge)

Author: Janitha Gunatilake
**/
#include "basis.hpp"

double wedge::eval_phi(const int* I,const double* X, const int n) const{
  return evalTriPhi(I[0], I[1], X[0], X[1], n)*evalPhi(I[2], X[2], n); 
}

double wedge::eval_dphidx(const int* I,const double* X, const int n) const{
  return evalTriDPhidx(I[0], I[1], X[0], X[1], n)*evalPhi(I[2], X[2], n); 
}

double wedge::eval_dphidy(const int* I,const double* X, const int n) const{
  return evalTriDPhidy(I[0], I[1], X[0], X[1], n)*evalPhi(I[2], X[2], n); 
}

double wedge::eval_dphidz(const int* I,const double* X, const int n) const{
  return evalTriPhi(I[0], I[1], X[0], X[1], n)*evalDphi(I[2], X[2], n); 
}

int wedge::ElementDOF(const int p) const{
  return pow(p+1,2)*(p+2)/2;;
}

void wedge::RefinedNodes(const int p, double** x) const{
return;
}


