/**
The module "element" generates the global projection matrix needed for p-hierarchical to nodal projection in visualization.
Thus, the Lagrange nodes of order p, and the entries of the projection matrix are generated, enabling calculation of the approximate nodal (Lagrange) solution for the problem.

Author: Janitha Gunatilake
**/


#include "element.hpp"
#include <cmath>

const short unsigned typeSize[6] = {2,2,6,1,1,1};
element::~element(){

  for(int i=0;i<numOfNodes;i++)
    delete[] I[i];
  delete [] I;

  for(int i=0;i<numOfShapeFuns;i++)
    delete[] IND[i];
  delete [] IND;

  for(int i=0;i<typeSize[solid];i++)
    delete pt_h[i];

  delete pt_lag;

};

element::element(const int solid_, const int p, const int orderOfMapping){ //1-based index

  hierarchicP = p;
  if(p<1) error("basis function order p < 1.");
  solid=solid_;

  basish* pt_[6];
  pt_[0] = new hexh(p,0);
  pt_[1] = new teth(p,0);
  pt_[2] = new wedgeh(p,0);
  pt_[3] = new quadh(p);
  pt_[4] = new trih(p);
  pt_[5] = new lineh(p);

  pt_h[0] = pt_[solid];

  if(solid==0) pt_h[1] = new hexh(p,1);
  if(solid==1) pt_h[1] = new teth(p,1);
  if(solid==2) {
    for(int i=1;i<6;i++)
    pt_h[i] = new wedgeh(p,i);
  }

  for(int i=0;i<6;i++){
	if(i != solid) delete pt_[i];
	}

  numOfShapeFuns = pt_h[0]->dimOfBasis();
  dim = pt_h[0]->spatialDim();

  switch(solid){

  case 0: //hex
    projectionMatrix = &element::ProjectionMatrix3D;
    break;

  case 1: //teth
    projectionMatrix = &element::ProjectionMatrix3D;
    break;

  case 2: //wedgeh
    projectionMatrix = &element::ProjectionMatrix3D;
    break;

  case 3: //quadh
    projectionMatrix = &element::ProjectionMatrix2D;
    break;

  case 4: //trih
    projectionMatrix = &element::ProjectionMatrix2D;
    break;

  case 5: //lineh
    projectionMatrix = &element::ProjectionMatrix2D;
    break;

  default: error("Invalid input for solid!");
  }

  IND = new int* [numOfShapeFuns];
  for(int i=0;i<numOfShapeFuns;i++)
    IND[i]=new int[3];

  int* pt = pt_h[0]->getModes();
  int count=0;

  /* Element DOF map: I[type][location][index] listing to I[count]  */
  // order: 1-vertices, 2-edges, 3-quad faces (if any), 4-tri faces (if any), 5- interior
  int interior=0;
  int quadFace=0;
  int triFace=0;

  for(int i=0; i < pt[0]; i++){ //vertex modes
    IND[count][0]=0;
    IND[count][1]=i;
    IND[count][2]=0;
    count++;
  }

  if (dim != 3){ // for 1D, 2D elements.
  for(int i=0; i<p-1;i++){ // p: 1-based
    for(int location=0;location<pt[1];location++){
      IND[count][0]=1;
      IND[count][1]=location;
      IND[count][2]=i;
      count++;
    }
    if (solid == 5) continue;

    if(i>1) //interior modes for tri/quad
      for(int j=0; j<i-1; j++){
	IND[count][0]=2;
	IND[count][1]=0;
	IND[count][2]=interior;
	count++;interior++;
      }

    if((i>0)&&(solid==4)){ //additional tri interior mode (if element=triangle)
      IND[count][0]=2;
      IND[count][1]=0;
      IND[count][2]=interior;
      count++; interior++;
    }
  }
  }

  else { 
  for(int i=0; i<p-1;i++){
    for(int location=0;location<pt[1];location++){//edge modes
      IND[count][0]=1;
      IND[count][1]=location;
      IND[count][2]=i;
      count++;
    }

    if(i>1) // quad + tri face modes
    for(int j=0;j<i-1;j++){
    for(int location=0; location<pt[2]; location++){
      IND[count][0]=2;
      IND[count][1]=location;
      IND[count][2]=quadFace;
      count++;
    }
    quadFace++;

    for(int location=pt[2]; location<pt[2]+pt[3]; location++){
      IND[count][0]=2;
      IND[count][1]=location;
      IND[count][2]=triFace;
      count++;
    }
    triFace++;

    }

  if(i>0){//additional tri face mode
      for(int location=pt[2]; location<pt[2]+pt[3]; location++){
	IND[count][0]=2;
	IND[count][1]=location;
	IND[count][2]=triFace;
	count++;
      }triFace++;
  }

    if((i>1)&&(solid==1)) // internal modes for tet
    for(int j=0;j<(i-1)*i/2;j++){
      IND[count][0]=3;
      IND[count][1]=0;
      IND[count][2]=interior;
      count++;interior++;
    }

    else if((i>2)&&(solid==2)) // internal modes for wedge
    for(int j=0;j<(i-1)*(i-2)/2;j++){
      IND[count][0]=3;
      IND[count][1]=0;
      IND[count][2]=interior;
      count++;interior++;
    }

    else if((i>3)&&(solid==0)) // internal modes for hex
    for(int j=0;j<(i-2)*(i-3)/2;j++){
      IND[count][0]=3;
      IND[count][1]=0;
      IND[count][2]=interior;
      count++;interior++;
    }

  }
  }


  /* Lagrage basis: Only affine and quadratic mappings are allowed */
  int p1 = orderOfMapping;
  if(p1>2)error("Order of mapping should be linear or quadratic only!");
  basis* ptLag[6];
  ptLag[0] = new hex();
  ptLag[1] = new tet();
  ptLag[2] = new wedge();
  ptLag[3] = new quad();
  ptLag[4] = new tri();
  ptLag[5] = new line();

  pt_lag = ptLag[solid];

  for(int i=0;i<6;i++)
    if(i!=solid)
      delete ptLag[i];

  numOfNodes=pt_lag->ElementDOF(p1);

  /* Create indices for Lagrange basis upto order 2 */
  I = new int* [numOfNodes];
  for(int i=0;i<numOfNodes;i++)
    I[i]=new int[3];

  count = 0;

  /* Initialize indices to zero */

  for(int k=0;k<numOfNodes;k++){
    I[count][0]=0;
    I[count][1]=0;
    I[count][2]=0;
    count++;
  }

  count=0;
  switch(solid){
  case 0://hex
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
	I[count][0]=p1*((2+j+i)%2);
	I[count][1]=i*p1;
	count++;
      }
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
	I[count][0]=p1*((2+j+i)%2);
	I[count][1]=i*p1;
	I[count][2]=1;
	count++;
     }
     break;
  case 1:break;//tet
  case 2:break;//wedge

  case 3: 
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++){
	I[count][0]=p1*((2+j+i)%2);
	I[count][1]=i*p1;
	count++;
      }
    if(p1==1) break;
    else
      for(int i=0;i<2;i++)
	for(int j=0;j<2;j++){
	  I[count][0]=(4+j-i)%(3-i);
	  I[count][1]=(2+i+j)%2+i;
	  count++;
	}
    I[count][0]=1;I[count][1]=1; break;

  case 4: 
    for(int i=0;i<2;i++)
      for(int j=0;j<2-i;j++){
	I[count][0]=j*p1;
	I[count][1]=i*p1;
	count++;
      }
    if(p1==1) break;
    else
      for(int i=0;i<2;i++)
	for(int j=0;j<=i;j++){
	  I[count][0]=(j+1)%2;
	  I[count][1]=i;
	  count++;
	}
    break;

  case 5: 
    for(int i=0; i<=p1; i++){
	I[(p1-i+1)%(p1+1)][0]=i;
	}
    break;
  }

}

void element::ProjectionMatrix2D(const int p_, const int* sign, vector <vector <double> > vertex, double** P, vector <vector <double> > & node, const int* type) const{

  int numOfVertices = pt_lag->ElementDOF(1);

  double** x_; // refined Lagrange nodes of the element

  /* Build the local projection matrix */
  int lagrangeNodes = pt_lag->ElementDOF(p_);

  x_ = new double* [lagrangeNodes];
  for(int i=0;i<lagrangeNodes;i++)
    x_[i] = new double[3];
  pt_lag->RefinedNodes(p_,x_);

  /* Build the global projection matrix */
  for(int i=0;i<lagrangeNodes;i++)
    for(int j=0;j<numOfShapeFuns;j++)
      P[i][j]=sign[j]*pt_h[0]->eval_phi(IND[j],x_[i],hierarchicP);

  /* Global nodes of the refined mesh */ 
  for(int i=0;i<lagrangeNodes;i++)
    for(int coordinate=0;coordinate<3;coordinate++){
      node[i][coordinate]=0.;
      for(int k=0;k<numOfVertices;k++)
	node[i][coordinate] += vertex[k][coordinate]*pt_h[0]->eval_phi(IND[k],x_[i],hierarchicP);
    }

  for(int i=0;i<lagrangeNodes;i++)
      delete[] x_[i];

  delete[] x_;
}

/* 3D visualization is done by projecting to Lagrange hexahedral mesh */
void element::ProjectionMatrix3D(const int p_, const int* sign, vector <vector <double> > vertex, double** P, vector <vector <double> > & node, const int* type) const{

  hex hex1;
  int numOfVertices = hex1.ElementDOF(1); 
  double xi,eta,zeta;

  double** x_; // refined Lagrange nodes of hex
  double** X_; // refined Lagrange nodes mapped to hierarchic tet
  vector < vector < double > > v; // vertices of the spatial hexahedral element.

  /* Build the local projection matrix */
  int lagrangeNodes = hex1.ElementDOF(p_);

  x_ = new double* [lagrangeNodes]; X_ = new double* [lagrangeNodes];
  for(int i=0;i<lagrangeNodes;i++){
    x_[i] = new double[3];
    X_[i] = new double[3];
  }
  hex1.RefinedNodes(p_,x_);

  for(int i=0;i<lagrangeNodes;i++)
    for(int j=0;j<3;j++)
      X_[i][j]=x_[i][j];

  /* Map hexahedral Lagrange nodes to hierarchic reference tetrahedron */
  if(solid==1)
  for(int i=0;i<lagrangeNodes;i++){

    // hex to unit tet
    xi = (1+x_[i][0])*(1-x_[i][1])*(1-x_[i][2])/8.;
    eta = (1+x_[i][1])*(1-x_[i][2])/4.;
    zeta = (1+x_[i][2])/2.;

    // unit tet to hierarchic tet
    X_[i][0] = 2*xi+eta+zeta-1;
    X_[i][1] = sqrt(3)*eta+zeta/sqrt(3);
    X_[i][2] = sqrt(8./3)*zeta;

  }

  if(solid==2)
  for(int i=0;i<lagrangeNodes;i++){

    // hex to unit wedge
    xi = (1+x_[i][0])*(1-x_[i][1])/4.;
    eta = (1+x_[i][1])/2.;
    zeta = x_[i][2];

    // unit wedge to hierarchic wedge
    X_[i][0] = 2*xi+eta-1;
    X_[i][1] = sqrt(3)*eta;
    X_[i][2] = zeta;

  }

  /* Build the global projection matrix */
  for(int i=0;i<lagrangeNodes;i++)
    for(int j=0;j<numOfShapeFuns;j++)
      P[i][j]=sign[j]*pt_h[type[j]]->eval_phi(IND[j],X_[i],hierarchicP);

  /* Global nodes of the refined mesh */ 
  v.resize(numOfVertices);
  for(int i=0;i<numOfVertices;i++)
    v[i].resize(3);

  for(unsigned i=0;i<vertex.size();i++)
    for(int j=0;j<3;j++)
      v[i][j]=vertex[i][j];

  if(solid==1)
  for(int j=0;j<3;j++){
    v[0][j]=vertex[0][j];
    v[1][j]=vertex[1][j];
    v[2][j]=vertex[2][j];
    v[3][j]=vertex[2][j];
    v[4][j]=vertex[3][j];
    v[5][j]=vertex[3][j];
    v[6][j]=vertex[3][j];
    v[7][j]=vertex[3][j];
  }

  if(solid==2)
  for(int j=0;j<3;j++){
    v[0][j]=vertex[0][j];
    v[1][j]=vertex[1][j];
    v[2][j]=vertex[2][j];
    v[3][j]=vertex[2][j];
    v[4][j]=vertex[3][j];
    v[5][j]=vertex[4][j];
    v[6][j]=vertex[5][j];
    v[7][j]=vertex[5][j];
  }

 int ind[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};

  for(int i=0;i<lagrangeNodes;i++)
    for(int coordinate=0;coordinate<3;coordinate++){
     node[i][coordinate]=0.;
      for(int k=0;k<numOfVertices;k++)
	node[i][coordinate] += v[k][coordinate]*hex1.eval_phi(ind[k],x_[i],1);
    }

  for(int i=0;i<lagrangeNodes;i++){
      delete[] x_[i];
      delete[] X_[i];
  }
  delete[] X_;
  delete[] x_;
}
void element::error(const string& msg) const{
  throw std::runtime_error(msg);
}


