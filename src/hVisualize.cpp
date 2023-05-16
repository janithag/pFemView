/**
The module hVisualize reads a "hierarchical" vtk file, extracts the information of the mesh as well as the solution, 
and write a vtk file that can be visualized on ParaView via the projections: hierarchical-to-nodal projection, and p-to-h projection.
 

Author: Janitha Gunatilake
**/

#include "hVisualize.hpp"

hVisualize::hVisualize(const char hVtk_file[]){
  readHVtk(hVtk_file);
  setOrientationFlags();
 }

void hVisualize::readHVtk(const char hVtk_file[]){

  ifstream fin(hVtk_file);

  if(!fin.is_open()){
    cerr<<"Could not open "<<hVtk_file<<endl;
    fin.clear();
  }

  /* Read the hierarchical VTK file */
  string s;
  s="0";

  while (s.compare("VERTICES") != 0)
    fin >> s;

  fin >> numberOfVertices;

  vertex.resize(numberOfVertices);
  getline(fin,s,'\n');

  for(int i=0;i<numberOfVertices;i++){
    vertex[i].resize(3);
    fin >> vertex[i][0] >> vertex[i][1] >> vertex[i][2];
  }

  while (s.compare("ELEMENTS") != 0)
    fin >> s;

  fin >> numberOfElements;

  elemType.resize(numberOfElements);
  p.resize(numberOfElements);
  elemDof.resize(numberOfElements);
  globalDof.resize(numberOfElements);

  while (s.compare("ELEMENT_TYPES") != 0)
    fin >> s;
  getline(fin,s,'\n');

  for(int i=0;i<numberOfElements;i++)
    fin >> elemType[i];

  spatialDim = 0;
  if (elemType[0]<3)
    spatialDim = 2;
  else if (elemType[0]<5)
    spatialDim = 1;

  while (s.compare("POLYNOMIAL_ORDERS") != 0)
    fin >> s;
  getline(fin,s,'\n');

  for(int i=0;i<numberOfElements;i++)
    fin >> p[i];

  while (s.compare("GLOBAL_DOFS") != 0)
    fin >> s;
  getline(fin,s,'\n');

  for(int i=0;i<numberOfElements;i++){
    fin >> elemDof[i];
    globalDof[i].resize(elemDof[i]);

    for(int j=0;j<elemDof[i];j++)
      fin >> globalDof[i][j];
  }

  while (s.compare("SOLUTION") != 0)
    fin >> s;
  fin >> numberOfDofs;
  Sol.resize(numberOfDofs);

  while (s.compare("SCALARS") != 0)
    fin >> s;
  fin >> SolName;

  while (s.compare("LOOKUP_TABLE") != 0)
    fin >> s;
  getline(fin,s,'\n');

  for(int i=0;i<numberOfDofs;i++)
    fin >> Sol[i];

  fin.close();

  maxP = *std::max_element(p.begin(),p.end());
}

void hVisualize::setOrientationFlags()
{
    vector < vector <int> > referenceElementType;
    vector < vector <int> > o;
    
    sign.resize(numberOfElements);
    type.resize(numberOfElements);
  
   for(int iel=0;iel<numberOfElements;iel++){
    sign[iel].resize(elemDof[iel]);
    type[iel].resize(elemDof[iel]);
   }
  
  int minIndex, maxIndex;
  int nNodes;
  referenceElementType.resize(numberOfElements);
  o.resize(numberOfElements);

  // Some information needed to set face orientation
  const int iFace[6][4] = {{0,1,4,5}, {1,2,5,6}, {3,2,7,6}, {0,3,4,7}, {0,1,3,2}, {4,5,7,6}}; //hex faces
  const int wFace[5][4] = {{0,1,3,4},{1,2,4,5},{2,0,5,3},{0,1,2},{3,4,5}}; //wedge faces
  const int possibleCombinations[8][3] = {{0,1,2},{1,0,3},{2,3,0},{3,2,1},{0,2,1},{1,3,0},{2,0,3},{3,1,2}};
  const int signs[8][3] = {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1}};

  for(int iel=0;iel<numberOfElements;iel++){
    if (elemType[iel] == 0)
      referenceElementType[iel].resize(6);
    else if (elemType[iel] == 2)
      referenceElementType[iel].resize(5);
    else {
      referenceElementType[iel].resize(1);
      referenceElementType[iel][0]=0;
    }
    
    switch (elemType[iel]) {
      case 0: {// Hexahedron
	  nNodes=8;
	  o[iel].resize(24); // edge+face orientation flags
	  int globalIndex[nNodes];

	  for(int i=0;i<24;i++)
	    o[iel][i]=1;

	  for(int i=0;i<nNodes;i++)
	    globalIndex[i] = globalDof[iel][i];

	  /* Edges */
	  if(globalIndex[0] > globalIndex[1])
	    o[iel][0]=-1;

	  if(globalIndex[1] > globalIndex[2])
	    o[iel][1]=-1;

	  if(globalIndex[3] > globalIndex[2])
	    o[iel][2]=-1;

	  if(globalIndex[0] > globalIndex[3])
	    o[iel][3]=-1;

	  if(globalIndex[4] > globalIndex[5])
	    o[iel][4]=-1;

	  if(globalIndex[5] > globalIndex[6])
	    o[iel][5]=-1;

	  if(globalIndex[7] > globalIndex[6])
	    o[iel][6]=-1;

	  if(globalIndex[4] > globalIndex[7])
	    o[iel][7]=-1;

	  if(globalIndex[0] > globalIndex[4])
	    o[iel][8]=-1;

	  if(globalIndex[1] > globalIndex[5])
	    o[iel][9]=-1;

	  if(globalIndex[2] > globalIndex[6])
	    o[iel][10]=-1;

	  if(globalIndex[3] > globalIndex[7])
	    o[iel][11]=-1;

  /* Faces */

  int minIndex[6][3];
  int min;

  // Find A
  for(int i=0; i<6; i++){
    min=globalIndex[iFace[i][0]];
    minIndex[i][0]=0;
    for(int j=1;j<4;j++)
      if(globalIndex[iFace[i][j]]<min){
	minIndex[i][0]=j;
        min = globalIndex[iFace[i][j]];
      }
  }

  int k;

  for(int i=0;i<6;i++){
	k=minIndex[i][0];
	minIndex[i][1]=possibleCombinations[k][1];
	minIndex[i][2]=possibleCombinations[k][2];

	if(globalIndex[iFace[i][minIndex[i][1]]] > globalIndex[iFace[i][minIndex[i][2]]])
	  k+=4;
	o[iel][12+i*2]=signs[k][0];
	o[iel][12+i*2+1]=signs[k][1];
	referenceElementType[iel][i]=(signs[k][2]==1); // type 0: 03 = -1, type 1: 03 = +1//

  }
  
      break;
    
    }
  case 1:{ //tetrahedron
    nNodes = 10;
    int index[nNodes], temp[nNodes];
    int fIndex[4], fTemp[4]; // , adjElem[4]
    int globalIndex[nNodes];

    /* vertex global indices of the tetrahedron */
    for(int i=0;i<nNodes;i++){
      globalIndex[i] = globalDof[iel][i];
      index[i]=i;
    }

    for(int i=0;i<4;i++)
      fIndex[i]=i;

    minIndex = std::min_element(globalIndex, globalIndex + 4)-globalIndex;

    /* align local index 0 with minimum global index */
    if (minIndex==3){
      index[0]=3; //vertices
      index[1]=0;
      index[3]=1;
      index[4]=7; //edges
      index[5]=6;
      index[6]=9;
      index[7]=8;
      index[8]=4;
      index[9]=5;
      fIndex[0]=3;
      fIndex[2]=0;
      fIndex[3]=2;
    }

    else
      while(index[0]!=minIndex){
	for(int i=0;i<3;i++){
  	index[i]=(index[i]+1)%3; // vertices
	index[i+4]=(index[4+i])%3+4; //edges
	index[i+7]=(index[7+i])%3+7;
	fIndex[i+1]=(fIndex[i+1])%3+1;//faces
	}
      }

    maxIndex = std::max_element(globalIndex,globalIndex + 4)-globalIndex;
    int Index=0;
    while(index[Index]!=maxIndex) Index++;
    maxIndex=Index;

    /* align local index 3 with maximum global index by rotating the face opposite to vertex 0 */
    if (3 != maxIndex)
    {
      for(int i=0;i<nNodes;i++)
	temp[i]=index[i];

      for(int i=0;i<4;i++)
	fTemp[i]=fIndex[i];

      index[1]=temp[3];
      index[2]=temp[1];
      index[3]=temp[2];
      index[4]=temp[7];
      index[5]=temp[8];
      index[6]=temp[4];
      index[7]=temp[6];
      index[8]=temp[9];
      index[9]=temp[5];
      fIndex[0]=fTemp[1];
      fIndex[1]=fTemp[3];
      fIndex[3]=fTemp[0];
  }

  maxIndex = std::max_element(globalIndex,globalIndex + 4)-globalIndex;
  Index=0;
  while(index[Index]!=maxIndex) Index++;
    maxIndex=Index;

  if (3 != maxIndex)
  {
    for(int i=0;i<nNodes;i++)
      temp[i]=index[i];

    for(int i=0;i<4;i++)
	fTemp[i]=fIndex[i];

      index[1]=temp[3];
      index[2]=temp[1];
      index[3]=temp[2];
      index[4]=temp[7];
      index[5]=temp[8];
      index[6]=temp[4];
      index[7]=temp[6];
      index[8]=temp[9];
      index[9]=temp[5];
      fIndex[0]=fTemp[1];
      fIndex[1]=fTemp[3];
      fIndex[3]=fTemp[0];
  }

  if(globalIndex[index[1]]>globalIndex[index[2]])
    referenceElementType[iel][0]=1;

  for(int inode=0;inode<nNodes;inode++)
    globalDof[iel][inode]= globalIndex[index[inode]];
  
  break;
  }
  case 2:{//wedge
    nNodes=6;
    o[iel].resize(15);  //9 + 3x2

    for(int i=0;i<15;i++)
	o[iel][i]=1;

    int globalIndex[nNodes];
    
    for(int i=0;i<nNodes;i++) {
      globalIndex[i] = globalDof[iel][i];
    
  }
      if(globalIndex[0] > globalIndex[1])
	o[iel][0]=-1;

      if(globalIndex[1] > globalIndex[2])
	o[iel][1]=-1;

      if(globalIndex[2] > globalIndex[0])
	o[iel][2]=-1;

      if(globalIndex[3] > globalIndex[4])
	o[iel][3]=-1;

      if(globalIndex[4] > globalIndex[5])
	o[iel][4]=-1;

      if(globalIndex[5] > globalIndex[3])
	o[iel][5]=-1;

      if(globalIndex[0] > globalIndex[3])
	o[iel][6]=-1;

      if(globalIndex[1] > globalIndex[4])
	o[iel][7]=-1;

      if(globalIndex[2] > globalIndex[5])
	o[iel][8]=-1;

  /* Faces */

  int minIndex[5][3];
  int min;

  // quadrilateral faces
  for(int i=0; i<3; i++){
    min=globalIndex[wFace[i][0]];
    minIndex[i][0]=0;
    for(int j=1;j<4;j++)
      if(globalIndex[wFace[i][j]]<min){
	minIndex[i][0]=j;
        min = globalIndex[wFace[i][j]];
      }
  }

  int k;

  for(int i=0;i<3;i++){
	k=minIndex[i][0];
	minIndex[i][1]=possibleCombinations[k][1];
	minIndex[i][2]=possibleCombinations[k][2];

	if(globalIndex[wFace[i][minIndex[i][1]]] > globalIndex[wFace[i][minIndex[i][2]]])
	  k+=4;
	o[iel][9+i*2]=signs[k][0];
	o[iel][9+i*2+1]=signs[k][1];
	referenceElementType[iel][i]=(signs[k][2]==1); // type 1: 03 = +1, type 0: 03 = -1//
  }
  
  // triangular faces
  for(int i=3; i<5; i++){
    min=globalIndex[wFace[i][0]];
    minIndex[i][0]=0;
    for(int j=1;j<3;j++)
      if(globalIndex[wFace[i][j]]<min){
	minIndex[i][0]=j;
        min = globalIndex[wFace[i][j]];
      }

    referenceElementType[iel][i]=minIndex[i][0];
    if(globalIndex[wFace[i][(minIndex[i][0]+1)%3]] > globalIndex[wFace[i][(minIndex[i][0]+2)%3]])
      referenceElementType[iel][i]+=3;
    
  }

   
    break;
  }
  case 3:{ //quad
      nNodes=4;
      o[iel].resize(nNodes);
      for(int i=0;i<nNodes;i++)
	o[iel][i]=1;

      int globalIndex[nNodes];

      for(int i=0;i<nNodes;i++)
        globalIndex[i] = globalDof[iel][i];

      if(globalIndex[0] > globalIndex[1])
	o[iel][0]=-1;

      if(globalIndex[1] > globalIndex[2])
	o[iel][1]=-1;

      if(globalIndex[3] > globalIndex[2])
	o[iel][2]=-1;

      if(globalIndex[0] > globalIndex[3])
	o[iel][3]=-1;
      break;
    }
  case 4:{ //triangle
      nNodes=3;
      o[iel].resize(nNodes);
      for(int i=0;i<nNodes;i++)
	o[iel][i]=1;

      int globalIndex[nNodes];

      for(int i=0;i<nNodes;i++)
	globalIndex[i] = globalDof[iel][i];

      if(globalIndex[0] > globalIndex[1])
	o[iel][0]=-1;
      if(globalIndex[1] > globalIndex[2])
	o[iel][1]=-1;
      if(globalIndex[2] > globalIndex[0])
	o[iel][2]=-1;
      break;
  }
  case 5: break;
  }
  }


  basish* pt[maxP+1][6];

  for(int i=0; i<=maxP; i++){

    pt[i][0] = new hexh(i+1,0);
    pt[i][1] = new teth(i+1,0);
    pt[i][2] = new wedgeh(i+1,0);
    pt[i][3] = new quadh(i+1);
    pt[i][4] = new trih(i+1);
    pt[i][5] = new lineh(i+1);
  }


  
  for(int iel=0;iel<numberOfElements;iel++){
	pt[p[iel]][elemType[iel]]->eval_sign(&o[iel][0], sign[iel]);  
	pt[p[iel]][elemType[iel]]->eval_type(&referenceElementType[iel][0], type[iel]);
  }

  for(int i=0;i<=maxP;i++)
    for(int j=0;j<6;j++)
      delete pt[i][j];
}


void hVisualize::writeVtk()
{
  int count;

  if(spatialDim == 1){
    hierarchicalToNodalProjection2D();
    count = pTohProjection2D();
  }

  else{
    hierarchicalToNodalProjection3D();
    count = pTohProjection3D();
  }

  char* filename= new char[60];
  sprintf(filename,"../output/mesh.%d.vtk",maxP+1);

  std::ofstream fout;
  fout.open(filename);
  if (!fout) { 
    cout << "Output mesh file "<<filename<<" cannot be opened.\n";
    exit(0);
  }

  /* head */
  fout<<"# vtk DataFile Version 3.0\nAMR mesh\nASCII\n";
  fout<<"DATASET UNSTRUCTURED_GRID\n";

  /* nodes */
  fout<<"POINTS "<<LNode.size()<<" double\n";
  for(unsigned i=0;i<LNode.size();i++){
    fout<<LNode[i][0]<<" "<<LNode[i][1]<<" "<<LNode[i][2]<<endl;
  }

  /* cells */
  fout<<"\nCELLS "<<Cell.size()<<" "<<count<<endl; // counter: total number of integers required to represent the list.

  for(unsigned i=0;i<Cell.size();i++){
    fout<<Cell[i].size()<<" ";
    for(unsigned j=0;j<Cell[i].size();j++)
      fout<<Cell[i][j]<<" ";
    fout<<endl;
  }

  fout<<"CELL_TYPES "<<CellType.size()<<endl;

  for(unsigned i=0;i<CellType.size();i++)
    fout<<CellType[i]<<endl;

  /* solution */
  fout<<"\nPOINT_DATA "<<LNode.size()<<endl;
  fout<<"SCALARS "<<SolName<<" double 1"<<endl<<"LOOKUP_TABLE default\n";
  for(unsigned j=0;j<LSol.size();j++)
      fout<<LSol[j]<<endl;

  fout.close();
  delete [] filename;
}

 void hVisualize::hierarchicalToNodalProjection2D(){ 

 vector < vector <double> > vt; // coordinates of the vertices of each element

  int LElemDof;
  int HElemDof;
  int LSolSize = 0;
  vector < vector <double> > x; // global lagrange nodes

  basis* pt[2];

  quad quad_ = quad();
  tri tri_ = tri();

  pt[0] = &quad_;
  pt[1] = &tri_;

  for(int kel=0; kel < numberOfElements; kel++){
    LElemDof = pt[elemType[kel]-3]->ElementDOF(p[kel]+1); 
    LSolSize += LElemDof;
  }

  LSol.resize(LSolSize);
  LNode.resize(LSolSize);
  for(int i=0;i<LSolSize; i++)
    LNode[i].resize(3);

  int offset = 0;

  int QuadElemDof = pt[0]->ElementDOF(maxP+1); 
  double** phi = new double*[QuadElemDof]; // create projection matrix to the size of quad

  type_elem.resize(6);
  for(int i=0;i<6;i++){
      type_elem[i].resize(maxP+1);
      for(int j=0;j<=maxP;j++)
	type_elem[i][j] = new const element(i, j+1, 1);
    }

  for(int i=0; i < QuadElemDof; i++)
    phi[i] = new double[type_elem[3][maxP]->ElementDOF()]; //size of quad. 
  int nVertices[2]={4,3};

  for(int kel=0; kel < numberOfElements; kel++){

    LElemDof = pt[elemType[kel]-3]->ElementDOF(p[kel]+1);
    HElemDof = type_elem[elemType[kel]][p[kel]]->ElementDOF();
    x.resize(LElemDof);

    for(int i=0;i<LElemDof;i++)
      x[i].resize(3);

    vt.resize(nVertices[elemType[kel]-3]);

    for(unsigned i=0;i<vt.size();i++){
      vt[i].resize(3);
      vt[i][0]=vertex[globalDof[kel][i]][0];
      vt[i][1]=vertex[globalDof[kel][i]][1];
      vt[i][2]=vertex[globalDof[kel][i]][2];
    }

    (type_elem[elemType[kel]][p[kel]]->*type_elem[elemType[kel]][p[kel]]->projectionMatrix)(p[kel]+1, &sign[kel][0], vt, phi, x, &type[kel][0]);

    for(int i=0;i<LElemDof;i++){
      LSol[i+offset]=0.;

      for(int j=0;j<HElemDof;j++)
	LSol[i+offset] += Sol[globalDof[kel][j]]*phi[i][j];
    }

    for(int i=0;i<LElemDof;i++)
      for(int j=0;j<3;j++)
	LNode[i+offset][j] = x[i][j];

    offset += LElemDof;

  }

  for(int i=0;i<QuadElemDof;i++)
    delete[] phi[i];
  delete[] phi;

 for(int i=0; i<6; i++)
    for(int j=0; j<=maxP; j++)
      delete type_elem[i][j];

 }


void hVisualize::hierarchicalToNodalProjection3D(){ 

  int LElemDof;
  int HElemDof;
  int LSolSize=0;

  vector < vector <double> > vt; // coordinates of the vertices of each element
  vector < vector <double> > x; // global lagrange nodes

  basis* pt[3];

  hex hex_ = hex();
  tet tet_ = tet();
  wedge wedge_ = wedge();

  pt[0]=&hex_;
  pt[1]=&tet_;
  pt[2]=&wedge_;

  for(int kel=0; kel < numberOfElements; kel++){
    LElemDof = pt[0]->ElementDOF(p[kel]+1); 
    LSolSize += LElemDof;
  }

  LSol.resize(LSolSize);
  LNode.resize(LSolSize);
  for(int i=0;i<LSolSize; i++)
    LNode[i].resize(3);

  int offset=0;

  type_elem.resize(6);
  for(int i=0;i<6;i++){
      type_elem[i].resize(maxP+1);
      for(int j=0;j<=maxP;j++)
	type_elem[i][j] = new const element(i, j+1, 1);
  }

  int nVertices[3]={8,4,6};
  int HexElemDof = pt[0]->ElementDOF(maxP+1); 

  double** phi = new double*[HexElemDof]; // create projection matrix to the size of hex

  for(int i=0; i < HexElemDof; i++)
    phi[i] = new double[type_elem[0][maxP]->ElementDOF()]; //size of hexDeposits



  for(int kel=0; kel < numberOfElements; kel++){
    LElemDof = pt[0]->ElementDOF(p[kel]+1);
    HElemDof = type_elem[elemType[kel]][p[kel]]->ElementDOF();

  x.resize(LElemDof);
  for(int i=0;i<LElemDof;i++)
	x[i].resize(3);

  vt.resize(nVertices[elemType[kel]]);

  for(unsigned i=0;i<vt.size();i++){
      vt[i].resize(3);
      vt[i][0]=vertex[globalDof[kel][i]][0];
      vt[i][1]=vertex[globalDof[kel][i]][1];
      vt[i][2]=vertex[globalDof[kel][i]][2];
    }

   (type_elem[elemType[kel]][p[kel]]->*type_elem[elemType[kel]][p[kel]]->projectionMatrix)(p[kel]+1, &sign[kel][0], vt, phi, x, &type[kel][0]);

    for(int i=0;i<LElemDof;i++){
      LSol[i+offset]=0.;

      for(int j=0;j<HElemDof;j++)
	LSol[i+offset] += Sol[globalDof[kel][j]]*phi[i][j];
    }

    for(int i=0;i<LElemDof;i++)
      for(int j=0;j<3;j++)
	LNode[i+offset][j] = x[i][j];

    offset += LElemDof;
  }

  for(int i=0;i<HexElemDof;i++)
    delete[] phi[i];
  delete[] phi;

  for(int i=0; i<6; i++)
    for(int j=0; j<=maxP; j++)
      delete type_elem[i][j];

 }

int hVisualize::pTohProjection2D()
{
  // for line: Number of splitted elements per single element = (maxP)
  // for quad, tri: number of splitted elements per single element = (maxP+1)^2

  int row, col;
  int drow;
  const int cellType[6] = {12,10,13,9,5,3};
  int elemOffset = 0;
  int nodeOffset = 0;
  int rowChanged;

  int CellSize = numberOfElements*pow((maxP+1),2);
  Cell.resize(CellSize);
  CellType.resize(CellSize);

  basis* pt[2];

  quad quad_ = quad();
  tri tri_ = tri();

  pt[0] = &quad_;
  pt[1] = &tri_;
  const int nVertices[2] = {4,3};
  int nQuads=0;
  int nTris=0;

  for(int kel=0; kel < numberOfElements; kel++){

    p_l = p[kel]+1;

    if(elemType[kel]==3)
      nQuads += pow(p_l,2);
    else
      nTris += pow(p_l,2);

    row = 0; drow=0;
    for(int i=0;i<pow((p_l),2);i++){
      Cell[i+elemOffset].resize(nVertices[elemType[kel]-3]);

      if(elemType[kel] == 3){ //quadrilateral
	row = i/p_l; col=i%p_l;
	Cell[i+elemOffset][0] = row*(p_l+1)+col+nodeOffset;
	Cell[i+elemOffset][1] = row*(p_l+1)+col+1+nodeOffset;
	Cell[i+elemOffset][2] = (row+1)*(p_l+1)+col+1+nodeOffset;
	Cell[i+elemOffset][3] = (row+1)*(p_l+1)+col+nodeOffset;
      }

      else{ //triangle
	rowChanged = false;

	if (i == 0){
	  Cell[i+elemOffset][0] = 0+nodeOffset;
	  Cell[i+elemOffset][1] = 1+nodeOffset;
	  Cell[i+elemOffset][2] = 1+p_l+nodeOffset;
	}

	else if(i < p_l*(p_l+1)/2){// upright triangles
	  int sum = 0;
	  for(int k=0;k<row+1;k++)
	    sum += p_l-k;

	  if(i == sum){rowChanged=1; row++;}

	  Cell[i+elemOffset][0] = Cell[i-1+elemOffset][1]+rowChanged;
	  Cell[i+elemOffset][1] = Cell[i+elemOffset][0]+1;
	  Cell[i+elemOffset][2] = Cell[i+elemOffset][1]+p_l-row;
	}
	else {// upside down triangles

	  int k = i-(p_l)*(p_l+1)/2;
	  if (k == 0){
	    Cell[i+elemOffset][0] = 1+nodeOffset;
	    Cell[i+elemOffset][1] = 2+p_l+nodeOffset;
	    Cell[i+elemOffset][2] = Cell[i+elemOffset][1]-1;
	  }
	  else{
	    int sum = 0;
	    for(int kk=0;kk<drow+1;kk++)
	      sum += p_l-1-kk;

	    if(k == sum){rowChanged=2; drow++;}

	    Cell[i+elemOffset][0] = Cell[i-1+elemOffset][0]+1+rowChanged;
	    Cell[i+elemOffset][1] = Cell[i+elemOffset][0]+p_l-drow+1;
	    Cell[i+elemOffset][2] = Cell[i+elemOffset][1]-1;
	  }

	}
      }
      CellType[i+elemOffset] = cellType[elemType[kel]];

    }
    elemOffset += pow(p_l,2);
    nodeOffset += pt[elemType[kel]-3]->ElementDOF(p_l);
  }
   Cell.resize(elemOffset);
   CellType.resize(elemOffset);

  return nQuads*5 + nTris*4;
}

int hVisualize::pTohProjection3D()
{

  int short nVertices=8;
  int row, col, height;

  int elemOffset = 0;
  int nodeOffset = 0;
  int nHex = 0;

  int CellSize = numberOfElements*pow((maxP+1),3);
  Cell.resize(CellSize);
  CellType.resize(CellSize);

  basis* pt;
  hex hex_ = hex();
  pt = &hex_;

  for(int kel=0; kel < numberOfElements; kel++){

    p_l = p[kel]+1;
    nHex += pow(p_l,3);

    for(int i=0;i<pow((p_l),3);i++){
      Cell[i+elemOffset].resize(nVertices);

	height=i/pow(p_l,2);  row = (i-pow(p_l,2)*height)/p_l; col=i%p_l;

	Cell[i+elemOffset][0]=height*pow(p_l+1,2)+row*(p_l+1)+col+nodeOffset;
	Cell[i+elemOffset][1]=height*pow(p_l+1,2)+row*(p_l+1)+col+1+nodeOffset;
	Cell[i+elemOffset][2]=height*pow(p_l+1,2)+(row+1)*(p_l+1)+col+1+nodeOffset;
	Cell[i+elemOffset][3]=height*pow(p_l+1,2)+(row+1)*(p_l+1)+col+nodeOffset;
	Cell[i+elemOffset][4]=(height+1)*pow(p_l+1,2)+row*(p_l+1)+col+nodeOffset;
	Cell[i+elemOffset][5]=(height+1)*pow(p_l+1,2)+row*(p_l+1)+col+1+nodeOffset;
	Cell[i+elemOffset][6]=(height+1)*pow(p_l+1,2)+(row+1)*(p_l+1)+col+1+nodeOffset;
	Cell[i+elemOffset][7]=(height+1)*pow(p_l+1,2)+(row+1)*(p_l+1)+col+nodeOffset;

      CellType[i+elemOffset]=12;

    }
    elemOffset += pow((p_l),3);
    nodeOffset += pt->ElementDOF(p_l);
  }
  Cell.resize(elemOffset);
  CellType.resize(elemOffset);
  return nHex*9; // 9 = number of vertices of a linear element +1*/
}
