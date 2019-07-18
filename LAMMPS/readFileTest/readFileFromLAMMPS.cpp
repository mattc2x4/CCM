/*
 * readFileFromLAMMPS.cpp
 *
 *  Created on: Jun 17, 2019
 *      Author: mattcohe
 *
 *      the purpose of this file is to test the readFile method that ethan uses in lammps
 *      and to verify the values are what we expect.
 */
/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/
/*
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
#include <math.h>
void rvec_Add(double [], double []);
void rvec_ScaledAdd( double [], double c, double [] );
int main(){

  int i, j, pj, natoms;
  int start_i, end_i, flag;
  double p_vdW1, p_vdW1i;
  double powr_vdW1, powgi_vdW1;
  double tmp, r_ij, fn13, exp1, exp2;
  double Tap, dTap, dfn13, CEvd, CEclmb, de_core;
  double dr3gamij_1, dr3gamij_3;
  double e_ele, e_vdW, e_core, SMALL = 0.0001;
  double e_lg, de_lg, r_ij5, r_ij6, re6;


  double temp [3];
  //these were taken from
  //D:\000-Matt-Research\newTests\caviness-mpi-4activeCandO-AnirFF-NewPython-testLammps-resultsWithNewCompilation\caviness-mpi-test9-300K-4activeCandO-AnirFF-NewPython-testLammps/rest-data
  //some randoms thrown in
  //values from most recent search
  int atomID[10] = {44,39,103,20,8,54,98,49,14,9};
  //[atom][x,y,z],
  //taken fwrom D:\000-Matt-Research\newTests\caviness-mpi-4activeCandO-AnirFF-NewPython-testLammps-resultsWithNewCompilation\caviness-mpi-test9-300K-4activeCandO-AnirFF-NewPython-testLammps/coord.txt
  //second timestep
  //TO-DO: get a simulation that prints the co
  double atomX[10][3] = {{13.5921967813,15.7346236032,9.06225696455},{11.7283171,6.53646606893,13.0772211394},{14.822720201,3.9696829713,14.1422560945},{14.142256094,15.244162851,5.18390515099},{14.9966880196,15.9580134837,7.06705922518},{14.1558199751,3.21334159215,14.5178846876},{15.3442837022 ,8.05884301693 ,14.9966880196},{3.60498603696,14.2173052574,14.602216034},{15.7846341654 ,14.5509467137,5.37967313358},{15.9580134837,7.06705922518,15.7288751389  }};
  double f[10][3] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}};



  // Tallying variables:
  double pe_vdw, f_tmp, delij[3];

  natoms = 5;
  p_vdW1i = 1.0 / p_vdW1;
  e_core = 0;
  e_vdW = 0;
  e_lg = de_lg = 0.0;
  /*-------------------------------- read in restraint data --------------------------------*/
/*
  double F1rest, F2rest, R12rest, Erest, dErestdr;
  int id1rest, id2rest;
  double Rijrest, dxrest, dyrest, dzrest, Rminrest, Rmaxrest, dErestx, dEresty, dErestz;
  int nRowsrest;
  double** restMatrix; // pointer to pointer (dynamic memory allocation); is deallocated at end of file
  //std::ofstream testFile;		//this is where we initialize our test file, which we will be writing to
  //testFile.open("nonbonded_TEST.txt", std::ios_base::app);
  // read in the restraint parameter data file
  std::ifstream restfile ("rest-data.txt");

  if (restfile.is_open()) {
	  restfile >> nRowsrest;
	  // allocate
	  restMatrix = new double*[nRowsrest];
	  for (i = 0; i < nRowsrest; i++) {
		  restMatrix[i] = new double[7];
	  }
	  // fill 2d array with infile parameters
	  for (i = 0; i < nRowsrest; i++) {
		  for (j = 0; j < 7; j++) {
			  restfile >> restMatrix[i][j];
		  }
	  }
	  restfile.close();
	  for (i = 0; i < nRowsrest; i++) {
	  		for (j = 0; j < 7; j++) {
	  			cout << restMatrix[i][j] << "  ";
	  		}
	  		cout << endl;
	  	}
  }
  else {
    std::cout << "Unable to open restraint parameter file" << std::endl;
	restMatrix = new double*[1]; // wouldn't be used in this case anyway?
  }
  /*-------------------------------- end of restraint read --------------------------------*/
  /*
  for (int restk = 0; restk < nRowsrest; restk++) {
  		  id1rest = restMatrix[restk][0];
  		  id2rest = restMatrix[restk][1];
  		  R12rest = restMatrix[restk][2];
  		  F1rest = restMatrix[restk][3];
  		  F2rest = restMatrix[restk][4];
  		  Rminrest = restMatrix[restk][5];
  		  Rmaxrest = restMatrix[restk][6];
  		  //testFile <<"\n\nCACLULATING RESTRAINT FORCE";
  		  //testFile << "id1rest: "<<id1rest<< " id2rest: "<<id2rest<< " R12rest: "<<R12rest<< " F1rest: "<<F1rest<< " F2rest: "<<F2rest<< " Rminrest: "<<Rminrest<< " Rmaxrest: "<<Rmaxrest<<"\n";
  		  for ( i = 0; i < (natoms-1); i++ ) {
  			  //testFile<<"id1 compared to: "<< system->my_atoms[i].orig_id<<" = "<<(id1rest == system->my_atoms[i].orig_id)<<"   id2 compared to "<<system->my_atoms[i].orig_id<<" = "<<(id2rest == system->my_atoms[i].orig_id)<<"\n";
  			  if ( (id1rest == atomID[i]) or (id2rest == atomID[i]) ) {
  				  for ( j = (i+1); j < natoms; j++) {
  					  //testFile<<"id1 compared to: "<< system->my_atoms[j].orig_id<<" = "<<(id1rest == system->my_atoms[j].orig_id)<<"   id2 compared to "<<system->my_atoms[j].orig_id<<" = "<<(id2rest == system->my_atoms[j].orig_id)<<"\n";
  					  if ( (id1rest == atomID[j]) or (id2rest == atomID[j]) ) {
  						  dxrest = atomX[i][0] - atomX[j][0];
  						  dyrest = atomX[i][1] - atomX[j][1];
  						  dzrest = atomX[i][2] - atomX[j][2];
  						  Rijrest = sqrt(dxrest*dxrest + dyrest*dyrest + dzrest*dzrest);
  						  //testFile<< "dxrest: "<<system->my_atoms[i].x[0]<<" - "<<system->my_atoms[j].x[0]<<" = "<<dxrest<<"\n";
  						  //testFile<< "dyrest: "<<system->my_atoms[i].x[1]<<" - "<<system->my_atoms[j].x[1]<<" = "<<dyrest<<"\n";
  						  //testFile<< "dzrest: "<<system->my_atoms[i].x[2]<<" - "<<system->my_atoms[j].x[2]<<" = "<<dzrest<<"\n";
  						  //testFile<< "Rijrest: "<<Rijrest<<"\n";
  						  //testFile<<"is Rijrest < Rmaxrest, and Rijrest > Rminrest "<<((Rijrest < Rmaxrest) and (Rijrest > Rminrest));
  						  if ((Rijrest < Rmaxrest) and (Rijrest > Rminrest)) {
  							  Erest = F1rest*(1.0 - exp(-1.0*F2rest*(Rijrest - R12rest)*(Rijrest - R12rest) ));
  							  dErestdr = 2.0*F1rest*F2rest*(Rijrest - R12rest)*exp(-1.0*F2rest*(Rijrest - R12rest)*(Rijrest - R12rest));
  							  dErestx = dErestdr*dxrest/Rijrest;
  							  dEresty = dErestdr*dyrest/Rijrest;
  							  dErestz = dErestdr*dzrest/Rijrest;
  							  //testFile<< "Erest: "<<Erest<<"\n";
  							  //testFile<< "dErestdr: "<<dErestdr<<"\n";
  							  //testFile<< "dErestx: "<<dErestx<<"\n";
  							  //testFile<< "dEresty: "<<dEresty<<"\n";
  							  //testFile<< "dErestz: "<<dErestz<<"\n";
  							  temp[0] = dErestx;
  							  temp[1] = dEresty;
  							  temp[2] = dErestz;
  							  rvec_Add( f[i], temp);
  							  rvec_ScaledAdd( f[j], -1.0 ,temp);
  							  //testFile<<"dErestx, dEresty, dErestz added to workspace->f["<< i<<"]";
  						  }
  					  }
  				  }
  			  }
  		  }
  	  }
}

  //testFile.close();
  /* End of added code for restraint force */
/*
void rvec_Add(  double f[], double temp[] ){
    f[0] += temp[0], f[1] += temp[1], f[2] += temp[2];
}
void rvec_ScaledAdd( double ret[], double c, double v[] ){
  ret[0] += c * v[0], ret[1] += c * v[1], ret[2] += c * v[2];
}


*/


