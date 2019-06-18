/*
 * readFileFromLAMMPS.cpp
 *
 *  Created on: Jun 17, 2019
 *      Author: mattcohe
 *
 *      the purpose of this file is to test the readFile method that ethan uses in lammps
 *      and to verify the values are what we expect.
 */
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main() {
	double F1rest, F2rest, R12rest, Erest;
	int id1rest, id2rest;
	double Rijrest, dxrest, dyrest, dzrest;
	int nRowsrest;
	double dub = 0;
	double** restMatrix; // pointer to pointer (dynamic memory allocation); is deallocated at end of file
	std::ifstream restfile("rest-data.txt");
	if (restfile.is_open()) {
		restfile >> nRowsrest;		//grabs number of atoms, listed by python
		// allocate
		//cout<<nRowsrest<<endl;
		restMatrix = new double*[nRowsrest];//creates an array of arrays on the heap, so it does not go out of scope
		for (int i = 0; i < nRowsrest; i++) {
			restMatrix[i] = new double[5];//5? this is very odd.  I guess he only cares about the first five things in each row, which makes sense
		}
		// fill 2d array with infile parameters
		for (int i = 0; i < nRowsrest; i++) {
			for (int j = 0; j < 5; j++) {
				restfile >> restMatrix[i][j];
			}
		}
		restfile.close();
	} else {
		std::cout << "Unable to open restraint parameter file" << std::endl;
		restMatrix = new double*[1]; // wouldn't be used in this case anyway?
	}
	 for (int i = 0; i < nRowsrest; i++){
			  for (int j = 0; j < 5; j++){
				  cout << restMatrix[i][j] << "  ";
			  }
			  cout << endl;
		  }
	 dub = double(restMatrix[1][0]);
	 cout<<"\n"<< dub<<endl;
	 for (int i = 0; i < nRowsrest; i++) {
	 			delete restMatrix[i];
	 		}
	return 0;
}



