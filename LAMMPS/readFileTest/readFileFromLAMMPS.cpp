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

	//this portion used to test file reading in reaxc_nonbonded.cpp

	/*
	double F1rest, F2rest, R12rest, Erest, dErestdr;
	int i, j, pj, natoms;
	int id1rest, id2rest;
	double Rijrest, dxrest, dyrest, dzrest, Rminrest, Rmaxrest, dErestx,
			dEresty, dErestz;
	int nRowsrest;
	double** restMatrix; // pointer to pointer (dynamic memory allocation); is deallocated at end of file
	// read in the restraint parameter data file
	std::ifstream restfile("rest-data.txt");
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
	} else {
		std::cout << "Unable to open restraint parameter file" << std::endl;
		restMatrix = new double*[1]; // wouldn't be used in this case anyway?
	}
	for (int q = 0; q < nRowsrest; q++) {
		for (int w = 0; w < 7; w++) {
			cout << restMatrix[q][w] << "  ";
		}
		cout << endl;
	}
	for (int q = 0; q < nRowsrest; q++) {
		delete[] restMatrix[q];
	}
	delete[] restMatrix;
	*/

	//this portion used to test file appending, which will be used in 00_reaxc_nonbonded_TEST.cpp


	ofstream testFile;
	testFile.open("testWrite.txt", std::ios_base::app);
	testFile << "Data";
	testFile << (1>0);
	return 0;

}
