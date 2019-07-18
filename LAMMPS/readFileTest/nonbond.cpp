
/*
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
#include <math.h>
int main() {
	double F1rest, F2rest, R12rest, Erest, dErestdr;
	int id1rest, id2rest;
	double Rijrest, dxrest, dyrest, dzrest, Rminrest, Rmaxrest, dErestx, dEresty, dErestz;
	int nRowsrest;
	int i;
	int j;
	double** restMatrix; // pointer to pointer (dynamic memory allocation); is deallocated at end of file
	ofstream testFile;//this is where we initialize our test file, which we will be writing to
	testFile.open("nonbonded_TEST.txt", std::ios_base::app);
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
	}
	else {
		std::cout << "Unable to open restraint parameter file" << std::endl;
		restMatrix = new double*[1]; // wouldn't be used in this case anyway?
	}
	-------------------------------- end of restraint read --------------------------------
	for (int restk = 0; restk < nRowsrest; restk++) {
		id1rest = restMatrix[restk][0];
		id2rest = restMatrix[restk][1];
		R12rest = restMatrix[restk][2];
		F1rest = restMatrix[restk][3];
		F2rest = restMatrix[restk][4];
		Rminrest = restMatrix[restk][5];
		Rmaxrest = restMatrix[restk][6];
		testFile <<"\n\nCACLULATING RESTRAINT FORCE";
		testFile << "id1rest: "<<id1rest<< " id2rest: "<<id2rest<< " R12rest: "<<R12rest<< " F1rest: "<<F1rest<< " F2rest: "<<F2rest<< " Rminrest: "<<Rminrest<< " Rmaxrest: "<<Rmaxrest<<"\n";
		for ( i = 0; i < (natoms-1); i++ ) {
			testFile<<"id1 compared to: "<< system->my_atoms[i].orig_id<<" = "<<(id1rest == system->my_atoms[i].orig_id)<<"   id2 compared to "<<system->my_atoms[i].orig_id<<" = "<<(id2rest == system->my_atoms[i].orig_id)<<"\n";
			if ( (id1rest == system->my_atoms[i].orig_id) or (id2rest == system->my_atoms[i].orig_id) ) {
				for ( j = (i+1); j < natoms; j++) {
					testFile<<"id1 compared to: "<< system->my_atoms[j].orig_id<<" = "<<(id1rest == system->my_atoms[j].orig_id)<<"   id2 compared to "<<system->my_atoms[j].orig_id<<" = "<<(id2rest == system->my_atoms[j].orig_id)<<"\n";
					if ( (id1rest == system->my_atoms[j].orig_id) or (id2rest == system->my_atoms[j].orig_id) ) {
						dxrest = system->my_atoms[i].x[0] - system->my_atoms[j].x[0];
						dyrest = system->my_atoms[i].x[1] - system->my_atoms[j].x[1];
						dzrest = system->my_atoms[i].x[2] - system->my_atoms[j].x[2];
						Rijrest = sqrt(dxrest*dxrest + dyrest*dyrest + dzrest*dzrest);
						testFile<< "dxrest: "<<system->my_atoms[i].x[0]<<" - "<<system->my_atoms[j].x[0]<<" = "<<dxrest<<"\n";
						testFile<< "dyrest: "<<system->my_atoms[i].x[1]<<" - "<<system->my_atoms[j].x[1]<<" = "<<dyrest<<"\n";
						testFile<< "dzrest: "<<system->my_atoms[i].x[2]<<" - "<<system->my_atoms[j].x[2]<<" = "<<dzrest<<"\n";
						testFile<< "Rijrest: "<<Rijrest<<"\n";
						testFile<<"is Rijrest < Rmaxrest, and Rijrest > Rminrest "<<((Rijrest < Rmaxrest) and (Rijrest > Rminrest));
						if ((Rijrest < Rmaxrest) and (Rijrest > Rminrest)) {
							Erest = F1rest*(1.0 - exp(-1.0*F2rest*(Rijrest - R12rest)*(Rijrest - R12rest) ));
							dErestdr = 2.0*F1rest*F2rest*(Rijrest - R12rest)*exp(-1.0*F2rest*(Rijrest - R12rest)*(Rijrest - R12rest));
							dErestx = dErestdr*dxrest/Rijrest;
							dEresty = dErestdr*dyrest/Rijrest;
							dErestz = dErestdr*dzrest/Rijrest;
							testFile<< "Erest: "<<Erest<<"\n";
							testFile<< "dErestdr: "<<dErestdr<<"\n";
							testFile<< "dErestx: "<<dErestx<<"\n";
							testFile<< "dEresty: "<<dEresty<<"\n";
							testFile<< "dErestz: "<<dErestz<<"\n";
							temp[0] = dErestx;
							temp[1] = dEresty;
							temp[2] = dErestz;
							rvec_Add( workspace->f[i], temp);
							rvec_ScaledAdd( workspace->f[j], -1.0 ,temp);
							testFile<<"dErestx, dEresty, dErestz added to workspace->f["<< i<<"]";
						}
					}
				}
			}
		}
	}
	testFile.close();
	 End of added code for restraint force
}

*/
