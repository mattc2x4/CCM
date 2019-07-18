#include <iostream>
#include <fstream>
#include <string>
using namespace std;


int main() {/*
  string line;
  int nRows, nColumns, i, j;
  nColumns = 5;
  ifstream myfile ("rest-data.txt");
  cout << "here" << endl;
  if (myfile.is_open()) {
	  myfile >> nRows;
	  cout<<nRows<<endl;
      double restMatrix [nRows][nColumns]; //Contains id1, id2, R12, F1, F2 for each restrained pair of atoms
	  for (i = 0; i < nRows;i++) {
		  for (j = 0; j < nColumns; j++) {
			  myfile >> restMatrix[i][j]; // line;
			  // cout << line << endl;
		  }
	  }
	  for (i = 0; i < nRows; i++){
		  for (j = 0; j < nColumns; j++){
			  cout << restMatrix[i][j] << "  ";
		  }
		  cout << endl;
	  }
      myfile.close();
  }
  else{
	  cout << "Unable to open file" << endl;
  }
  return 0;
}

*/

//6-27
//I used this section to test something dumb and pointless, but to prove that you cannot modify one dimension only of a 2d array, related to the force modification in
//reaxc_nonbonded, where f is one dimensional, but otherwise is two dimensional.  they must do something to turn it into a 2d array, or maybe they don't and do something else
// its very confusing but i believe that is the root of the issue
/*
	int nRows, nColumns, i, j;
	double vec [3] = {0,1,2};
	nRows = 5;
	double restMatrix [nRows][nColumns];
	nColumns = 5;
	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nColumns; j++) {
			restMatrix[i][j] = 0; // line;
			// cout << line << endl;
		}
	}
	restMatrix[1] = vec;
	for (i = 0; i < nRows; i++) {
		for (j = 0; j < nColumns; j++) {
			cout << restMatrix[i][j] << "  ";
		}
		cout << endl;
	}
*/
	// 7-18 testing for removal of double force application. its really dumb but gotta write to a file, i is timestep basically.
	//this should only add force 5 times, even though this while loop should run 10.
	int pastStep = 0;
	int currStep = 0;
	fstream doublefile;		//initialize file.  can't open because we must open in read/write
	//cout<<doublefile.is_open()<<endl;
	int i = 0;
	bool flag = true;
	while (i<5){
		//cout<<"loop"<<endl;
		cout<< "i is: "<<i<<endl;

		currStep = i;				//read current time step to currStep
		doublefile.open("doublefile.txt",ios::in);		//read past step from file
		doublefile >> pastStep;
		doublefile.close();
		cout<<"past step :"<<pastStep<<endl;
		if(currStep != pastStep){				//compare current step to past step
			cout<<"FORCE IS ADDED!!!!!!!!!!!!! !!!!!!!!!!!"<<endl;		//add forces
		}
		  //doublefile << currStep;
		if(flag){
			i++;
			flag = false;
		}
		else if (!flag){
			flag = true;
		}
		doublefile.open("doublefile.txt",ios::out);			//update value in file to be currStep
        doublefile << currStep ;
		doublefile.close();
	}


}
