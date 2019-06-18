#include <iostream>
#include <fstream>
#include <string>
using namespace std;

/*
int main() {
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
