/*
 * Atom.cpp
 *
 *  Created on: Jan 30, 2020
 *      Author: mattcohe
 */
#include "Atom.hpp"
using namespace std;

Atom::Atom(int idVal, double xloc, double yloc, double zloc, int atomtype){
	id = idVal;
	x = xloc;
	y = yloc;
	z = zloc;
	type = atomtype;
}



