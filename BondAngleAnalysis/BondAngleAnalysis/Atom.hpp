/*
 * Atom.hpp
 *
 *  Created on: Jan 30, 2020
 *      Author: mattcohe
 */

#ifndef ATOM_HPP_
#define ATOM_HPP_
#include <string>
using namespace std;
class Atom{
	friend class Angle;
	double x;
	double y;
	double z;
	int id;
	int type;

	Atom(int idVal, double xloc, double yloc, double zloc, int atomtype);
};





#endif /* ATOM_HPP_ */
