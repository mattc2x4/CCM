/*
 * angle.hpp
 *
 *  Created on: Jan 30, 2020
 *      Author: mattcohe
 */

#ifndef ANGLE_HPP_
#define ANGLE_HPP_

#include <string>
using namespace std;

class Atom;

class Angle{
	friend class Atom;
public:
	Atom *end1;
	Atom *end2;
	Atom *vert;
	int measure;		//degrees

	Angle(Atom *vert);

};




#endif /* ANGLE_HPP_ */
