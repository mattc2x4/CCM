/*
 * findBondAngles.cpp
 *
 *  Created on: Jan 30, 2020
 *      Author: mattcohe
 */
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>




struct Angle{
	struct Atom *end1;
	struct Atom *end2;
	struct Atom *vert;
	double measure;
};

struct Atom{
	double x;
	double y;
	double z;
	int id;
	int type;
};
