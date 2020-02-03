/*
 * angle.cpp
 *
 *  Created on: Jan 30, 2020
 *      Author: mattcohe
 */
#include "Angle.hpp"
#include <string>
using namespace std;
#include<stdlib.h>

Angle::Angle(Atom *vertex, Atom *End1, Atom *End2){
	vert = vertex;
	end1 = End1;
	end2 = End2;
	measure = 0;
}



