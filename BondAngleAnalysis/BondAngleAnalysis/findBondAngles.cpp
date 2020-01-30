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
#include "Angle.hpp"
#include "Atom.hpp"

//http://www.open-std.org/jtc1/sc22/wg21/docs/TR18015.pdf
//checking!!

string dump = "dump_final.lammps";
string bonds = "MD_bonds_final.reaxc";
string markedXYZ = "marked_dump.xyz";
string markedlammps = "marked_dump.lammps";
