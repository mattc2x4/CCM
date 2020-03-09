/*
 * findBondAngles.cpp
 *
 *  Created on: Jan 30, 2020
 *      Author: mattcohe
 */
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Angle.hpp"
#include "Atom.hpp"

//http://www.open-std.org/jtc1/sc22/wg21/docs/TR18015.pdf

void getSimData(int &firstFrame, int &lastFrame, int &step, const char* dump){
	//take in vars: first frame, last frame, step size, modify by call by reference. will modify values.

	//would reading it in as a vector speed things up?
	bool readStep = false;
	int secondFrame;
	int stepFound = 0;
	int currFrame = 0;
	FILE* fp = fopen(dump, "r");
	if (fp == NULL)
	    exit(EXIT_FAILURE);
	char *currWord;
	char* line = NULL;
	size_t len = 0;
	//int i = 0;
	while ((getline(&line, &len, fp)) != -1) {
	    currWord = strtok(line, " ");
	    //cout<<"currWord: "<<currWord<<endl;
	    if(strcmp(currWord, "ITEM:") == 0){
	    	//cout<<"found ITEM header"<<endl;
	    	currWord = strtok(NULL, " ");
	    	//cout<<currWord<<endl;
	    	if (strcmp(currWord,"TIMESTEP\n") == 0){
	    		cout<<"found timestep header"<<endl;
	    		stepFound++;
	    		readStep = true;
	    		continue;
	    	}
	    }
	    if(readStep && (stepFound == 1)){
	    	sscanf(currWord, "%d", &firstFrame);
	    	readStep = false;
	    	cout<<"first step found"<<endl;
	    }
	    else if(readStep && stepFound == 2){
	    	sscanf(currWord, "%d", &secondFrame);
	    	step = secondFrame - firstFrame;
	    	readStep = false;
	    	cout<<"second step found"<<endl;
	    }
	    else if (readStep){
	    	cout<<"step found"<< endl;
	    	sscanf(currWord, "%d", &currFrame);
	    	readStep = false;
	    }
	    free(currWord);
	    //i++;
	}
	fclose(fp);
	if (line)
	    free(line);
	if(currWord){
		free(currWord);
	}
	lastFrame = currFrame;
	cout<<"done sim read"<< endl;
}

int main(){
	string dump = "dump_final.lammps";
	string bonds = "MD_bonds_final.reaxc";
	string markedXYZ = "marked_dump.xyz";
	string markedlammps = "marked_dump.lammps";
	const char* dumpChar = dump.c_str();
	int firstFrame = 0;
	int lastFrame = 0;
	int step = 0;
	getSimData(firstFrame, lastFrame, step, dumpChar);
	printf("Getting Sim Data: \n\tfirst frame: %d\n\tLast frame: %d\n\tStep: %d", firstFrame,lastFrame,step);

	return 0;
}


