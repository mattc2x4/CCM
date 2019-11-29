# -*- coding: utf-8 -*-
import math
import numpy as np
"""
Created on Tue Oct 29 15:45:57 2019

@author: mattcohe
"""

"""Constants"""
#this is where you imput the 3 points of the angle. 

vertexType = 4   #this is where you should put the type of the vertex
endType1 = 3    #the other two types you want to calculate angle of
endType2 = 3
#vert = 1, end1 = 3, end2 = 5 for TEST case

currStep = 4080000   #put the first step here
incrSize = -1   ## of timesteps inbetween each print in dump/ bond file

"""end Constants"""


"""don't touch these"""


atomList = []   #this will hold all angle data types. will be in order, with the 0th atom being the atom with ID 1.

angList = []  #this will hold all angle vals calculated below, in the format [[vertID,endID1,endID2, ANGLE],...] endID1 and 
#end2ID interchangable locations. 

boxdim = [20,20,20]     #[x,y,z] lengths

def main():
    #cosLaw([[0,0,0],[0,0,21]],[[0,0,0],[0,21,0]])
    dump = open("dump_final_SHORT.txt","r")
    bonds = open("MD_bonds_TEST.txt", "r")
    
    fillAtomList(dump,currStep)
    getAngleID(bonds,currStep)
    print(angList[0])
    calcAngles(currStep)
    print(angList[0])
    
def fillAtomList(dump,currStep):
    #read in the atom data from timestep "currStep"
    #finds meantions of timestep in header. compares with input (master) step. if it is greater than, breaks loop, allowing
    #main code to run instead. when it finds the correct 
    #i am sorry to all my computerscience teachers, but i have to use continue and break because of this file format.
    startRead = False
    dumpLine = dump.readlines()
    for i in range(len(dumpLine)):
        dumpWord = dumpLine[i].split()
        if(dumpWord[0] == "ITEM:" and dumpWord[1] == "TIMESTEP"):    #read timestep, if its correct continue. otherwise, break loop.
            myStep = int(dumpLine[i+1].split()[0])
            print(myStep)
            if(myStep > currStep):
                #print("breaking at " + str(myStep))
                break
            elif(myStep == currStep):
                step = True
        if(dumpWord[0] == "ITEM:" and dumpWord[1] == "ATOMS"):    #when we see this header we want to skip this iteration, then continue on the next line, hence continue. 
            startRead = True
            continue
        if (startRead and step):
            atomList.append(Atom(float(dumpWord[3]), float(dumpWord[4]), float(dumpWord[5]), int(dumpWord[0]), int(dumpWord[1])))
            #this adds an atom object in atomList. x,y,z,ID,TYPE
            
    

def distance(atom1, atom2):
    # this function is used to calculate the distance between 2 atoms.
    #boxdim has x, y, z dimensions of box as array, 0 = x, 1 = y, 2 = z
    #see reac.f "dista2" function
    dx = atom1.x - atom2.x
    dy =  atom1.y - atom2.y
    dz = atom1.z - atom2.z
    dx = dx - round(dx/boxdim[0]) * boxdim[0]
    dy = dy - round(dy/boxdim[1]) * boxdim[1]
    dz = dz - round(dz/boxdim[2]) * boxdim[2]
    dr = (dx*dx + dy*dy + dz*dz)**(0.5)
    return dr

def cosLaw(vertAtom,endAtom1, endAtom2):
    #accepts Atom objects as input
   #Returns the angle in radians between vectors 'v1' and 'v2'::
   a = distance(vertAtom,endAtom2)
   #print(a)
   b = distance(vertAtom, endAtom1)
   #print(b)
   c = distance(endAtom1,endAtom2)
   #print(c)
   return(math.degrees(math.acos((a**2+b**2-c**2)/(2*a*b))))
   
def getAngleID(bonds,currStep):
    # This function will be called on each timestep. It will collect all angles in the timestep, and then add them to angList. this list
    #will be modified to calculate the angle later. 
    #call once to pull data as ID
    #TODO: make angle data actually contain Atom objects
    lineFile = bonds.readlines()
    read = False
    for line in lineFile:       #loops through each line of the file
        wordList = line.split()
        endCount = 0
        if (len(wordList) > 2):
            if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # this checks timestep. if timestep in file matches input, read lines.upon encountering another, break loop.
                if(currStep == int(wordList[2])):
                    read = True
                else:
                    print("breaking\n")
                    break
            if (read == True and wordList[0] != '#'):
                print(str(atomList[int(wordList[0]) - 1].TYPE) + str(atomList[int(wordList[0]) - 1].TYPE == vertexType))
                if (atomList[int(wordList[0]) - 1].TYPE == vertexType):
                    #access the corresponding atom data. subtract one to translate to index. Make sure ID is correct.  and int(wordList[0]) == atomList[int(wordList[0]) - 1].ID
                    print("vertex type found")
                    print(wordList)
                    bondnum = int(wordList[2])       #gets number of bond this atom has.
                    for i in range(bondnum):
                        try: 
                            atomList[int(wordList[3 + i]) - 1].TYPE  == endType1 or atomList[int(wordList[3 + i]) - 1].TYPE == endType2
                        except IndexError:
                            print(wordList[3 + i] + " " + wordList[3 + i])
                        if (atomList[int(wordList[3 + i]) - 1].TYPE  == endType1 or atomList[int(wordList[3 + i]) - 1].TYPE == endType2):      
                           endCount+=1
                           print("found end")
                           if (endCount == 1):     #this is the first one we found, store in temp
                                firstEndID = int(wordList[3+i])
                           elif(endCount == 2):    #this is the second one we found, add time, reset counter.
                                angList.append(Angle(int(wordList[0]),firstEndID,int(wordList[3+i]),currStep,-1))
                                endCount = 0

def calcAngles(currStep):
    #this function will take the info in angList, and use it to calculate angle values based on the atom data in atomList.
    #called on angLIst once per timestep, though angList will contain all steps. may be problematic due to data sizes. 
    #TODO: modify this to decrease compexity based on what yoon wants. this dude needs to respond faster.
    for ang in angList:
        if (ang.timestep == currStep):
            ang.angle = cosLaw(atomList[ang.vertID-1],atomList[ang.end1ID-1],atomList[ang.end2ID-1])
        
            

class Atom:
    #this is an atom data structure. This is used to store each atom's x,y,z vals, Id, and Type. default values are -1.
    def __init__(self, x,y,z,ID,TYPE):
        self.x = x
        self.y = y
        self.z = z
        self.ID = ID
        self.TYPE = TYPE
   # def printData():
       # print("Atom ID: " + str(self.ID) + "\n" + "TYPE: " + str(self.TYPE) + "\n" + "X: " + str(self.x) + "\n" + "Y: " + str(self.y) + "\n" + "Z: " + str(self.z) + "\n") 
    def __str__(self):      #printing atom object will yield the following
        return ("Atom ID: " + str(self.ID) + "\n" + "TYPE: " + str(self.TYPE) + "\n" + "X: " + str(self.x) + "\n" + "Y: " + str(self.y) + "\n" + "Z: " + str(self.z) + "\n")

class Angle:
    #this is an angle data structure. This is used to store each angle's data as seen below.
    def __init__(self, vertID,end1ID,end2ID,timestep,angle):
        self.vertID = vertID
        self.end1ID = end1ID
        self.end2ID = end2ID
        self.timestep = timestep
        self.angle = angle
    def __str__(self):      #printing angle object will yield the following
        return ("Vertex ID: " + str(self.vertID) + "\n" + "End1ID: " + str(self.end1ID) + "\n" + "End2ID: " + str(self.end2ID) + "\n" + "Timestep: " + str(self.timestep) + "\n" + "angle: " + str(self.angle) + "\n")

    
main()