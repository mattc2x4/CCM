# -*- coding: utf-8 -*-
import math
import numpy as np
"""
Created on Tue Oct 29 15:45:57 2019

@author: mattcohe
"""

"""Constants"""
#this is where you imput the 3 points of the angle. 

vertexType = 0   #this is where you should put the type of the vertex
endType1 = 0    #the other two types you want to calculate angle of
endType2 = 0
currStep = 4080000   #put the first step here
incrSize = -1   ## of timesteps inbetween each print in dump/ bond file

"""end Constants"""


"""don't touch these"""
vertID = 0
endID1 = 0
endID2 = 0


end1V = []  #vector whose endpoint is end1, vertex is vertex stored as [x,y,z]
end2V = []

vertList = []   #lists containing ID's of each type. 
end1List = []
end2List = []

atomList = []   #this will hold all angle data types. will be in order, with the 0th atom being the atom with ID 1.

angList = []  #this will hold all angle vals calculated below, in the format [[vertID,endID1,endID2, ANGLE],...] endID1 and 
#end2ID interchangable locations. 

boxdim = [20,20,20]     #[x,y,z] lengths
#f = open("testfile.txt")

def main():
    #cosLaw([[0,0,0],[0,0,21]],[[0,0,0],[0,21,0]])
    dump = open("dump_final_SHORT.txt","r")
    dumpLine = dump.readlines()
    fillAtomList(dumpLine,currStep)
    print(len(atomList))
    print(atomList[0])
    atom1 = Atom(0,0,0,1,1)
    atom2 = Atom(0,0,1,1,1)
    atom3 = Atom(0,1,0,1,1)
    print(cosLaw(atom1,atom2,atom3))
    
def fillAtomList(dumpLine,currStep):
    #read in the atom data from timestep "currStep"
    #finds meantions of timestep in header. compares with input (master) step. if it is greater than, breaks loop, allowing
    #main code to run instead. when it finds the correct 
    #i am sorry to all my computerscience teachers, but i have to use continue and break because of this file format.
    startRead = False
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
   #Returns the angle in radians between vectors 'v1' and 'v2'::
   a = distance(vertAtom,endAtom2)
   #print(a)
   b = distance(vertAtom, endAtom1)
   #print(b)
   c = distance(endAtom1,endAtom2)
   #print(c)
   return(math.degrees(math.acos((a**2+b**2-c**2)/(2*a*b))))
   
def getAngleID():
    # This function will be called on each timestep. It will collect all angles in the timestep, and then add them to angList. this list
    #will be modified to calculate the angle later. 
    #call once to pull data as ID
    lineFile = f.readlines()
    endCount = 0
    for line in lineFile:       #loops through each line of the file
        wordList = line.split()
        if (len(wordList) > 2):
          if (wordList[0] == '#' and wordList[1] == 'Timestep' and currStep == int(wordList[2])):      # the hashtag starts the header area
            if (currStep == 0 and wordList[0] != '#'):
                if (int(wordList[0]) in vertList):
                    bondnum = int(wordList[2])       #gets number of bond this atom has.
                    for i in range(bondnum):
                        if (int(wordList[3 + i]) in end1List or int(wordList[3 + i]) in end2List):
                            endCount+=1
                            if (endCount == 1):     #this is the first one we found, store in temp
                                firstEndID = int(wordList[3+i])
                            elif(endCount == 2):    #this is the second one we found, add time, reset counter.
                                angList.append([int(wordList[0]),firstEndID,int(wordList[3 + i])])
                                endCount = 0
    f.close()

def calcAngles(currStep):
    #this function will take the info in angList, and use it to calculate angle values based on the atom data in atomList.
    #angList: [[vertID,endID1,endID2, ANGLE,timeStep],...]
    #called on angLIst once per timestep, though angList will contain all steps. may be problematic due to data sizes. 
    for ang in angList:
        if (ang[4] == currStep):
            cosLaw(atomList[ang[0] - 1], atomList[ang[1] - 1]. atomList[ang[2] - 1])
            

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
    def __str__(self):
        return ("Atom ID: " + str(self.ID) + "\n" + "TYPE: " + str(self.TYPE) + "\n" + "X: " + str(self.x) + "\n" + "Y: " + str(self.y) + "\n" + "Z: " + str(self.z) + "\n")

class Angle:
    #this is an angle data structure. This is used to store each angle's data as seen below.
     vertID = -1
     end1ID = -1
     end2ID = -1
     timeStep = -1
     angleVal = -1
    
main()