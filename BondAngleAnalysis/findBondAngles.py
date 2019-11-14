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

"""end Constants"""


"""don't touch these"""
vertID = 0
endID1 = 0
endID2 = 0
currStep = -1

end1V = []  #vector whose endpoint is end1, vertex is vertex stored as [x,y,z]
end2V = []

vertList = []   #lists containing ID's of each type. 
end1List = []
end2List = []

angList = []  #this will hold all angle vals calculated below, in the format [[vertID,endID1,endID2, ANGLE],...] endID1 and 
#end2ID interchangable locations. 

boxdim = [20,20,20]     #[x,y,z] lengths
#f = open("testfile.txt")

def main():
    cosLaw([[0,0,0],[0,0,1]],[[0,0,0],[0,1,0]])
    


def distance(end1V, end2V):
    # this function is used to calculate the distance between 2 atoms.
    #boxdim has x, y, z dimensions of box as array, 0 = x, 1 = y, 2 = z
    #see reac.f "dista2" function
    dx = end1V[0] - end2V[0]
    dy = end1V[1] - end2V[1]
    dz = end1V[2] - end2V[2]
    dx = dx - round(dx/boxdim[0]) * boxdim[0]
    dy = dy - round(dy/boxdim[1]) * boxdim[1]
    dz = dz - round(dz/boxdim[2]) * boxdim[2]
    dr = (dx*dx + dy*dy + dz*dz)**(0.5)
    return dr

def cosLaw(v1, v2):
   #v1,v2 in the format of [[start point], [end point]]
   #Returns the angle in radians between vectors 'v1' and 'v2'::
   a = distance(v1[0],v1[1])
   print(a)
   b = distance(v2[0],v2[1])
   print(b)
   c = distance(v1[1],v2[1])
   print(c)
   print(math.degrees(math.acos((a**2+b**2-c**2)/(2*a*b))))
   
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
#                
            if (currStep == 0 and wordList[0] != '#'):
                add = True
                if (int(wordList[0]) in vertList and add):
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


    
main()