# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 09:59:40 2019

@author: mattcohe

THis is the most recent and fixed version.  It could still use quite a few changes
check search more
remove second for loop in search responsible for appending to the rest-ALLData.txt. 
"""

# run this python script as master for simulations requiring restraint energy.
# Use the search() function at regular intervals
# to update rest-data.txt (restrained atoms)
# from mpi4py import MPI
# GLOBAL VARIABLES




# atom types
Ctype = 5;
Otype = 6;
Ntype = 3;
Htype = 7;
excludeC = []       # if 2 C's are bonded then we should skip them in the search, since C is principal.  Added in excludeN. 
# min and max pair distances for restraint criteria
OHdist = [1.5, 8.0]
NHdist = [0.9, 1.2]
NCdist = [3.0 ,8.0]
COdist = [1.3 ,1.6]

# for bond finding portion
Clist = []
Olist = []
Nlist = []
Hlist = []
NHlist = [] # this is for checking whether or not the N and H in the groups we are applying force to are from the same molecule. 
exclude = []    #modified in excludeN, containts N ID's which have 2 active C bonds and should be excluded from force addition
#should only be accessed by things that require input as ID, or subtract one. 
#new array: 2d.  restID[0][i] will be associated with a group that is valid to recieve restraint force. written in order [C,O,N,H]
#restOD[i][0] = C, restOD[i][1] = O, restOD[i][2] = N, restOD[i][3] = H
restID = []

# energy equation parameters (F1, F2)
# F1OC = 50
# F2OC = 0.5
# F1OH = 250
# F2OH = 0.75
# F1CN = 300
# F2CN = 0.75

boxdim = [0] *3

F1OC = 75
F2OC = 1.0
F1OH = 200
F2OH = 0.75
F1CN = 300
F2CN = 0.75

F1CO = '?'
F2CO = '?'

# equilibrium distances
# R12OC = 1.95
# R12OH = 1.1
# R12CN = 1.3

# equilibrium distances
R12OC = 1.95
R12OH = 1.05
R12CN = 1.4


#distances [min,max] within which we consider the bond to be made. +/- 15% of above equilibrium distances
CObondDist = [2.3 , 2.9]
OHbondDist = [.85 , 1.2]
CNbondDist = [1.25 , 1.61]
NHbondDist = [2.1 , 3.2]

#where we put the formed bonds found in findSuccess
bondID = []


#time step
timestep = 50000





#file which will contain every single instance where we apply restraint energy
#type is A, which means that rather than writing over the file, which is the case in rest-data.txt, it will log every instance of restraint case being met.
#TO-DO: find a way to include the exact time step in each instance of restraint.
newfile = open("rest-ALLdata.txt",'a')
coordFile = open("coord.txt",'a')
bondFile = open("pyBond.txt",'a')



def main():
    
    # import libraries
    # print "got here..."
    from mpi4py import MPI
    # print "and here"
    my_rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()
    from lammps import lammps

    if my_rank == 0:
        print "begin python script"

    # user input; global variables
    lammpsScript = "in.accelerated-test.txt"
    lammpsArgs = ["-echo","log"]

    #--------------END OF USER DEFINED PARAMETERS--------------------------------------

    # run lammps script
    lmp1 = lammps(comm=MPI.COMM_WORLD,cmdargs=lammpsArgs)  #lammps(cmdargs=lammpsArgs) 
    lmp1.file(lammpsScript)
    initialize()
    lmp1.command("run 0")
    #todo, grab box dimensions for periodic boundaries.
     #- lmp1.extract_global("boxzlo",0)
    #boxlo,boxhi,xy,yz,xz,periodicity,box_change = lmp1.extract_box()
    #coordFile.write("xdim: " + str(boxdim[0]) + "\nydim: " + str(boxdim[1]) + "\nzdim: " + str(boxdim[2]) + "\n\n")
    # Get lammps data as python variables
    natoms = lmp1.get_natoms()
    coordinates = lmp1.gather_atoms("x",1,3)
    atomType = lmp1.gather_atoms("type",0,1)
    if(my_rank == 0):
        getCandN(atomType)
        getNH()
    for i in range(100):
        newfile = open("rest-ALLdata.txt",'a')
        bondFile = open("pyBond.txt",'a')
        #gets global quantity ntimestep (current timestep) in lammps.  more extractable stuff can be viewed in library.cpp
        currentStep = lmp1.extract_global("ntimestep",0)
        natoms = lmp1.get_natoms()
        boxdim[0] = lmp1.extract_global("boxxhi",1)
        boxdim[1] = lmp1.extract_global("boxyhi",1)
        boxdim[2] = lmp1.extract_global("boxzhi",1)
        coordinates = lmp1.gather_atoms("x",1,3)
        atomType = lmp1.gather_atoms("type",0,1)
        if(my_rank == 0):
            coordFile.write("Time Step: " + str(currentStep) + "\n")
            search(natoms, atomType, coordinates,currentStep)
        lmp1.command("run " + str(timestep)) # lmp1.command("run 100000")
        newfile.close()
        bondFile.close()
 
    

    # add additional commands to lammps instance and run additional steps

# =============================================================================
#     for i in range(2):
# 		natoms = lmp1.get_natoms()
# 		coordinates = lmp1.gather_atoms("x",1,3)
# 		atomType = lmp1.gather_atoms("type",0,1)
# 		coordinatesLocal = lmp1.extract_atom("x",3)
# 		coordFile.write("\n" + str(coordinates(0)) + str(coordinates(1))+ str(coordinates(2)))
#         search(natoms, atomType, coordinates)        
#         lmp1.command("run " + str(timestep)) # lmp1.command("run 100000")
#         timeCurr += timestep
# =============================================================================

    lmp1.close()

    if my_rank == 0:
        print "End of run"
        MPI.Finalize()
    # End of python script
    



# supporting functions below
#--------------------------#

def search(natoms, atomType, c,currentStep): # c = coordinates
    # this function is used to search for atom pairs and to produce a
    # data file "rest-data.txt"
    # in nested for loops: i = Carbon, j = Oxygen, k = Nitrogen, m = Hydrogen; don't confuse address with id; address = id - 1
    coordFile = open("coord.txt",'a')
    restID = []
    for i in range(0,natoms):
        if (atomType[i] == Ctype ):
            for j in range(0,natoms):
                if ((atomType[j] == Otype) and (COdist[0] < distance(i,j,c)) and (COdist[1] > distance(i,j,c))):
                    for k in range(0,natoms):
                        if ( (atomType[k] == Ntype) and (NCdist[0] < distance(i,k,c)) and (NCdist[1] > distance(i,k,c))):
                            for m in range(0,natoms):
                                if ( (atomType[m]==Htype) and (OHdist[0] < distance(j,m,c)) and (OHdist[1] > distance(j,m,c)) and (NHdist[0] < distance(k,m,c)) and (NHdist[1] > distance(k,m,c))):
                                    # +1 means address converted to id
                                    #[[C,O,N,H],...]
                                    restID.append([i+1,j+1,k+1,m+1])
                                    #coordFile.write("\ntimestep:" + str(currentStep) + "\n")
                                    #coordFile.write("restID array: " + str(restID) + "\n")
                                    newfile.write("\nATOMS FOUND: timestep: " + str(currentStep) + " \nC\n ID: " + str(i+1) + " X: " + str(c[i*3]) + " Y: " + str(c[i*3+1]) + " Z: " + str(c[i*3+2]))
                                    newfile.write("\nO\n ID: " + str(j+1) + " X: " + str(c[j*3]) + " Y: " + str(c[j*3+1]) + " Z: " + str(c[j*3+2]))
                                    newfile.write("\nN\n ID: " + str(k+1) + " X: " + str(c[k*3]) + " Y: " + str(c[k*3+1]) + " Z: " + str(c[k*3+2]))
                                    newfile.write("\nH\n ID: " + str(m+1) + " X: " + str(c[m*3]) + " Y: " + str(c[m*3+1]) + " Z: " + str(c[m*3+2]))
    difFile = open("difFile.txt",'a')
    difFile.write("Current Step: " + str(currentStep) +"\n NHlist: " + str(NHlist) + "\n")
    if (NHlist):
        for i in restID:
             if(not validGroupNH(i)):
                 difFile.write("restID: " + str(restID) + "\n")
                 difFile.write(str(i) + " has N and H from different starting groups\n")
                 difFile.write(str(i) + " removing\n")
                 restID.remove(i)
                 
             else:
                 difFile.write(str(i) + " is good\n")
    difFile.close()         
    for i in range(0,len(restID)):      #you win python, this is such a stupid way to do this, and I hate you.  Calls findOptimal for the len, and if it doesn't return -1 it deletes at the index
        delInd = findOptimal(restID, c)
        if (delInd != -1):
            coordFile.write("deleting: " + str(restID[delInd]) + "\n")
            del restID[delInd]
        else:
            coordFile.write("none inoptimal\n")
            break
            
    coordFile.write(str(currentStep) + " restID array (post removal): " + str(restID) + "\n")

    restfile = open("rest-data.txt",'w')
    restfile.write(str(len(restID) * 3))        #number of groups we would apply force to TODO: must check for distances.
    #newfile.write("\n" + "timestep: " + str(currentStep))
    newfile.write(str(len(restID) * 3))
    for i in range(0,len(restID)):
        #O,C
        #O,H
        #C,N
        
        #if (my_rank == 0):
        restfile.write("\n" + str(restID[i][1]) + " " + str(restID[i][0]) + " " + str(R12OC) + " " + str(F1OC) + " " + str(F2OC) + " " + str(COdist[0]) + " " + str(COdist[1]))
        restfile.write("\n" + str(restID[i][1]) + " " + str(restID[i][3]) + " " + str(R12OH) + " " + str(F1OH) + " " + str(F2OH) + " " + str(OHdist[0]) + " " + str(OHdist[1]))
        restfile.write("\n" + str(restID[i][0]) + " " + str(restID[i][2]) + " " + str(R12CN) + " " + str(F1CN) + " " + str(F2CN) + " " + str(NCdist[0]) + " " + str(NCdist[1]))
        newfile.write("\n" + str(restID[i][1]) + " " + str(restID[i][0]) + " " + str(R12OC) + " " + str(F1OC) + " " + str(F2OC) + " " + str(COdist[0]) + " " + str(COdist[1]))
        newfile.write("\n" + str(restID[i][1]) + " " + str(restID[i][3]) + " " + str(R12OH) + " " + str(F1OH) + " " + str(F2OH) + " " + str(OHdist[0]) + " " + str(OHdist[1]))
        newfile.write("\n" + str(restID[i][0]) + " " + str(restID[i][2]) + " " + str(R12CN) + " " + str(F1CN) + " " + str(F2CN) + " " + str(NCdist[0]) + " " + str(NCdist[1]))

        
    restfile.close()
    coordFile.close()

  #  newfile.write("\n" + "timestep: " + str(currentStep))
#    if (str(len(restID)) == 0):
#        newfile.close()
#        coordFile.close()

        
    return 0


def distance(address1, address2,coordinates):
    # this function is used to calculate the distance between 2 atoms.  Input is the address (not id!) of the two atoms
    #boxdim has x, y, z dimensions of box as array, 0 = x, 1 = y, 2 = z
    #see reac.f "dista2" function
    dx = coordinates[3*address1] - coordinates[3*address2]
    dy = coordinates[3*address1+1] - coordinates[3*address2+1]
    dz = coordinates[3*address1+2] - coordinates[3*address2+2]
    dx = dx - round(dx/boxdim[0]) * boxdim[0]
    dy = dy - round(dy/boxdim[1]) * boxdim[1]
    dz = dz - round(dz/boxdim[2]) * boxdim[2]
    dr = (dx*dx + dy*dy + dz*dz)**(0.5)
    return dr

def initialize():
    restfile = open("rest-data.txt",'w')
    restfile.write("0")
    restfile.close()
    return 0

#return the perimeter of the atoms. only input 1D array, atom group.
#distance only accepts ID's, so gotta subtract one
def getPerim(atomGroup,coord):
    return distance(atomGroup[0]-1,atomGroup[1]-1,coord) + distance(atomGroup[0]-1,atomGroup[3]-1,coord) + distance(atomGroup[1]-1,atomGroup[2]-1,coord) + distance(atomGroup[2]-1,atomGroup[3]-1,coord)

#python really got me here.  Basically call this for as many items in the list. I really should think of a better way, but python. returns an index that is inoptimal (to be deleted) or -1 ( if none are inoptimal)
def findOptimal(restID, c):
    for i in range(0,len(restID)):
        for j in range(0,len(restID)):
            if (restID[i][0] == restID[j][0] and (i != j)):       #if groups share C
                if (getPerim(restID[j],c) <= getPerim(restID[i],c)):     #if the Perimeter of the group at j is less than Perimeter distance of the group at i, delete group at i. else delete group at j
                    coordFile.write("Passing group: " + str(restID[i]) + "\n")
                    coordFile.write(str(restID[i]) + "Perim = " + str(getPerim(restID[i],c)) + "\n")
                    coordFile.write(str(restID[j]) + "Perim = " + str(getPerim(restID[j],c)) + "\n")
                    return i       #adds i to delArray, which will be deleted after completion of loops.
                else:
                    coordFile.write("Passing group: " + str(restID[j]) + "\n")
                    coordFile.write(str(restID[i]) + "Perim = " + str(getPerim(restID[i],c)) + "\n")
                    coordFile.write(str(restID[j]) + "Perim = " + str(getPerim(restID[j],c)) + "\n")
                    return j
            if (restID[i][2] == restID[j][2] and (i != j)):       #if groups share N
                 if (getPerim(restID[j],c) <= getPerim(restID[i],c)):     #if the Perimeter of the group at j is less than Perimeter distance of the group at i, delete group at i. else delete group at j
                     coordFile.write("Passing group: " + str(restID[i]) + "\n")
                     coordFile.write(str(restID[i]) + "Perim = " + str(getPerim(restID[i],c)) + "\n")
                     coordFile.write(str(restID[j]) + "Perim = " + str(getPerim(restID[j],c)) + "\n")
                     return i      #deletes array at i
                 else:
                     coordFile.write("Passing group: " + str(restID[j]) + "\n")
                     coordFile.write(str(restID[i]) + "Perim = " + str(getPerim(restID[i],c)) + "\n")
                     coordFile.write(str(restID[j]) + "Perim = " + str(getPerim(restID[j],c)) + "\n")
                     return j
       
    return -1

#This function will search for successfully created bonds
#def findSuccessBonds(natoms, atomType, c,currentStep):
#    bondFile = open("pyBond.txt",'a')
#    for i in range(0,natoms):
#        if atomType[i] == Ctype:
#            for j in range(0,natoms):
#                if ((atomType[j] == Otype) and (CObondDist[0] < distance(i,j,c)) and (CObondDist[1] > distance(i,j,c))):
#                    for k in range(0,natoms):
#                        if ( (atomType[k] == Ntype) and (CNbondDist[0] < distance(i,k,c)) and (CNbondDist[1] > distance(i,k,c))):
#                            for m in range(0,natoms):
#                                if ( (atomType[m]==Htype) and (OHbondDist[0] < distance(j,m,c)) and (OHbondDist[1] > distance(j,m,c))):
#                                    # +1 means address converted to id
#                                    #[[C,O,N,H],...]
#                                    # and (NHbondDist[0] < distance(k,m,c)) and (NHbondDist[1] > distance(k,m,c)) 
#                                    if ([i+1,j+1,k+1,m+1] not in bondID ):
#                                        bondID.append([i+1,j+1,k+1,m+1])
#                                        bondFile.write("\ntimestep:" + str(currentStep) + "\n")
#                                        bondFile.write("bondID array: " + str(bondID) + "\n")
#                                        bondFile.write("\nATOMS FOUND: timestep: " + str(currentStep) + " \nC\n ID: " + str(i+1) + " X: " + str(c[i*3]) + " Y: " + str(c[i*3+1]) + " Z: " + str(c[i*3+2]))
#                                        bondFile.write("\nO\n ID: " + str(j+1) + " X: " + str(c[j*3]) + " Y: " + str(c[j*3+1]) + " Z: " + str(c[j*3+2]))
#                                        bondFile.write("\nN\n ID: " + str(k+1) + " X: " + str(c[k*3]) + " Y: " + str(c[k*3+1]) + " Z: " + str(c[k*3+2]))
#                                        bondFile.write("\nH\n ID: " + str(m+1) + " X: " + str(c[m*3]) + " Y: " + str(c[m*3+1]) + " Z: " + str(c[m*3+2]))
#                                    
def excludeN(atomType,currStep,my_rank):
    #this is used to anaylize the bonds file and exclude N's from search if they have 2 active Carbon bonds. 
    #will append as ID
    
    #START FROM THE BEGINNING :9=
    f = open("bonds.txt", "r")
    lineFile = f.readlines()
    for line in lineFile:       #loops through each line of the file
        wordList = line.split()
        if (len(wordList) > 2):
            if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
                thisStep = int(wordList[2])
            elif (currStep == thisStep):
                stored = False                  
                if (int(wordList[0]) in Nlist):                    
                    if (int(wordList[0]) in exclude):        #and we haven't already stored this O
                        stored = True
                    if (not stored):
                        bondnum = int(wordList[2])
                        Ccount = 0
                        tempC = []
                        for i in range(bondnum): 
                            if (int(wordList[3 + i]) in Clist):
                                Ccount+=1
                                tempC.append(int(wordList[3+i]))
                                if (Ccount>1):
                                    if (tempC.size() > 1):
                                        for C in tempC:
                                            excludeC.append(C)
                                    exclude.append(int(wordList[0])) 
    f.close()
    
def getCandN(atomType):
    #get the ID's of carbons and Nitrogens, for use in excludeN
    #stores as atom ID.
    #call at step zero, should never change.
    for i in range(len(atomType)):
        if (atomType[i] == Ctype):
            Clist.append(i+1)
        elif (atomType[i] == Ntype):
            Nlist.append(i+1)
        elif (atomType[i] == Htype):
            Hlist.append(i+1)

def validGroupNH(group):
    #this should see if the N and H are within the NHlist. 
    #called within findOptimal
    #returns true if it has found a valid group. 
    return group[2:] in NHlist

def getNH():
    # THis function gathers NH pairs. N should appear twice, with 2 different H. Data used in checkGroupNH
    #call once to pull data as ID
    difFile = open("difFile.txt",'a')
    f = open("bonds.txt", "r")
    currStep = -1
    lineFile = f.readlines()
    for line in lineFile:       #loops through each line of the file
        wordList = line.split()
        if (len(wordList) > 2):
            if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
                currStep = int(wordList[2])
            if (currStep == 0 and wordList[0] != '#'):
                add = True 
                if (int(wordList[0]) in Nlist and add):
                    bondnum = int(wordList[2])       #gets number of bond this atom has.
                    for i in range(bondnum):
                        if (int(wordList[3 + i]) in Hlist):
                            for NH in NHlist:
                                if (int(wordList[3+i]) in NH):
                                    add = False

                            if (add):
                                NHlist.append([int(wordList[0]),int(wordList[3 + i])])
                                #difFile.write("adding " + wordList[0] + " " + wordList[3+i]+ "\n")
    difFile.write("Nlist: " + str(Nlist) + " within getNH\n")
    difFile.write("NHlist: " + str(NHlist) + " within getNH\n")
    f.close() 
    difFile.close()
    
                                    
main()