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
Htype = 4;

# min and max pair distances for restraint criteria
OHdist = [1.5, 8.0]
NHdist = [0.9, 1.2]
NCdist = [3.0 ,8.0]
COdist = [1.3 ,1.6]

# id lists of atoms which participate in restraint.  Initially empty but
# will be populated by the search function.
# The ith element of each list is associated with the ith element of each other list
# example: Clist[3], Olist[3], Nlist[3], and Hlist[3] all coorespond to one restraint group
#Clist = []
#Olist = []
#Nlist = []
#Hlist = []
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

F1OC = 75
F2OC = 1.0
F1OH = 200
F2OH = 0.75
F1CN = 300
F2CN = 0.75

# equilibrium distances
# R12OC = 1.95
# R12OH = 1.1
# R12CN = 1.3

# equilibrium distances
R12OC = 1.95
R12OH = 1.05
R12CN = 1.4

#time step
timestep = 10000





#file which will contain every single instance where we apply restraint energy
#type is A, which means that rather than writing over the file, which is the case in rest-data.txt, it will log every instance of restraint case being met.
#TO-DO: find a way to include the exact time step in each instance of restraint.
newfile = open("rest-ALLdata.txt",'a')
coordFile = open("coord.txt",'a')



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
    lmp1 = lammps(cmdargs=lammpsArgs) # lammps(comm=MPI.COMM_WORLD,cmdargs=lammpsArgs)
    lmp1.file(lammpsScript)
    initialize()
    lmp1.command("run 0")

    # Get lammps data as python variables
    natoms = lmp1.get_natoms()
    coordinates = lmp1.gather_atoms("x",1,3)
    atomType = lmp1.gather_atoms("type",0,1)
    for i in range(100):
        #gets global quantity ntimestep (current timestep) in lammps.  more extractable stuff can be viewed in library.cpp
        currentStep = lmp1.extract_global("ntimestep",0)
        natoms = lmp1.get_natoms()
        coordinates = lmp1.gather_atoms("x",1,3)
        atomType = lmp1.gather_atoms("type",0,1)
        coordFile.write("Time Step: " + str(currentStep) + "\n")
        #for o in range(natoms):
            #coordFile.write("Atom Index: "+ str(o) + " X: "+ str(coordinates[3*o]) + " Y: "+ str(coordinates[3*o+1]) + " Z: " + str(coordinates[3*o+2]))
            #coordFile.write("\n")
        search(natoms, atomType, coordinates,currentStep)
        lmp1.command("run " + str(timestep)) # lmp1.command("run 100000")
 
    newfile.close()
    coordFile.close()
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
    # End of python script
    



# supporting functions below
#--------------------------#

def search(natoms, atomType, c,currentStep): # c = coordinates
    # this function is used to search for atom pairs and to produce a
    # data file "rest-data.txt"
    # in nested for loops: i = Carbon, j = Oxygen, k = Nitrogen, m = Hydrogen; don't confuse address with id; address = id - 1
    restID = []
    for i in range(0,natoms):
        if atomType[i] == Ctype:
            for j in range(0,natoms):
                if ((atomType[j] == Otype) and (COdist[0] < distance(i,j,c)) and (COdist[1] > distance(i,j,c))):
                    for k in range(0,natoms):
                        if ( (atomType[k] == Ntype) and (NCdist[0] < distance(i,k,c)) and (NCdist[1] > distance(i,k,c))):
                            for m in range(0,natoms):
                                if ( (atomType[m]==Htype) and (OHdist[0] < distance(j,m,c)) and (OHdist[1] > distance(j,m,c)) and (NHdist[0] < distance(k,m,c)) and (NHdist[1] > distance(k,m,c))):
                                    # +1 means address converted to id
                                    #[[C,O,N,H],...]
                                    restID.append([i+1,j+1,k+1,m+1])
                                    coordFile.write("restID array: " + str(restID) + "\n")
                                    newfile.write("\nATOMS FOUND: timestep: " + str(currentStep) + " \nC\n ID: " + str(i+1) + " X: " + str(c[i*3]) + " Y: " + str(c[i*3+1]) + " Z: " + str(c[i*3+2]))
                                    newfile.write("\nO\n ID: " + str(j+1) + " X: " + str(c[j*3]) + " Y: " + str(c[j*3+1]) + " Z: " + str(c[j*3+2]))
                                    newfile.write("\nN\n ID: " + str(k+1) + " X: " + str(c[k*3]) + " Y: " + str(c[k*3+1]) + " Z: " + str(c[k*3+2]))
                                    newfile.write("\nH\n ID: " + str(m+1) + " X: " + str(c[m*3]) + " Y: " + str(c[m*3+1]) + " Z: " + str(c[m*3+2]))
                                    
    # end of nested loops
    # this block checks to remove competing groups, and allows closest groups to pass. if C or N is same ID between 2 groups, they are competing
    #no need to worry about index out of bounds when deleting, seems for loops are modified when you use the del command.
    for i in range(0,len(restID)):
        for j in range(0,len(restID)):
            if (restID[i][0] == restID[j][0] and (i != j)):       #if groups share C
                if (getPerim(restID[j],c) <= getPerim(restID[i]),c):     #if the Perimeter of the group at j is less than Perimeter distance of the group at i, delete group at i. else delete group at j
                    coordFile.write("deleting group: " + str(restID[i]) + "\n")
                    coordFile.write(str(restID[i]) + "Perim = " + str(getPerim(restID[i])) + "\n")
                    coordFile.write(str(restID[j]) + "Perim = " + str(getPerim(restID[j])) + "\n")
                    del restID[i]       #deletes array at i
                else:
                    coordFile.write("deleting group: " + str(restID[j]) + "\n")
                    coordFile.write(str(restID[i]) + "Perim = " + str(getPerim(restID[i])) + "\n")
                    coordFile.write(str(restID[j]) + "Perim = " + str(getPerim(restID[j])) + "\n")
                    del restID[j]
            if (restID[i][2] == restID[j][2] and (i != j)):       #if groups share N
                 if (getPerim(restID[j],c) <= getPerim(restID[i]),c):     #if the Perimeter of the group at j is less than Perimeter distance of the group at i, delete group at i. else delete group at j
                     coordFile.write("deleting group: " + str(restID[i]) + "\n")
                     coordFile.write(str(restID[i]) + "Perim = " + str(getPerim(restID[i])) + "\n")
                     coordFile.write(str(restID[j]) + "Perim = " + str(getPerim(restID[j])) + "\n")
                     del restID[i]       #deletes array at i
                 else:
                     coordFile.write("deleting group: " + str(restID[j]) + "\n")
                     coordFile.write(str(restID[i]) + "Perim = " + str(getPerim(restID[i])) + "\n")
                     coordFile.write(str(restID[j]) + "Perim = " + str(getPerim(restID[j])) + "\n")
                     del restID[j]
    
    coordFile.write("restID array (post removal): " + str(restID) + "\n")

    restfile = open("rest-data.txt",'w')
    restfile.write(str(len(restID)))        #number of groups we would apply force to TODO: must check for distances.
    newfile.write("\n" + "timestep: " + str(currentStep))
    newfile.write(str(len(restID)))
    for i in range(0,len(restID)):
        #O,C
        #O,H
        #C,N
        restfile.write("\n" + str(restID[i][1]) + " " + str(restID[i][0]) + " " + str(R12OC) + " " + str(F1OC) + " " + str(F2OC) + " " + str(COdist[0]) + " " + str(COdist[1]))
        restfile.write("\n" + str(restID[i][1]) + " " + str(restID[i][3]) + " " + str(R12OH) + " " + str(F1OH) + " " + str(F2OH) + " " + str(OHdist[0]) + " " + str(OHdist[1]))
        restfile.write("\n" + str(restID[i][0]) + " " + str(restID[i][2]) + " " + str(R12CN) + " " + str(F1CN) + " " + str(F2CN) + " " + str(NCdist[0]) + " " + str(NCdist[1]))
        newfile.write("\n" + str(restID[i][1]) + " " + str(restID[i][0]) + " " + str(R12OC) + " " + str(F1OC) + " " + str(F2OC) + " " + str(COdist[0]) + " " + str(COdist[1]))
        newfile.write("\n" + str(restID[i][1]) + " " + str(restID[i][3]) + " " + str(R12OH) + " " + str(F1OH) + " " + str(F2OH) + " " + str(OHdist[0]) + " " + str(OHdist[1]))
        newfile.write("\n" + str(restID[i][0]) + " " + str(restID[i][2]) + " " + str(R12CN) + " " + str(F1CN) + " " + str(F2CN) + " " + str(NCdist[0]) + " " + str(NCdist[1]))

        
    restfile.close()
    
    newfile.write("\n" + "timestep: " + str(currentStep))
        
    return 0


def distance(address1, address2,coordinates):
    # this function is used to calculate the distance between 2 atoms.  Input is the address (not id!) of the two atoms
    dx = coordinates[3*address1] - coordinates[3*address2]
    dy = coordinates[3*address1+1] - coordinates[3*address2+1]
    dz = coordinates[3*address1+2] - coordinates[3*address2+2]
    dr = (dx*dx + dy*dy + dz*dz)**(0.5)
    return dr

def initialize():
    restfile = open("rest-data.txt",'w')
    restfile.write("0")
    restfile.close()
    return 0

#return the perimeter of the atoms. only input 1D array, atom group.
def getPerim(atomGroup,coord):
    return distance(atomGroup[0],atomGroup[1],coord) + distance(atomGroup[0],atomGroup[3],coord) + distance(atomGroup[1],atomGroup[2],coord) + distance(atomGroup[2],atomGroup[3],coord)


main()