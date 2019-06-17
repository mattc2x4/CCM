# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 09:59:40 2019

@author: mattcohe
"""

# run this python script as master for simulations requiring restraint energy.
# Use the search() function at regular intervals
# to update rest-data.txt (restrained atoms)
# from mpi4py import MPI
# GLOBAL VARIABLES

# atom types
Ctype = 1;
Otype = 2;
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
Clist = []
Olist = []
Nlist = []
Hlist = []

# energy equation parameters (F1, F2)
F1OC = 50
F2OC = 0.5
F1OH = 250
F2OH = 0.75
F1CN = 300
F2CN = 0.75

# equilibrium distances
R12OC = 1.95
R12OH = 1.1
R12CN = 1.5

#time step
timestep = 10000
#current time step
#this has to be a list to get around variable scope issues present in python2.7
timeCurr = [0]



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
    for i in range(2):
        natoms = lmp1.get_natoms()
        coordinates = lmp1.gather_atoms("x",1,3)
        atomType = lmp1.gather_atoms("type",0,1)
        coordFile.write("Time Step: " + str(timeCurr) + "\n")
        for j in range(natoms):
            coordFile.write("Atom Index: "+ str(j) + " X: "+ str(coordinates[j]) + " Y: "+ str(coordinates[j+1]) + " Z: " + str(coordinates[j+2]))
            coordFile.write("\n")
        search(natoms, atomType, coordinates)
        lmp1.command("run " + str(timestep)) # lmp1.command("run 100000")
        timeCurr[0] = timeCurr[0] + timestep

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
    newfile.close()



# supporting functions below
#--------------------------#

def search(natoms, atomType, c): # c = coordinates
    # this function is used to search for atom pairs and to produce a
    # data file "rest-data.txt"
    # in nested for loops: i = Carbon, j = Oxygen, k = Nitrogen, m = Hydrogen; don't confuse address with id; address = id - 1
    for i in range(0,natoms):
        nf = True # nf = 'not found'; we want to break the loops as soon as the carbon (i) is matched with a restraint set of atoms
        if atomType[i] == Ctype:
            for j in range(0,natoms):
                if ((atomType[j] == Otype) and (COdist[0] < distance(i,j,c)) and (COdist[1] > distance(i,j,c)) and ((j+1) not in Olist) and nf ):
                    for k in range(0,natoms):
                        if ( (atomType[k] == Ntype) and (NCdist[0] < distance(i,k,c)) and (NCdist[1] > distance(i,k,c)) and ((k+1) not in Nlist) and nf ):
                            for m in range(0,natoms):
                                if ( (atomType[m]==Htype) and (OHdist[0] < distance(j,m,c)) and (OHdist[1] > distance(j,m,c)) and (NHdist[0] < distance(k,m,c)) and (NHdist[1] > distance(k,m,c)) and ((m+1) not in Hlist) ):
                                    Clist.append(i+1) # +1 means address converted to id
                                    Olist.append(j+1)
                                    Nlist.append(k+1)
                                    Hlist.append(m+1)
                                    nf = False
                                    break
    # end of nested loops
    restfile = open("rest-data.txt",'w')
    restfile.write(str(len(Clist)*3))
    for i in range(0,len(Clist)):
        restfile.write("\n" + str(Olist[i]) + " " + str(Clist[i]) + " " + str(R12OC) + " " + str(F1OC) + " " + str(F2OC) + " " + str(COdist[0]) + " " + str(COdist[1]))
        restfile.write("\n" + str(Olist[i]) + " " + str(Hlist[i]) + " " + str(R12OH) + " " + str(F1OH) + " " + str(F2OH) + " " + str(OHdist[0]) + " " + str(OHdist[1]))
        restfile.write("\n" + str(Clist[i]) + " " + str(Nlist[i]) + " " + str(R12CN) + " " + str(F1CN) + " " + str(F2CN) + " " + str(NCdist[0]) + " " + str(NCdist[1]))
    restfile.close()
    newfile.write(str(len(Clist)*3))
    newfile.write("\n" + "timestep: " + str(timeCurr[0]))
    for i in range(0,len(Clist)):
        newfile.write("\n" + str(Olist[i]) + " " + str(Clist[i]) + " " + str(R12OC) + " " + str(F1OC) + " " + str(F2OC) + " " + str(COdist[0]) + " " + str(COdist[1]))
        newfile.write("\n" + str(Olist[i]) + " " + str(Hlist[i]) + " " + str(R12OH) + " " + str(F1OH) + " " + str(F2OH) + " " + str(OHdist[0]) + " " + str(OHdist[1]))
        newfile.write("\n" + str(Clist[i]) + " " + str(Nlist[i]) + " " + str(R12CN) + " " + str(F1CN) + " " + str(F2CN) + " " + str(NCdist[0]) + " " + str(NCdist[1]))
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

main()