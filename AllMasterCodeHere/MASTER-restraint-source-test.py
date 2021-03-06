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
nonH = []
nonC = []       #non: nonactive of __ type
nonCH = []
nonCO = []
atom_id_dict = {} #dictionary to store atom and molecule id 
#should only be accessed by things that require input as ID, or subtract one. 
#new array: 2d.  restID[0][i] will be associated with a group that is valid to recieve restraint force. written in order [C,O,N,H]
#restOD[i][0] = C, restOD[i][1] = O, restOD[i][2] = N, restOD[i][3] = H
restID = []

# energy equation parameters (F1, F2)
# F1OC = 50
# F2OC = 0.5
# F1OH = 50
# F2OH = 0.75
# F1CN = 300
# F2CN = 0.75

boxdim = [0] *3

F1OC = 75
F2OC = 1.0
F1OH = 50
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

#true if you want all C's and H's to have a force applied. 
CHtag = True

#file which will contain every single instance where we apply restraint energy
#type is A, which means that rather than writing over the file, which is the case in rest-data.txt, it will log every instance of restraint case being met.
coordFile = open("coord.txt",'a')



def main():
    
    # import libraries
    from mpi4py import MPI
    my_rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()
    from lammps import lammps

    if my_rank == 0:
        print ("begin python script")

    # user input; global variables
    lammpsScript = "in.accelerated-test.txt"
    lammpsArgs = ["-echo","log"]

    #--------------END OF USER DEFINED PARAMETERS--------------------------------------

    # run lammps script
    lmp1 = lammps(comm=MPI.COMM_WORLD,cmdargs=lammpsArgs)  #lammps(cmdargs=lammpsArgs) 
    lmp1.file(lammpsScript)
    initialize()
    lmp1.command("run 0")
    natoms = lmp1.get_natoms()
    coordinates = lmp1.gather_atoms("x",1,3)
    atomType = lmp1.gather_atoms("type",0,1)
    if(my_rank == 0):
        getTypes(atomType)
        getNH()
        if(CHtag):
            getCH()
        getCO()      #CHANGE
        writeatomid(atom_id_dict)
    for i in range(100):
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
            search(natoms, atomType, coordinates,currentStep,atom_id_dict)
        MPI.COMM_WORLD.Barrier()
        lmp1.command("run " + str(timestep)) # lmp1.command("run 100000")
    lmp1.close()
    if my_rank == 0:
        print ("End of run")
        MPI.Finalize()
    # End of python script
    



# supporting functions below
#--------------------------#

def search(natoms, atomType, c,currentStep,atom_id_dict): # c = coordinates
    # this function is used to search for atom pairs and to produce a
    # data file "rest-data.txt"
    # in nested for loops: i = Carbon, j = Oxygen, k = Nitrogen, m = Hydrogen; don't confuse address with id; address = id - 1
    #makes calls to validGroup, to check N and H items.
    #makes call to findOptimal, to remove competing groups. 
    newfile = open("rest-ALLdata.txt", 'a')
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
                                    #newfile.write("\nATOMS FOUND: timestep: " + str(currentStep) + " \nC\n ID: " + str(i+1) + " X: " + str(c[i*3]) + " Y: " + str(c[i*3+1]) + " Z: " + str(c[i*3+2]))
                                    #newfile.write("\nO\n ID: " + str(j+1) + " X: " + str(c[j*3]) + " Y: " + str(c[j*3+1]) + " Z: " + str(c[j*3+2]))
                                    #newfile.write("\nN\n ID: " + str(k+1) + " X: " + str(c[k*3]) + " Y: " + str(c[k*3+1]) + " Z: " + str(c[k*3+2]))
                                    #newfile.write("\nH\n ID: " + str(m+1) + " X: " + str(c[m*3]) + " Y: " + str(c[m*3+1]) + " Z: " + str(c[m*3+2]))
    exclude_C(atom_id_dict,c,restID) #CHANGE
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
    
    if (CHtag):
        restfile.write(str(len(restID) * 3 + len(nonCH)))
    else:
        restfile.write(str(len(restID) * 3))        #number of groups we would apply force to 
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
        newfile.write( "\n" + str(restID) + "TimeStep: " + str(currentStep) + "\n")
        newfile.write("\n" + str(restID[i][1]) + " " + str(restID[i][0]) + " " + str(R12OC) + " " + str(F1OC) + " " + str(F2OC) + " " + str(COdist[0]) + " " + str(COdist[1]))
        newfile.write("\n" + str(restID[i][1]) + " " + str(restID[i][3]) + " " + str(R12OH) + " " + str(F1OH) + " " + str(F2OH) + " " + str(OHdist[0]) + " " + str(OHdist[1]))
        newfile.write("\n" + str(restID[i][0]) + " " + str(restID[i][2]) + " " + str(R12CN) + " " + str(F1CN) + " " + str(F2CN) + " " + str(NCdist[0]) + " " + str(NCdist[1]))
    if(CHtag):
        writeCH(newfile,restfile)
    difFile = open("difFile.txt",'a')
    writeCO(difFile,restfile)     
    restfile.close()
    difFile.close()
    newfile.close()
    coordFile.close()
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
    
def getTypes(atomType):
    #get the ID's of carbons and Nitrogens, for use in excludeN
    #gets nonactive C and H for restforce app.
    #stores as atom ID.
    #call at step zero, should never change.
    for i in range(len(atomType)):
        if (atomType[i] == Ctype):
            Clist.append(i+1)
        elif (atomType[i] == Ntype):
            Nlist.append(i+1)
        elif (atomType[i] == Htype):
            Hlist.append(i+1)
        elif(atomType[i] == 4):     #non active H
            nonH.append(i+1)
        elif(atomType[i] == 1):     #non active C
            nonC.append(i+1)
        elif (atomType[i] == Otype):
            Olist.append(i+1)

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
                                if (int(wordList[3+i]) in NH):       #each H should only appear once
                                    add = False

                            if (add):
                                NHlist.append([int(wordList[0]),int(wordList[3 + i])])
                                #difFile.write("adding " + wordList[0] + " " + wordList[3+i]+ "\n")
    difFile.write("Nlist: " + str(Nlist) + " within getNH\n")
    difFile.write("NHlist: " + str(NHlist) + " within getNH\n")
    f.close() 
    difFile.close()
    
def getCH():
    # THis function gathers CH pairs. C can appear 3 times (?) and should have every H bonded to it listed
    #this will be added to nonCH
    #call once to pull data as ID
    #part of applying forces to all C-H pairs.
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
                if ((int(wordList[0]) in Clist or int(wordList[0]) in nonC) and add):       #if carbon is active or non active, TODO check this. should we do active C? i think yes. 
                    bondnum = int(wordList[2])       #gets number of bond this atom has.
                    for i in range(bondnum):
                        if (int(wordList[3 + i]) in nonH):      
                            for CH in nonCH:
                                if (int(wordList[3+i]) in CH):      #each H should only appear once
                                    add = False

                            if (add):
                                nonCH.append([int(wordList[0]),int(wordList[3 + i])])
                                #difFile.write("adding " + wordList[0] + " " + wordList[3+i]+ "\n")
    difFile.write("CHlist: " + str(nonCH) + " within getNH\n")
    f.close() 
    difFile.close()
    
def writeCH(difFile,restfile):
    #will write all CH pairs to the rest-data file. meant to be toggleable. 
    #difFile = open("rest-ALLdata.txt",'a')
    #restfile = open("rest-data.txt",'w')
    for CH in nonCH:
        difFile.write("\n" + str(CH[0]) + " " + str(CH[1]) + " " + str(1.1) + " " + str(200) + " " + str(.75) + " " + str(COdist[0]) + " " + str(COdist[1]))
        restfile.write("\n" + str(CH[0]) + " " + str(CH[1]) + " " + str(1.1) + " " + str(200) + " " + str(.75) + " " + str(COdist[0]) + " " + str(COdist[1]))
    #difFile.close()
    #restfile.close()




def getCO():
    #pls remove most of the write statements which were only there for checking purposes
    difFile=open("difFile.txt",'a')
    f=open("bonds.txt",'r')
    lineFile=f.readlines()
    for line in lineFile:
        wordList=line.split()
        if (len(wordList) > 2):
            if ((wordList[0] == '#') and (wordList[1] == 'Timestep')):
                currStep=int(wordList[2])
            if ((currStep == 0) and (wordList[0] != '#')):
                add = True
                if ((int(wordList[0]) in Clist) and  add):
                    bondnum=int(wordList[2])
                    (Oxycheck,NonCcheck)=(False,False)
                    for i in range(bondnum):
                        if int(wordList[3+i]) in Olist:
                            Oxycheck=True
                            Oxyindex=3+i 
                        if int(wordList[3+i]) in nonC:
                            NonCcheck=True
                            NonCindex=3+i
                    if Oxycheck and NonCcheck:
                        nonCO.append([int(wordList[NonCindex]),int(wordList[Oxyindex])])  
    difFile.write("COlist: " + str(nonCO) + " within getNH\n")
    f.close() 
    difFile.close()
                        
                                   
def writeCO(difFile,restfile):
    tempnew=open('CO_debug_testing.txt','w')
    for CO in nonCO:
        difFile.write("\n" + str(CO[0]) + " " + str(CO[1]) + " " + str(1.1) + " " + str(200) + " " + str(.75) + " " + str(COdist[0]) + " " + str(COdist[1]))
        tempnew.write("\n" + str(CO[0]) + " " + str(CO[1]) + " " + str(1.1) + " " + str(200) + " " + str(.75) + " " + str(COdist[0]) + " " + str(COdist[1]))
        restfile.write("\n" + str(CO[0]) + " " + str(CO[1]) + " " + str(1.1) + " " + str(200) + " " + str(.75) + " " + str(COdist[0]) + " " + str(COdist[1]))
    tempnew.close()

def writeatomid(atom_id_dict):
    #Writes the atom ids inot the dictionary to be used in the Exclude C function
    
    f=open("8Epon-4DETDA-H.txt","r")
    g=open("Temp_Atom_id.txt","w")
    currStep=-1
    lineFile = f.readlines() #makes a list with each element of the list pertaining to a line in the file
    for line in lineFile:       #loops through each line of the file
        wordList=line.split()
        if (len(wordList) > 4): 
            atom_id_dict[wordList[0]]=wordList[1] #assigning atomid to atom name 
    for i in atom_id_dict:
        g.write("Atom ID:"+" "+str(i)+" "+"Molecule Id:"+" "+str(atom_id_dict[i])+'\n') #just for testing purporses
    g.write(str(atom_id_dict)+"\n") 
    f.close()
    g.close()


def exclude_C(atom_id_dict,coord,restID):
    g=open("Temp_Atom_id.txt","a")
    g.write(str(restID)+"\n")
    delete_index=[] # will contain the indices which will be deleted
    for i in range(len(restID)):
        for j in range(len(restID)):
            if (atom_id_dict[str(restID[i][0])]==atom_id_dict[str(restID[j][0])]) and  (i != j) and (restID[i][0] != restID[j][0]):  #checks whether two C's share the same molecule #id's in form of string
                if getPerim(restID[i],coord) < getPerim(restID[j],coord): #if Perimeter is greater for the second group, delete it
                    if j not in delete_index:
                        delete_index.append(j)
                        g.write("Groups to be deleted:"+str(restID[j])+'\n')
                else:
                    if i not in delete_index:
                        delete_index.append(i)
                        g.write("Groups to be deleted:"+str(restID[i])+'\n')
    g.write(str(delete_index)+'\n')
    g.close()
    count=0
    for i in range(len(delete_index)):
        del restID[delete_index[i]]
        count=count+1
        if i==len(delete_index)-1:
            continue
        else:
            delete_index[i+1]-=count
    
main()