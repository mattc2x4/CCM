#Constants involved in force app
F1OH = 50
F2OH = 0.75
R12OH = 1.05
OHdist = [1.5, 8.0]

#Simulation based data
restID = []     #[[gpsO, silH, silO, gpsH],...]
timestep = 50000
boxdim = [0] *3
gpsOType = 4
gpsHType = 7
gpsSiType = 5
silHType = 3
silOType = 2

#dictionaries which contain special atom ID's
#active denoting that they will recieve force
activeSilO = {}
activeGpsO = {}
activeGpsH = {}
SilHList = {}
gpsSiList = {}
gpsOList = {}


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
    lmp1.command("run 0")
    natoms = lmp1.get_natoms()
    coordinates = lmp1.gather_atoms("x",1,3)
    atomType = lmp1.gather_atoms("type",0,1)
    getTypes(atomType)
    getActiveAtoms()
    for i in range(1):
        currentStep = lmp1.extract_global("ntimestep",0)
        natoms = lmp1.get_natoms()
        boxdim[0] = lmp1.extract_global("boxxhi",1)
        boxdim[1] = lmp1.extract_global("boxyhi",1)
        boxdim[2] = lmp1.extract_global("boxzhi",1)
        coordinates = lmp1.gather_atoms("x",1,3)
        atomType = lmp1.gather_atoms("type",0,1)
        if(my_rank == 0):
            #coordFile.write("Time Step: " + str(currentStep) + "\n")
            search(natoms, atomType, coordinates,currentStep)
        MPI.COMM_WORLD.Barrier()
        lmp1.command("run " + str(timestep))
    lmp1.close()
    if my_rank == 0:
        print ("End of run")
        MPI.Finalize()


def search(natoms, atomType, c,currentStep):
    #goes throught the active type dictionaries and gets the groups which satisfy distance criteria. not 100%sure what to use for that yet though.
    for gpsO in activeGpsO:
        for silO in activeSilO:
            silOgpsHdist = distance(activeGpsO[gpsO]-1, silO-1, c)
            silHgpsOdist = distance(activeSilO[silO] - 1, gpsO-1,c)
            if(OHdist[0] < silHgpsOdist  and silHgpsOdist < OHdist[1] and OHdist[0] < silOgpsHdist and silOgpsHdist < OHdist[1]):
                restID.append([gpsO, activeSilO[silO], silO, activeGpsO[gpsO]])
    newfile = open("rest-ALLdata.txt", 'a')
    newfile.write("restID: " + str(restID))

         


def removeInoptimalGroups():
    #removes groups that contain the same silica or gps OH, based on which are the fatherst away from each other.
    return 0

def deleteWater():
    return 0

def getActiveAtoms():
    #gets active O's in silica, and active OH in gps.
    difFile = open("difFile.txt",'a')
    #difFile.write("entering getActiveAtoms\n")
    f = open("bonds.txt", "r")
    currStep = -1
    lineFile = f.readlines()
    for line in lineFile:       #loops through each line of the file
        wordList = line.split()
        if (len(wordList) > 2):
            if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
                currStep = int(wordList[2])
            if (currStep == 0 and wordList[0] != '#'):
                if (int(wordList[0]) in SilHList):
                    #difFile.write("found SilH\n")
                    bondnum = int(wordList[2])       #gets number of bond this atom has.
                    for i in range(bondnum):
                        activeSilO[int(wordList[3 + i])] = int(wordList[0])     #adds every atom bonded to H in silicone to activeSilO, should only be oxygens.
                        #difFile.write("adding SilO\n")

                elif(int(wordList[0]) in gpsSiList):        #this is for getting gps O's
                    bondnum = int(wordList[2])       #gets number of bond this atom has.
                    for i in range(bondnum):
                        if(int(wordList[3+i]) in gpsOList):
                            activeGpsO[int(wordList[3 + i])] = 1     #adds every atom bonded to H in silicone to activeGpsO, should only be oxygens.
    
    #this section has to have the gps O loaded to work, so it goes after. It gets the activeGpsH, which is bonded to the O.
    #look for activeGpsO and record things they're bonded to, as long as it isnt a Si.
    for line in lineFile:       #loops through each line of the file
        wordList = line.split()
        if (len(wordList) > 2):
            if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
                currStep = int(wordList[2])
            if (currStep == 0 and wordList[0] != '#'):
                if (int(wordList[0]) in activeGpsO):
                    bondnum = int(wordList[2])       #gets number of bond this atom has.
                    for i in range(bondnum):
                        if(int(wordList[3+i]) not in gpsSiList):
                            activeGpsH[int(wordList[3 + i])] = 1
                            activeGpsO[int(wordList[0])] = int(wordList[3 + i])  #makes value of activeGpsO the H bonded with it. 
    difFile.write("activeSilO: " + str(activeSilO) + "\n")
    difFile.write("activeGpsH: " + str(activeGpsH) + "\n")
    difFile.write("activeGpsO: " + str(activeGpsO) + "\n")
    f.close()
    difFile.close()  


def getTypes(atomType):
    #gets H's in silica slab, used to find active O's.
    #gets Si's in gps to find active O's
    #gets all O's in gps to differentiate from the C thats bonded to a silicone.
    #stores as atom ID.
    #call at step zero, should never change.
    difFile = open("difFile.txt",'a')
    for i in range(len(atomType)):
        if (atomType[i] == silHType):
            SilHList[i+1] = 1
        elif (atomType[i] == gpsSiType):
            gpsSiList[i+1] = 1
        elif (atomType[i] == gpsOType):
                gpsOList[i+1] = 1
    difFile.write("SilHList: " + str(SilHList) + "\n")
    difFile.write("gpsOList: " + str(gpsOList) + "\n")
    difFile.write("gpsSiList: " + str(gpsSiList) + "\n")

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

main()
