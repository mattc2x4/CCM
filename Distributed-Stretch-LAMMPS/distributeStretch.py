from lammps import lammps
from mpi4py import MPI

def main():
    mobileAtoms = [] #contains the ID's of the atoms in the mobile region. 
    my_rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()
    if my_rank == 0:
        print ("begin python script")

    # user input; global variables
    lammpsScript = "in.accelerated-test.txt"
    lammpsArgs = ["-echo","log"]

    # run lammps script
    lmp1 = lammps(comm=MPI.COMM_WORLD,cmdargs=lammpsArgs)  #lammps(cmdargs=lammpsArgs) 
    timestep = 0
    lmp1.file(lammpsScript)
    lmp1.command("run 0")
    natoms = lmp1.get_natoms()
    coordinates = lmp1.gather_atoms("x",1,3)
    atomType = lmp1.gather_atoms("type",0,1)
    sim = getSimData(coordinates,natoms,mobileAtoms)
    currL = sim[0]
    v = sim[1]
    dt = sim[2]
    if(my_rank == 0):
        stretchMobile(coordinates)
        currL = updateL(currL, v,dt)
        lmp1.command("run " + str(1))

def stretchMobile(coordinates):
    #this function should modify the list, and then be dispersed in the main function
    return 0

def updateL(pastL,v,dt):
    #returns the new length value. should be pretty tiny. 
    #to be called every timestep 
    return pastL + v * dt

def getSimData(x,numAtoms,mobileAtoms):
    #this should create an array/dictionary which tells us which atom ID's are in the mobile region and which we are interested in modifying. 
    #returns initial L,v,dt.

    inFile = open("in.tension.txt", "r")
    inLines = inFile.readlines()
    mobileBottom = 0
    mobileTop = 0
    v = 0
    bot = False
    top = False
    rate = False
    dtFound = False
    for line in inLines:
        wordList = line.split()
        if(wordList[0] == "region"):
            if (wordList[1] == "bottom"):
                mobileBottom = float(wordList[8])
                bot = True
                print("bottom gotten: " + str(mobileBottom))
            if (wordList[1] == "top"):
                mobileTop = float(wordList[7])
                print("top gotten: " + str(mobileTop))
                top = True
        if (wordList[0] == "fix" and wordList[2] == "top" and wordList[3] == "move"):
            v = float(wordList[7])
            rate = True
            print("v gotten: " + str(v))
        if (wordList[0] == "timestep"):
            dt = float(wordList[1])
            dtFound = True
            print("dt gotten: " + str(dt))
        if(bot and top and rate and dtFound):
            break
    if(not(bot and top and rate and dtFound)):
        raise ValueError("bottom, top, rate or timestep not found in file.\nExpected bottom syntax: region bottom block INF INF INF INF INF 11.5 units box\nexpected top syntax: region top block INF INF INF INF 35.0 INF units box\nExpected rate syntax: fix 102 top move linear 0.0 0.0 0.01 units box\nexpected timestep syntax: timestep 0.25")
    inFile.close()
    for i in range(numAtoms):
        if(mobileBottom <= x[3*i+2] <= mobileTop):
            mobileAtoms.append(i+1)
    return mobileTop - mobileBottom,v,dt

main()