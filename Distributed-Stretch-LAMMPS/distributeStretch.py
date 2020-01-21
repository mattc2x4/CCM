from lammps import lammps
from mpi4py import MPI

def main():
    mobileAtoms = [] #contains the ID's of the atoms in the mobile region. 
    my_rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()
    if my_rank == 0:
        print ("begin python script")

    # user input; global variables
    lammpsScript = "in.tension.txt"
    lammpsArgs = ["-echo","log"]

    # run lammps script
    lmp1 = lammps(comm=MPI.COMM_WORLD,cmdargs=lammpsArgs)  #lammps(cmdargs=lammpsArgs) 
    lmp1.file(lammpsScript)
    lmp1.command("run 0")
    natoms = lmp1.get_natoms()
    coordinates = lmp1.gather_atoms("x",1,3)
    sim = getSimData(coordinates,natoms,mobileAtoms)
    currL = sim[0]
    v = sim[1]
    dt = sim[2]
    numSteps = 10
    currStep = 0
    debug = open("debug.txt", "a")
    if(my_rank == 0):
        for i in range(numSteps):
            debug = open("debug.txt", "a")
            debug.write("in for loop lmao")
            debug.close()
            lmp1.command("fix 102 top move linear 0.0 0.0 0.01 units box")
            coordinates = lmp1.gather_atoms("x",1,3)
            stretchMobile(coordinates, mobileAtoms,dt,v,currL)
            currL = updateL(currL, v,dt)
            lmp1.scatter_atoms("x",1,3,coordinates)
            lmp1.command("fix 60 mobile nvt temp 300.0 300.0 200")
            lmp1.command("run 1")
        lmp1.command("unfix 60")
        lmp1.close()
    debug.close()
    if (my_rank == 0):
        print ("End of run")
        MPI.Finalize()



def stretchMobile(coordinates, mobileAtoms,dt,v,L):
    #this function should modify the coord list, and then be dispersed in the main function
    debug = open("debug.txt", "a")
    debug.write("stretching")
    for ID in mobileAtoms:      #go throught the list of atoms in the mobile region
        i = ID -1       #find the index of the atom in the coordinates array
        coordinates[3*i+2] = coordinates[3*i+2] + (dt*v)/L      #access the z value for the ith atom (x[3*i+2]), and modify based on formula from sanjib
    debug.close()

def updateL(pastL,v,dt):
    #returns the new length value. should be pretty tiny. 
    #to be called every timestep 
    debug = open("debug.txt", "a")
    debug.write("updating L")
    debug.close()
    return pastL + v * dt

def getSimData(x,numAtoms,mobileAtoms):
    #this should create an array/dictionary which tells us which atom ID's are in the mobile region and which we are interested in modifying. 
    #returns initial L,v,dt.
    debug = open("debug.txt", "a")
    debug.write("getting SIM Data")
    inFile = open("in.tension.txt", "r")
    inLines = inFile.readlines()
    mobileBottom = 0
    mobileTop = 0
    v = 0
    numSteps = 0
    bot = False
    top = False
    dtFound = False
    for line in inLines:
        wordList = line.split()
        if(wordList[0] == "region"):
            if (wordList[1] == "bottom"):
                mobileBottom = float(wordList[8])
                bot = True
                debug.write("bottom gotten: " + str(mobileBottom))
            if (wordList[1] == "top"):
                mobileTop = float(wordList[7])
                debug.write("top gotten: " + str(mobileTop))
                top = True
        if (wordList[0] == "timestep"):
            dt = float(wordList[1])
            dtFound = True
            debug.write("dt gotten: " + str(dt))
        if(bot and top and dtFound):
            break
    if(not(bot and top and dtFound)):
        raise ValueError("bottom, top, rate or timestep or number of steps not found in file.\nExpected bottom syntax: region bottom block INF INF INF INF INF 11.5 units box\nexpected top syntax: region top block INF INF INF INF 35.0 INF units box\nExpected rate syntax: fix 102 top move linear 0.0 0.0 0.01 units box\nexpected timestep syntax: timestep 0.25")
    inFile.close()
    for i in range(numAtoms):
        if(mobileBottom <= x[3*i+2] <= mobileTop):
            mobileAtoms.append(i+1)
    debug.write("mobile Atoms: " + str(mobileAtoms))
    debug.close()
    inFile.close()
    return mobileTop - mobileBottom,v,dt, numSteps

main()