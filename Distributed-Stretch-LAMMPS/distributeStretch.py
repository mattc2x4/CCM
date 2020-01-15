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
    natoms = lmp1.get_natoms()
    coordinates = lmp1.gather_atoms("x",1,3)
    atomType = lmp1.gather_atoms("type",0,1)
    if(my_rank == 0):
        stretchAtoms()
        updateL(currL)
        lmp1.command("run " + str(timestep))

def stretchAtoms():
    #this function should modify the list, and then be dispersed in the main function
    return 0

def updateL(pastL):
    #returns the new length value. should be pretty tiny. 
    return 0

def getRegion(inputFile):
    #this should create an array/dictionary which tells us which atom ID's are in the mobile region and which we are interested in modifying. 
    #might be able to do stuff on a group basis, where sanjib already identified the mobile region. 
    return 0

def storeMobileIDs():
    return 0

def getSimVals(inputFile):
    #get the initial L, dt, v, ztop, zbottom
    return 0

main()