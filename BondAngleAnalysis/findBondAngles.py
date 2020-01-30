# -*- coding: utf-8 -*-
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.stats
import seaborn as sns
import timeit
"""
Created on Tue Oct 29 15:45:57 2019

@author: mattcohe
"""

"""Constants"""
#this is where you imput the 3 points of the angle. 
#atom type 1 = Al 
#atom type 2 = Mg
#atom type 3 = O
#atom type 4 = Si

vertexType = 1   #VERTEX AND END TYPES HERE.
endType1 = 3    #the other two types you want to calculate angle of
endType2 = 3

"""end Constants"""


atomList = []   #this will hold all angle data types. will be in order, with the 0th atom being the atom with ID 1. Per timestep value
#should only contain atoms for the current timestep
angleVals = [] #containts angle values in raw form, to be graphed
angList = []  #this will hold all angle objects. 
#clear angList every timestep.
#end2ID interchangable locations. 

boxdim = [0,0,0]     #[x,y,z] lengths

def main():
    #CONFIGURE FILES
    #INSERT DUMP_FINAL FILES TO DUMP
    #INSERT BOND FILES TO BONDS
    #i = 0
    dump = "dump_final.lammps"      #WRITE FILE NAMES HERE  store file names as strings, to be opened and closed by each subroutine. 
    bonds = "MD_bonds_final.reaxc"
    markedXYZ = "marked_dump.xyz"
    markedlammps = "marked_dump.lammps"
    analyzeSimAndMark(dump,bonds, markedXYZ,markedlammps)
 

def analyzeSimAndMark(dump, bonds, markedXYZ, markedlammps):
    #TODO: add different file name spots and make clear which are which, like for xyz translation and marked file names. 
    simData = getSimData(dump)
    currStep = simData[0]
    finalStep = simData[1]
    #print(currStep)
    #print(finalStep)
    incrSize = simData[2]
    # i = 0
    #print (simData)
    while(currStep <= finalStep):
        print("currStep: " + str(currStep))
        fillAtomList(dump, currStep)
        #print("AtomList[0]: " + str(atomList[0]))
        getAngleID(bonds,currStep)
        calcAngles(currStep)
        #print("average value: " + str(sum(angleVals) / len(angleVals)))
        markAtomsDumpAll(func, dump, markedlammps, currStep)
        currStep += incrSize
        angList.clear()
        atomList.clear()
        # i+=1
        # if(i>1):
        #     break
    plot(angleVals,"sim")
    #print(len(angleVals))
    lammpsToXYZ(markedlammps, markedXYZ, {1:"Al", 2:"Mg", 3:"O", 4:"Si"})


def getSimData(dump):
    #this function returns a list [firstFrame, finalFrame, stepSize]
    #print("getting Sim Data")
    dump = open(dump,"r")
    numStepFound = 0
    step = []
    dumpLine = dump.readlines()
    i = 0
    lastFound = False
    j = len(dumpLine) - 1
    while(numStepFound<2):
        dumpWord = dumpLine[i].split()
        if(dumpWord[0] == "ITEM:" and dumpWord[1] == "TIMESTEP"):
            numStepFound+=1
            step.append(int(dumpLine[i+1].split()[0]))
            #print("timestep found")
        i+=1

    #print("out of first while")
    while(not lastFound):
        dumpWord = dumpLine[j].split()
        if(dumpWord[0] == "ITEM:" and dumpWord[1] == "TIMESTEP"):
            step.append(int(dumpLine[j+1].split()[0]))
            lastFound = True
            #print("last found")
        j-=1
    stepSize = step[1] - step[0]
    #print(step)
    #print(stepSize)
    del step[1]
    step.append(stepSize)
    dump.close()
    return step


def fillAtomList(dump,currStep):
    #read in the atom data from timestep "currStep"
    #finds meantions of timestep in header. compares with input (master) step. if it is greater than, breaks loop, allowing
    #main code to run instead. when it finds the correct 
    #updates Boxdim, which changes the box dimensions. this affects finding distance with periodic boundaries. 
    #print("filling atom list for " + str(currStep) + " step.")
    startRead = False
    readBox = False
    step = False
    flag = True
    dump = open(dump,"r")
    dumpLine = dump.readlines()
    for i in range(len(dumpLine)):
        dumpWord = dumpLine[i].split()
        if(dumpWord[0] == "ITEM:" and dumpWord[1] == "TIMESTEP"):    #read timestep, if its correct continue. otherwise, break loop.
            # print("found TIMESTEP header")
            # print(dumpWord)
            myStep = int(dumpLine[i+1].split()[0])
            if(myStep > currStep):
                # print("breaking at " + str(myStep))
                break
            elif(myStep == currStep):
                step = True
                # print("found correct Timestep: " + str(myStep))
                # print(dumpWord)
            #print("myStep : " + str(myStep) + " currStep: " + str(currStep))
        if(dumpWord[0] == "ITEM:" and dumpWord[1] == "BOX" and step and not readBox):    #UPDATING BOXDIM.  havent read box yet this step,and have read in the timestep.
            # print("found BOX header.")
            # print(dumpWord)
            readBox = True
            currLine = i
            box1 = dumpLine[currLine + 1].split()   #getting three lines which contain X lo hi, y lo hi, z lo hi for boxdim.
            box2 = dumpLine[currLine + 2].split()   #Y
            box3 = dumpLine[currLine + 3].split()   #Z
            boxdim[0] = float(box1[1]) - float(box1[0])
            boxdim[1] = float(box2[1]) - float(box2[0])
            boxdim[2] = float(box3[1]) - float(box3[0])
            #print("updating boxdim: " + str(boxdim))
        if(dumpWord[0] == "ITEM:" and dumpWord[1] == "ATOMS" and readBox):    #when we see this header we want to skip this iteration, then continue on the next line, hence continue. 
            # print("found ATOMS header")
            # print(dumpWord)
            startRead = True
            continue
        if (startRead and step):
            if (flag):
                #print("loading atoms")
                # print(dumpWord)
                flag = False
            atomList.append(Atom(float(dumpWord[3]), float(dumpWord[4]), float(dumpWord[5]), int(dumpWord[0]), int(dumpWord[1])))
            #this adds an atom object in atomList. x,y,z,ID,TYPE
    dump.close()
            

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
    #accepts Atom objects as input
   #Returns the angle in radians between vectors 'v1' and 'v2'::
   a = distance(vertAtom,endAtom2)
   #print(a)
   b = distance(vertAtom, endAtom1)
   #print(b)
   c = distance(endAtom1,endAtom2)
   #print(c)
   return(math.degrees(math.acos((a**2+b**2-c**2)/(2*a*b))))
   
def getAngleID(bonds,currStep):
    # This function will be called on each timestep. It will collect all angles in the timestep, and then add them to angList. this list
    #will be modified to calculate the angle later. 
    #call once to pull data as ID
    #print("getting angles for " + str(currStep) + " step.")
    bonds = open(bonds, "r")
    lineFile = bonds.readlines()
    read = False
    for line in lineFile:       #loops through each line of the file
        wordList = line.split()
        endCount = 0
        if (len(wordList) > 2):
            if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # this checks timestep. if timestep in file matches input, read lines.upon encountering another, break loop.
                #print("\tCurrStep: " + str(currStep) + " Step found: " + wordList[2])
                if(currStep == int(wordList[2])):
                    #print("\tfound timestep")
                    read = True
                elif(currStep > int(wordList[2])):
                    continue
                else:
                    #print("\tbreaking\n")
                    break
            if (read == True and wordList[0] != '#'):
                #print(str(atomList[int(wordList[0]) - 1].TYPE) + str(atomList[int(wordList[0]) - 1].TYPE == vertexType))
                if (atomList[int(wordList[0]) - 1].TYPE == vertexType):
                    #access the corresponding atom data. subtract one to translate to index. Make sure ID is correct.  and int(wordList[0]) == atomList[int(wordList[0]) - 1].ID
                    #print("\tfound vertex")
                    bondnum = int(wordList[2])       #gets number of bond this atom has.
                    for i in range(bondnum):
                        try: 
                            atomList[int(wordList[3 + i]) - 1].TYPE  == endType1 or atomList[int(wordList[3 + i]) - 1].TYPE == endType2
                        except IndexError:
                            print("Files Incomplete: an atom exists in bonds file that does not exist in dump this atom is: " + str(int(wordList[3 + i])))
                        if (atomList[int(wordList[3 + i]) - 1].TYPE  == endType1 or atomList[int(wordList[3 + i]) - 1].TYPE == endType2):      
                           endCount+=1
                           #print("\tfound end")
                           if (endCount == 1):     #this is the first one we found, store in temp
                                firstEndID = int(wordList[3+i])
                           elif(endCount == 2):    #this is the second one we found, add time, reset counter.
                                angList.append(Angle(int(wordList[0]),firstEndID,int(wordList[3+i]),currStep,-1))
                                endCount = 0
    #print("\t" + str(len(angList)))
    bonds.close()

def calcAngles(currStep):
    #this function will take the info in angList, and use it to calculate angle values based on the atom data in atomList.
    #called on angLIst once per timestep, though angList will contain all steps. may be problematic due to data sizes. 
    #TODO: modify this to decrease compexity based on what yoon wants. this dude needs to respond faster.
    #print("calcing Angles for: " + str(currStep) + " step.")
    for ang in angList:
        if (ang.timestep == currStep):
            ang.angle = cosLaw(atomList[ang.vertID-1],atomList[ang.end1ID-1],atomList[ang.end2ID-1])
            angleVals.append(ang.angle)

        
def plot(arr,currStep):
    # sns.distplot(arr, hist = False, kde = True, rug = True,
    #          color = 'darkblue', 
    #          kde_kws={'linewidth': 3},
    #          rug_kws={'color': 'black'})

    sns.distplot(arr, hist=True, kde=True, 
             bins=int(180/5), color = 'darkblue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4})
    # x = ax.lines[0].get_xdata() # Get the x data of the distribution
    # y = ax.lines[0].get_ydata() # Get the y data of the distribution
    # maxid = np.argmax(y) # The id of the peak (maximum of y data)
    # plt.plot(x[maxid],y[maxid], 'bo', ms=10)

    # Plot formatting
    plt.title("Density Plot: "  + currStep)
    plt.xlabel('Angle (deg)')
    plt.ylabel('Density')
    plt.savefig("Plot-" + currStep + ".png")
    plt.show()
           
def plotNorm(arr,currStep):
    x_min = 0.0
    x_max = 180.0

    mean = sum(arr) / len(arr)
    std = 2.0

    y = scipy.stats.norm.pdf(arr,mean,std)

    plt.plot(arr,y, color='coral')

    plt.grid()

    plt.xlim(x_min,x_max)
    plt.ylim(0,0.25)

    plt.title('How to plot a normal distribution in python with matplotlib',fontsize=10)

    plt.xlabel('x')
    plt.ylabel('Normal Distribution')

    plt.savefig("normal_distribution.png")
    plt.show()

def markAtomsDumpAll(func, dump, markedlammps, currStep):
    #this function takes func: which should be a function that evaluates to a boolean, and should take an angle object as an input
    #to be called once per timestep. it should copy lines that don't evauluate to true into the file, and for lines that do evaluate to true it should copy lines, 
    # and add a column with some value. 
    #this should mark all atoms in the angle, for this specific timestep. 
    markDict = {}       #dictionary containing the ID of all atoms in angle that satisfies func condition.
    marks = {1:10,2:20,3:30,4:40}          #if the angle fits the criteria described in func, it will be marked with the value here. Otherwise, it will be marked by its type.  
    markedFile = open(markedlammps, "a")
    dump = open(dump,"r")
    #print("Marking")
    step = False
    atomHeaderSeen = False
    for ang in angList:
        #print("checking: " + str(ang))
        if (func(ang)):
            #print("ADDING")
            markDict[ang.vertID] = 1
            markDict[ang.end1ID] = 1
            markDict[ang.end2ID] = 1
    #print(str(markDict))
    dumpLine = dump.readlines()
    for i in range(len(dumpLine)):
        dumpWord = dumpLine[i].split()
        #print("DUMPWORD =" + dumpWord)
        if(dumpWord[0] == "ITEM:" and dumpWord[1] == "TIMESTEP"):    #read timestep, if its correct continue. otherwise, break loop.
            myStep = int(dumpLine[i+1].split()[0])
            #print("\tfound timestep header")
            if(myStep > currStep):
                #print("\tMARKING breaking at " + str(myStep))
                break
            elif(myStep == currStep):
                #print("\tstep found MARKING")
                step = True
                markedFile.write(dumpLine[i])
                continue
        if(step):
            if(dumpWord[0] == "ITEM:" and dumpWord[1] == "ATOMS"):
                atomHeaderSeen = True
                markedFile.write(dumpLine[i] + " mark1")
                continue
                #print("\tatomHeaderSeen")
            if(atomHeaderSeen):
                line = dumpLine[i].rstrip('\n')
                splitLine = line.split()
                #print(line)
                if(int(dumpWord[0]) in markDict):
                    #print("writing")
                    line += str(marks[int(splitLine[1])])
                    print(line, file=markedFile)
                else:
                    line += str(splitLine[1])
                    print(line, file=markedFile)
            else:
                line = dumpLine[i].rstrip('\n')
                print(line, file=markedFile)
    dump.close()
    markedFile.close()

def lammpsToXYZ(inputFile, outputFile,transDict):
    #takes input of strings for names of input and output files, and a dictionary described below. 
    #trans Dict should be in the following format:
    #transDict = {1:Al,2:O} where 1 is the type that corresponds to Aluminum atoms. etc, fill for all types available. 
    #print("converting lammps to xyz")
    inp = open(inputFile,"r")
    out = open(outputFile,"w")
    inpLine = inp.readlines()
    xlen = 0
    ylen = 0
    zlen = 0
    for i in range(len(inpLine)):
        wordList = inpLine[i].split()
        if(len(wordList) >= 2):
            if (wordList[0] == "ITEM:"):
                if(wordList[1] == "TIMESTEP"):
                    myStep = int(inpLine[i+1].split()[0])
                elif(wordList[1] == "NUMBER" and wordList[3] == "ATOMS"):
                    numAtoms = int(inpLine[i+1].split()[0])
                elif(wordList[1] == "BOX"):
                    box1 = inpLine[i + 1].split()   #getting three lines which contain X lo hi, y lo hi, z lo hi for boxdim.
                    box2 = inpLine[i + 2].split()   #Y
                    box3 = inpLine[i + 3].split()   #Z
                    xlen = round(float(box1[1]) - float(box1[0]),6)
                    ylen = round(float(box2[1]) - float(box2[0]),6)
                    zlen = round(float(box3[1]) - float(box3[0]),6)
                    print(str(numAtoms), file=out)
                    print("Atoms. Timestep: " + str(myStep) + "          " + str(xlen) + "  " + str(ylen) + "  " + str(zlen) + "  90.000000  90.000000  90.000000",file=out) #type cell system angles here
            elif(len(wordList) > 5):
                xyzline = transDict[int(wordList[1])] + "  " + str(wordList[3]) + "  " + str(wordList[4]) + "  " + str(wordList[5]) + "  "
                for mark in range(6,len(wordList)):
                    xyzline += wordList[mark] + "  "
                print(xyzline, file=out)
    inp.close()
    out.close()
    #print("\tDone.")






    
    

def func(ang):
    #this is the function that users should modify. it takes in one angle object, and makes some comparison with angle fields, to return a boolean.
    #the example code included will mark any atom that is below 60 degrees.
    if(ang.angle < 90):
        return True
    else:
        return False

class Atom:
    #this is an atom data structure. This is used to store each atom's x,y,z vals, Id, and Type. default values are -1.
    def __init__(self, x,y,z,ID,TYPE):
        self.x = x
        self.y = y
        self.z = z
        self.ID = ID
        self.TYPE = TYPE
    def __str__(self):      #printing atom object will yield the following
        return ("Atom ID: " + str(self.ID) + "\n" + "TYPE: " + str(self.TYPE) + "\n" + "X: " + str(self.x) + "\n" + "Y: " + str(self.y) + "\n" + "Z: " + str(self.z))

class Angle:
    #this is an angle data structure. This is used to store each angle's data as seen below.
    def __init__(self, vertID,end1ID,end2ID,timestep,angle):
        self.vertID = vertID
        self.end1ID = end1ID
        self.end2ID = end2ID
        self.timestep = timestep
        self.angle = angle
    def __str__(self):      #printing angle object will yield the following
        return ("Vertex ID: " + str(self.vertID) + "\n" + "End1ID: " + str(self.end1ID) + "\n" + "End2ID: " + str(self.end2ID) + "\n" + "Timestep: " + str(self.timestep) + "\n" + "angle: " + str(self.angle) + "\n")

    
main()