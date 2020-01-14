# -*- coding: utf-8 -*-
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.stats
import seaborn as sns
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

vertexType = 1   #this is where you should put the type of the vertex
endType1 = 3    #the other two types you want to calculate angle of
endType2 = 3

#O-Al-O


currStep  = 4080000   #put the first step here
finalStep = 4580000    #put the final Step here
incrSize = 5000   ## of timesteps inbetween each print in dump/ bond file

"""end Constants"""


atomList = []   #this will hold all angle data types. will be in order, with the 0th atom being the atom with ID 1.
angleVals = [] #containts angle values in raw form, to be graphed
angList = []  #this will hold all angle vals calculated below, in the format [[vertID,endID1,endID2, ANGLE],...] endID1 and 
#clear angList every timestep.
#end2ID interchangable locations. 

boxdim = [0,0,0]     #[x,y,z] lengths

def main():
    currStep  = 4080000
    #CONFIGURE FILES
    #INSERT DUMP_FINAL FILES TO DUMP
    #INSERT BOND FILES TO BONDS
    i = 0
    while(currStep < finalStep):
        dump = open("dump_final.lammps","r") 
        bonds = open("MD_bonds_final.reaxc", "r")
        print("\ncurrStep: " + str(currStep))
        fillAtomList(dump,currStep)
        print("AtomList[0]: " + str(atomList[0]))
        getAngleID(bonds,currStep)
        calcAngles(currStep)
        print("average value: " + str(sum(angleVals) / len(angleVals)))
        currStep += incrSize
        angList.clear()
        i+=1
        if(i>1):
            break
        dump.close()
        bonds.close()
    plot(angleVals,"sim")
    print(len(angleVals))
    # print("currStep: " + str(currStep))
    # fillAtomList(dump,currStep)
    # getAngleID(bonds,currStep)
    # calcAngles(currStep)
    # print("average value: " + str(sum(angleVals) / len(angleVals)))
    # currStep += incrSize
    # angList.clear()



def fillAtomList(dump,currStep):
    #read in the atom data from timestep "currStep"
    #finds meantions of timestep in header. compares with input (master) step. if it is greater than, breaks loop, allowing
    #main code to run instead. when it finds the correct 
    #updates Boxdim, which changes the box dimensions. this affects finding distance with periodic boundaries. 
    print("filling atom list for " + str(currStep) + " step.")
    startRead = False
    readBox = False
    step = False
    dumpLine = dump.readlines()
    print(len(dumpLine))
    for i in range(len(dumpLine)):
        dumpWord = dumpLine[i].split()
        if(dumpWord[0] == "ITEM:" and dumpWord[1] == "TIMESTEP"):    #read timestep, if its correct continue. otherwise, break loop.
            print("found TIMESTEP header")
            myStep = int(dumpLine[i+1].split()[0])
            if(myStep > currStep):
                print("breaking at " + str(myStep))
                break
            elif(myStep == currStep):
                step = True
                print("found correct Timestep: " + str(myStep))
            #print("myStep : " + str(myStep) + " currStep: " + str(currStep))
        if(dumpWord[0] == "ITEM:" and dumpWord[1] == "BOX" and step and not readBox):    #UPDATING BOXDIM.  havent read box yet this step,and have read in the timestep.
            print("found BOX header.")
            readBox = True
            currLine = i
            box1 = dumpLine[currLine + 1].split()   #getting three lines which contain X lo hi, y lo hi, z lo hi for boxdim.
            box2 = dumpLine[currLine + 2].split()   #Y
            box3 = dumpLine[currLine + 3].split()   #Z
            boxdim[0] = float(box1[1]) - float(box1[0])
            boxdim[1] = float(box2[1]) - float(box2[0])
            boxdim[2] = float(box3[1]) - float(box3[0])
            print("updating boxdim: " + str(boxdim))
        if(dumpWord[0] == "ITEM:" and dumpWord[1] == "ATOMS"):    #when we see this header we want to skip this iteration, then continue on the next line, hence continue. 
            print("found ATOMS header")
            startRead = True
            continue
        if (startRead and step):
            atomList.append(Atom(float(dumpWord[3]), float(dumpWord[4]), float(dumpWord[5]), int(dumpWord[0]), int(dumpWord[1])))
            #this adds an atom object in atomList. x,y,z,ID,TYPE
            

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
    print("getting angles for " + str(currStep) + " step.")
    lineFile = bonds.readlines()
    read = False
    for line in lineFile:       #loops through each line of the file
        wordList = line.split()
        endCount = 0
        if (len(wordList) > 2):
            if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # this checks timestep. if timestep in file matches input, read lines.upon encountering another, break loop.
                if(currStep == int(wordList[2])):
                    read = True
                else:
                    #print("breaking\n")
                    break
            if (read == True and wordList[0] != '#'):
                #print(str(atomList[int(wordList[0]) - 1].TYPE) + str(atomList[int(wordList[0]) - 1].TYPE == vertexType))
                if (atomList[int(wordList[0]) - 1].TYPE == vertexType):
                    #access the corresponding atom data. subtract one to translate to index. Make sure ID is correct.  and int(wordList[0]) == atomList[int(wordList[0]) - 1].ID
                    bondnum = int(wordList[2])       #gets number of bond this atom has.
                    for i in range(bondnum):
                        try: 
                            atomList[int(wordList[3 + i]) - 1].TYPE  == endType1 or atomList[int(wordList[3 + i]) - 1].TYPE == endType2
                        except IndexError:
                            print("Files Incomplete: an atom exists in bonds file that does not exist in dump this atom is: " + str(int(wordList[3 + i])))
                        if (atomList[int(wordList[3 + i]) - 1].TYPE  == endType1 or atomList[int(wordList[3 + i]) - 1].TYPE == endType2):      
                           endCount+=1
                           #print("found end")
                           if (endCount == 1):     #this is the first one we found, store in temp
                                firstEndID = int(wordList[3+i])
                           elif(endCount == 2):    #this is the second one we found, add time, reset counter.
                                angList.append(Angle(int(wordList[0]),firstEndID,int(wordList[3+i]),currStep,-1))
                                endCount = 0

def calcAngles(currStep):
    #this function will take the info in angList, and use it to calculate angle values based on the atom data in atomList.
    #called on angLIst once per timestep, though angList will contain all steps. may be problematic due to data sizes. 
    #TODO: modify this to decrease compexity based on what yoon wants. this dude needs to respond faster.
    print("calcing Angles for: " + str(currStep) + " step.")
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
    #plt.show()
           
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