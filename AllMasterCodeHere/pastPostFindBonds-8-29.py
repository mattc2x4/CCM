# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 14:00:41 2019

@author: mattcohe
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 15:39:59 2019
@author: mattcohe
"""
"""
Call Heirarchy:
    #todo
brief explanation: 
    first we get pairs.  this records NH and CO groups from model, calls getCO and getNH. this is important to make sure the OH and CN bonds are with the correct atoms.
    then we getFormed. this finds OH and CN pairs and adds them to a list. 
    then we call getH2O, which finds H2O on last timestep lol. also calls checkCO to save another loop and for convenience. this finds secondary C bonds with formed OH groups. 
    similiar function to get H2O but is very indicative of incorrect bonds.
    mergeCONH mkaes checks to insure correct H bonds with correct O, and correct C N. insures CO-NH has same atoms as CN-OH. also removes things in checkCO.
    
    t
"""

"""
TODO: add unmade bond counter: NH and CO still fully intact. in this case only more time needed/better locations
"""


NHlist = []
COlist = []
restID = []
natoms = 10
atomType = []
Clist = []
Olist = []
Hlist = []
Nlist = []
bondList = []
CNlist = []
OHlist = []
bondnum = 0
currStep = -1
Ctype = 5
Otype = 6
Htype = 7
Ntype = 3
pastStep = -1
H2Olist = []
H1list = [] #non active H, for checking for H2O
H1type = 4
Hcount = 0
C1type = 1
C1list = [] #every nonactive carbon
realC1list = [] # the non active carbons bonded to active O's
greatestStep = [-1]
OHremove = []
removeGroup = []
origNC = []
currNC = []



def main():
    getCONH_fromMODEL('8Epon-4DETDA-H.txt')
    #getCONH_fromMODEL('4Epon-2DETDA-Packmol-H.txt')
    #getCONH_fromMODEL('2Epon-1DETDA-H.txt')

    
    #stores all values by index not by ID
    f = open("bonds.txt", "r")
    lineFile = f.readlines()
    getPairs(lineFile,Clist,Olist,Nlist,Hlist)
#    print("Clist: " + str(Clist))
#    print("Olist: " + str(Olist) + " len " + str(len(Olist)))
#    print("Hlist: " + str(Hlist) + " len " + str(len(Hlist)))
#    print("Nlist: " + str(Nlist) + " len " + str(len(Nlist)))
    #print("C1list: " + str(C1list))
    #print("NHList: " + str(NHlist) + " len " + str(len(NHlist)))
    #print("COlist: " + str(COlist) + " len " + str(len(COlist)))
    #print("realC1list: " + str(realC1list) + "len " + str(len(realC1list)))
    getFormed(lineFile)
    getH2O(lineFile)
    print('OHlist: ' + str(OHlist) + " len: "  + str(len(OHlist)))
    print('CNlist: ' + str(CNlist) + str(len(CNlist)))
    print("H2Olist: " + str(H2Olist) + " len: " + str(len(H2Olist)))
    print("OHremove: " + str(OHremove)  + " len " + str(len(OHremove)))
    mergeCONH()
    print("BondList: " + str(bondList)+ " len: " + str(len(bondList)))
    print("bonds: "  + str(len(bondList) ))
    print("H2O: "  + str(len(H2Olist) ))
    print("OH: "  + str(len(OHlist) ))
    print("CN: "  + str(len(CNlist) ))
    #print("origNC: "  + str(origNC))
    #print("currNC: "  + str(currNC))
    print("OHremove: "  + str(len(OHremove) ))
    crossCheckNC()
    crossCheckOHNC()

    
    



    #print(H1list)
    
    


          
def getCO (Clist,Olist,wordList):  
    #print('in getCO on ' + str(wordList))
    #this should pull CO pairs from timestep 0.
    #to be referenced and called in getBonds to link NC to OH once bonds formed.
    #mods COlist  (stored [C,O] by index)
    #also adds all non active C in epoxide to C1list.
    add = True
    if (int(wordList[0]) in Clist):
        for CO in COlist:
            if (int(wordList[0]) in CO):
                add = False
      #  print("found C " + str(wordList))
        bondnum = int(wordList[2])       #gets number of bond this atom has.
        for i in range(bondnum):
         #   print("bondnum = " + str(bondnum))
            if (int(wordList[3 + i]) in Olist and add):
                COlist.append([int(wordList[0]),int(wordList[3 + i])])
                #print("CO added: " + str([wordList[0],wordList[3 + i]]))
            if((int(wordList[3+i]) in C1list) and (int(wordList[3+i]) not in realC1list)):
                realC1list.append(int(wordList[3+i]))

def getNC(wordList,currStep):     #get origional C that N is bonded to, so that we can make sure they're still bonded later on. 
    add = True
    if (int(wordList[0]) in Nlist):
        for NC in origNC:
            if (int(wordList[0]) in NC):
                add = False
      #  print("found C " + str(wordList))
        bondnum = int(wordList[2])       #gets number of bond this atom has.
        for i in range(bondnum):
         #   print("bondnum = " + str(bondnum))
            if (int(wordList[3 + i]) in C1list and add):
                origNC.append([int(wordList[0]),int(wordList[3 + i])])
               
    

def getNH(Nlist,Hlist,wordList):
    #print('in getNH on ' + str(wordList))
    #this should pull NH pairs from timestep 0
    add = True
    
    if (int(wordList[0]) in Nlist and add):
        #print("found N " + str(wordList))
        bondnum = int(wordList[2])       #gets number of bond this atom has.
        for i in range(bondnum):
            if (int(wordList[3 + i]) in Hlist):
                    for NH in NHlist:
                        if (int(wordList[3+i]) in NH):
                            add = False
                    if (add):
                        NHlist.append([int(wordList[0]),int(wordList[3 + i])])
                                #print("NH added: " + str([wordList[0],wordList[3 + i]]))
                
def getPairs(lineFile,Clist,Olist,Nlist,Hlist):
    #calls getCO and getNH on every line in the first timestep
    currStep = -1
    print("in findPairs")
    for line in lineFile:       #loops through each line of the file
        wordList = line.split()
        if (len(wordList) > 2):
            if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
                currStep = int(wordList[2]) 
            elif (currStep == 0 and wordList[0] != '#'):
                getCO(Clist,Olist,wordList)
                getNH(Nlist,Hlist,wordList)
                getNC(wordList,currStep)
            elif(currStep > 0):
                break


def mergeCONH():
    print('in mergeCONH')
    for CN in CNlist:       #get CN
        for CO in COlist:       #find Carbon initial bond
            if (CN[0] == CO[0]):        #Match C between CN and CO
                for NH in NHlist:       #looking for H linked with N
                    if (CN[1] == NH[0]): #Match N between CN and NH
                        for OH in OHlist:                            
                            if (NH[1] == OH[1] and CO[1] == OH[0]):     #Match H between NH and OH, and O between OH and CO
                                bondList.append([CN[0], OH[0],CN[1],OH[1], OH[2]])     # to get here, C same between CN and CO, N same between NH and CN, H same between NH and OH, and O same between OH and CO.  Includes Timestep.
    #print(bondList)
    for group in bondList:
        if (group[1] in OHremove or group[1] in H2Olist ):
            removeGroup.append(group)
    #print(bondList)
    for group in removeGroup:
        bondList.remove(group)
        #print("removing " + str(group) + " because " + str(group[1]) + " in OHremove")

def getCONH_fromSIM():
    for i in range(natoms):     #gets all carbon locations.  won't be ID till end. 
        if (atomType[i] == Ctype):
            Clist.append(i)
        elif (atomType[i] == Otype):
            Olist.append(i)
        elif (atomType[i] == Htype):
            Hlist.append(i)
        elif (atomType[i] == Ntype):
            Nlist.append(i)
        elif (atomType[i] == H1type):
            H1list.append(i)
        elif (atomType[i] == C1type):
            H1list.append(i)
            
def getCONH_fromMODEL(modelName):
    f = open(modelName, "r")
    lineModel = f.readlines()
    for line in lineModel:       #loops through each line of the file
        wordList = line.split() 
        #print(wordList)
        #todo: skip header:: meh yay done
        if (len(wordList) == 7): 
            #print(wordList[2])
            if (int(wordList[2]) == Ctype):
                Clist.append(int(wordList[0]))
            elif (int(wordList[2]) == Otype):
                Olist.append(int(wordList[0]))
            elif (int(wordList[2]) == Htype):
                Hlist.append(int(wordList[0]))
            elif (int(wordList[2]) == Ntype):
                Nlist.append(int(wordList[0]))
            elif (int(wordList[2]) == H1type):
                H1list.append(int(wordList[0]))
            elif (int(wordList[2]) == C1type):
                C1list.append(int(wordList[0]))
            
            
def getFormed(lineFile):
    mostCurrent = getMostCurrentTimeStep(lineFile)
    for line in lineFile:       #loops through each line of the file
       wordList = line.split()     #splits the line into a list of words dilineated by any space
       #print(wordList)
       if (len(wordList) > 2):
           if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
               currStep = int(wordList[2])      #grabs and stores timestep
           elif (wordList[0] != '#'):
                 #this should go through entire file and get all formed OH and NC bonds. 
                 if (mostCurrent == currStep):
                     getCN(wordList,currStep)
                     getOH(wordList,currStep)
                     checkNC(wordList)


    
    
def getCN(wordList,currStep):
    stored = False
    if (int(wordList[0]) in Clist):      #THis section checks for N-C bonds, if first is C
        for group in CNlist:
            if (int(wordList[0]) in group):        #and we haven't already stored this C
                stored = True
                break
        if (not stored):
            bondnum = int(wordList[2])       #gets number of bond this atom has.
            for i in range(bondnum):
                if (int(wordList[3 + i]) in Nlist):
                    CNlist.append([int(wordList[0]),int(wordList[3+i]),currStep])        #adds C and N to an array as a group.  
                            
def getOH(wordList,currStep):
    stored = False
    if (int(wordList[0]) in Olist):
        for group in OHlist:
            if (int(wordList[0]) in group):        #and we haven't already stored this O
                stored = True
                break
        if (not stored):
            
            bondnum = int(wordList[2])
            for i in range(bondnum):
                if (int(wordList[3 + i]) in Hlist):
                    OHlist.append([int(wordList[0]),int(wordList[3+i]), currStep])

def getMostCurrentTimeStep(lineFile):
    currStep = -1
    greatestStep = 0
    for line in lineFile:       #loops through each line of the file
        wordList = line.split() 
        if (len(wordList) > 2):
            if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
                currStep = int(wordList[2])
                if (greatestStep < currStep):
                    greatestStep = currStep
    return greatestStep
   

def getH2O(lineFile):
    #records the O in H2O.  this O should be removed from the bondlist in mergeCONH
    latestStep = getMostCurrentTimeStep(lineFile)
    for line in lineFile:       #loops through each line of the file
        wordList = line.split()     #splits the line into a list of words dilineated by any space
        Hcount = 0
        if (len(wordList) > 2):
            if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
                currStep = int(wordList[2])      #grabs and stores timestep
            elif (wordList[0] != '#' and latestStep == currStep):
                stored = False                  
                if (int(wordList[0]) in Olist):
                    checkCO(wordList)
                    
                    if (int(wordList[0]) in H2Olist):        #and we haven't already stored this O
                        stored = True
                        break
                    if (not stored):
                        bondnum = int(wordList[2])
                        Hcount = 0
                        for i in range(bondnum): 
                            if (int(wordList[3 + i]) in Hlist or int(wordList[3 + i]) in H1list):
                                Hcount+=1
                                if (Hcount>1):
                                    H2Olist.append(int(wordList[0]))

def crossCheckOHNC():     #will print any OH and CN group that is not recorded in bondList.  NC part not so good becuase its hard to know how many bonds it should have.  
    for OH in OHlist:
        for group in bondList:
            stored = False
            if (group[1] == OH[0]):
                stored = True
                break
        if (not stored and OH[0] not in OHremove):
            print(str(OH) + " OH not Correct")
    for CN in CNlist:
        for group in bondList:
            stored1 = False
            if (group[0] == CN[0]):
                stored1 = True
                break
        if (not stored1):
            print(str(CN) + " CN not Correct")
            
            
def crossCheckNC(): #prints the bonds in origNC that are no longer bonded correctly. 
    for NC in origNC:
        for group in currNC:
            stored = False
            if (group[0] == NC[0]):
                stored = True
                break
        if (not stored):
            print(str(NC) + " CN no longer bonded")

def checkNC(wordList):  #checks afterwards that all the NC bonds are still present, appends to currNC which is then checked in crossCheckNC
    remain = False
    if (int(wordList[0]) in Nlist):
        for NC in origNC:
            if (NC[0] == int(wordList[0])):
                bondnum = int(wordList[2])
                for i in range(bondnum): 
                    if (int(wordList[3 + i]) == NC[1]):
                        remain = True
                if (remain):
                    currNC.append(NC)
        
                
    

def checkCO(wordList):
    #should be called on oxygen lines being added to OHlist.
    stored = False
    append = True
    if (int(wordList[0]) in OHremove):      
        stored = True
    if (not stored):
        bondnum = int(wordList[2])
        for i in range(bondnum): 
            if (int(wordList[3 + i]) in C1list and int(wordList[3 + i]) in OHremove):
                OHremove.remove(int(wordList[3+i]))
                append = False
            elif(int(wordList[3 + i]) in C1list):
                append = False

        if (append):
            OHremove.append(int(wordList[0]))
           # print(str(wordList[0]) + " not bonded to C1" )
    
main()