# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 15:39:59 2019

@author: mattcohe
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

def main():
    print('in main')
    getCONH_fromMODEL('4Epon-2DETDA-Packmol-H.txt')
    print("Clist: " + str(Clist))
    print("Olist: " + str(Olist))
    print("Hlist: " + str(Hlist))
    print("Nlist: " + str(Nlist))
    #stores all values by index not by ID
    f = open("bonds.txt", "r")
    lineFile = f.readlines()
    getPairs(lineFile,Clist,Olist,Nlist,Hlist)
    print("NHList: " + str(NHlist))
    print("COlist: " + str(COlist))
    getFormed(lineFile)
    print('OHlist: ' + str(OHlist))
    print('CNlist: ' + str(CNlist))
    mergeCONH()
    print("BondList: " + str(bondList))

    

          
def getCO (Clist,Olist,wordList):  
    #print('in getCO on ' + str(wordList))
    #this should pull CO pairs from timestep 0.
    #to be referenced and called in getBonds to link NC to OH once bonds formed.
    #mods COlist  (stored [C,O] by index)
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
            
def getCONH_fromMODEL(modelName):
    print('in get from model')
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
            if (int(wordList[2]) == Otype):
                Olist.append(int(wordList[0]))
            if (int(wordList[2]) == Htype):
                Hlist.append(int(wordList[0]))
            if (int(wordList[2]) == Ntype):
                Nlist.append(int(wordList[0]))
            
            
def getFormed(lineFile):
       for line in lineFile:       #loops through each line of the file
        wordList = line.split()     #splits the line into a list of words dilineated by any space
        #print(wordList)
        if (len(wordList) > 2):
            if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
                currStep = int(wordList[2])      #grabs and stores timestep
            elif (wordList[0] != '#'):
                stored = False
                #this should go through entire file and get all formed OH and NC bonds. 
                if (len(wordList) > 2):
                    if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
                        currStep = int(wordList[2])
                if (wordList[0] != '#'):      #into the body of the file
                    if (currStep != 0):
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

                        elif (int(wordList[0]) in Olist):
                            for group in OHlist:
                                if (int(wordList[0]) in group):        #and we haven't already stored this O
                                    stored = True
                                    break
                            if (not stored):
                                bondnum = int(wordList[2])
                                for i in range(bondnum):
                                    if (int(wordList[3 + i]) in Hlist):
                                        OHlist.append([int(wordList[0]),int(wordList[3+i]), currStep])

    
def getMostCurrentTimeStep(bondFile):
    f = open("testbonds.txt", "r")
    lineFile = f.readlines()
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
main()