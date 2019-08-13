# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 09:15:43 2019
        this is just so bad but within what I have it has to be
@author: mattcohe
"""
#stores all values by index not by ID
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
bondnum = 0
stored = False
for i in range(natoms):     #gets all carbon locations.  won't be ID till end. 
    if (atomType[i] == Ctype):
        Clist.append(i)
    elif (atomType[i] == Otype):
        Olist.append(i)
    elif (atomType[i] == Htype):
        Hlist.append(i)
    elif (atomType[i] == Ntype):
        Nlist.append(i)
totalBond = 0
f = open("testbonds.txt", "r")
lineFile = f.readlines()
formedBonds = []
currStep = 0
for line in lineFile:       #loops through each line of the file
    stored = False
    wordList = line.split()     #splits the line into a list of words dilineated by any space
    if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
        currStep = wordList[2]      #grabs and stores timestep
    elif (wordList[0] != '#'):      #into the body of the file
        findPairs()
        if (wordList[0] in Clist):      #THis section checks for N-C bonds, if first is C
            for group in bondList:
                if (wordList[0] in bondList[group]):        #and we haven't already stored this C
                    stored = True
                    break
            if (not stored):
                bondnum = wordList[3]       #gets number of bond this atom has.
                for i in range(bondnum):
                    if (wordList[3 + i] in Nlist):
                        bondList.append([wordList[0],-1,wordList[3+i],-1])        #adds C and N to an array as a group.  -1 are placeholders for O and H, and should be in the form [C,O,N,H] to follow previous group notation.
        elif (wordList[0] in Nlist):
            for group in bondList:
                if (wordList[0] in bondList[group]):        #and we haven't already stored this C
                    stored = True
                    break
            if (not stored):
                bondnum = wordList[3]
                for i in range(bondnum):
                    
            
            
                        
def getCO (Clist,Olist,wordList):    
    #this should pull CO pairs from timestep 0.
    #to be referenced and called in getBonds to link NC to OH once bonds formed.
    #mods COlist  (stored [C,O] by index)
    if (wordList[0] in Clist):
        bondnum = wordList[3]       #gets number of bond this atom has.
        for i in range(bondnum):
            if (wordList[3 + i] in Olist):
                COlist.append([wordList[0],wordList[3 + i]])
        
    

def getNH(Nlist,Hlist,wordList):
    #this should pull NH pairs from timestep 0
    if (wordList[0] in Nlist):
        bondnum = wordList[3]       #gets number of bond this atom has.
        for i in range(bondnum):
            if (wordList[3 + i] in Hlist):
                NHlist.append([wordList[0],wordList[3 + i]])
                
def findPairs(lineFile,Clist,Olist,Nlist,Hlist,wordList):
    #calls getCO and getNH on every line in the first timestep
    for line in lineFile:       #loops through each line of the file
        stored = False
        wordList = line.split()
        if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
            currStep = wordList[2] 
        elif (currStep == 0):
            getCO(Clist,Olsit,wordList)
            getNH(Nlist,Hlist,wordList)
            
            
def findCH(Clist, CHlist):
    #this uses special H identifier to help keep the H on the active C by applying restraint force. We still need to think of the parameters.
    for line in lineFile:       #loops through each line of the file
        wordList = line.split()     #splits the line into a list of words dilineated by any space
        if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
            currStep = wordList[2] 
        elif (currStep == 0):
            