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
CNlist = []
OHlist = []
bondnum = 0
start = False

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
currStep = -1
for line in lineFile:       #loops through each line of the file
    stored = False
    wordList = line.split()     #splits the line into a list of words dilineated by any space
    if (wordList[0] == '#' and wordList[1] == 'Timestep' and int(wordList[2]) == currentStep):      # the hashtag starts the header area
        int(currStep) = wordList[2]      #grabs and stores timestep
        start = True
    if (wordList[0] != '#' and start):      #into the body of the file
        if (currStep == 0):
            findPairs(lineFile,Clist,Olist,Nlist,Hlist)
        if (wordList[0] in Clist):      #THis section checks for N-C bonds, if first is C
            for group in bondList:
                if (wordList[0] in bondList[group]):        #and we haven't already stored this C
                    stored = True
                    break
            if (not stored):
                bondnum = wordList[3]       #gets number of bond this atom has.
                for i in range(bondnum):
                    if (wordList[3 + i] in Nlist):
                        CNlist.append([wordList[0],wordList[3+i],currStep])        #adds C and N to an array as a group.  
        elif (wordList[0] in Olist):
            for group in bondList:
                if (wordList[0] in bondList[group]):        #and we haven't already stored this N
                    stored = True
                    break
            if (not stored):
                bondnum = wordList[3]
                for i in range(bondnum):
                    if (wordList[3 + i] in Hlist):
                        OHlist.append([wordList[0],wordList[3+i], currStep])
            
start = False           
                        
def getCO (Clist,Olist,wordList):    
    #this should pull CO pairs from timestep 0.
    #to be referenced and called in getBonds to link NC to OH once bonds formed.
    #mods COlist  (stored [C,O] by index)
    if (int(wordList[0]) in Clist):
        bondnum = wordList[3]       #gets number of bond this atom has.
        for i in range(bondnum):
            if (wordList[3 + i] in Olist):
                COlist.append([wordList[0],wordList[3 + i]])
        
    

def getNH(Nlist,Hlist,wordList):
    #this should pull NH pairs from timestep 0
    if (int(wordList[0]) in Nlist):
        bondnum = wordList[3]       #gets number of bond this atom has.
        for i in range(bondnum):
            if (wordList[3 + i] in Hlist):
                NHlist.append([wordList[0],wordList[3 + i]])
                
def findPairs(lineFile,Clist,Olist,Nlist,Hlist):
    #calls getCO and getNH on every line in the first timestep
    for line in lineFile:       #loops through each line of the file
        wordList = line.split()
        if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
            currStep = int(wordList[2]) 
        if (currStep == 0):
            getCO(Clist,Olsit,wordList)
            getNH(Nlist,Hlist,wordList)
        elif(currStep > 0):
            break
badReac = False;           
def mergeCONH():
    flag = True
    for CN in CNlist:
        for CO in COlist:
            if (CN[0] == CO[0]):
                for NH in NHlist:
                    if (CN[1] == NH[0]):
                        
                
Clist = [71,121,120,72]

Olist = [122,74,121,73]

Nlist = [9,8,89,90
        
            
#def findCH(Clist, CHlist):
#    #this uses special H identifier to help keep the H on the active C by applying restraint force. We still need to think of the parameters.
#    for line in lineFile:       #loops through each line of the file
#        wordList = line.split()     #splits the line into a list of words dilineated by any space
#        if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
#            currStep = wordList[2] 
#        elif (currStep == 0):
#            CHlist.append([wordList[0],wordList[3 + i]])
