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
currStep = -1
Ctype = 5
Otype = 6
Htype = 7
Ntype = 3
currentStep = 0
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
        currStep = int(wordList[2])      #grabs and stores timestep
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
    for CN in CNlist:       #get CN
        for CO in COlist:       #find Carbon initial bond
            if (CN[0] == CO[0]):        #Match C between CN and CO
                for NH in NHlist:       #looking for H linked with N
                    if (CN[1] == NH[0]): #Match N between CN and NH
                        for OH in OHlist:                            
                            if (NH[1] == OH[1] and CO[1] == OH[0]):     #Match H between NH and OH, and O between OH and CO
                                bondList.append([CN[0], OH[0],CN[1],OH[1]])     # to get here, C same between CN and CO, N same between NH and CN, H same between NH and OH, and O same between OH and CO.
                            
                
                
Clist = [71,121,120,72]

Olist = [122,74,121,73]

Nlist = [9,8]
        
def main():
    NH = [[8,19],[8,20]]
    CO = [[1,2],[3,4]]
    CN = [[1,8],[3,8]]
    OH = [[2,19],[4,20]]
    mergeCOHN()
main()           
            
#def findCH(Clist, CHlist):
#    #this uses special H identifier to help keep the H on the active C by applying restraint force. We still need to think of the parameters.
#    for line in lineFile:       #loops through each line of the file
#        wordList = line.split()     #splits the line into a list of words dilineated by any space
#        if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
#            currStep = wordList[2] 
#        elif (currStep == 0):
#            CHlist.append([wordList[0],wordList[3 + i]])