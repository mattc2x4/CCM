# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 09:15:43 2019
        this is just so bad but within what I have it has to be lmao wtf
@author: mattcohe
"""
natoms = 10
atomType = []
Clist = []
Olist = []
Hlist = []
Nlist = []
bondList = []
bondnum = 0
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
    wordList = line.split()     #splits the line into a list of words dilineated by any space
    if (wordList[0] == '#' and wordList[1] == 'Timestep'):      # the hashtag starts the header area
        currStep = wordList[2]      #grabs and stores timestep
    elif (wordList[0] != '#'):
        if (wordList[0] in Clist):      #THis section checks for N-C bonds, if first is C
            for group in bondList:
                if (wordList[0] not in bondList[group]):        #and we haven't already stored this group
                    bondnum = wordList[3]
                    for i in range(bondnum):
                        if (wordList[3 + i] in Nlist):
                            bondList.append([wordList[0],0,wordList[3+i],0])        #adds C and N to an array as a group.  zeroes are placeholders for H and O, and should be in the form [C,O,N,H] to follow previous group notation.
            
           
                    