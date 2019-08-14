# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 10:24:50 2019

@author: mattcohe
"""

#tar are ID's of first molecule atoms to change. Must include number of tar, starting at 1.
#also must modify section of ifs to add to mol size.
tar0 = 506
tar1 = 507
tar2 = 520
tar3 = 521
numTar = 4
newfArray = []
#numMol is how many molecules you want to change.  
#molsize is how many atoms in a molecule
#change to is what you want the type to be changed to.
numMol = 32
molsize = 49
currMol = 0
end = numMol * numTar
curr = 0
changeTo0 = 5    
changeTo1 = 6    
changeTo2 = 5    
changeTo3 = 6    

f = open("32Epon-16DETDA-MOD0.txt", "r")
lineFile = f.readlines()
newFile = open("32Epon-16DETDA-MOD1.txt","w+")
tar0-=1
tar1-=1
tar2-=1
tar3-=1

for line in lineFile:       #loops through each line of the file
    wordList = line.split()
    if ((curr == tar0 or curr == tar1 or curr == tar2 or curr == tar3) and (currMol < end)):
        newfArray.append(wordList)
        if (curr == tar0):
            wordList[2] = str(changeTo0)
            tar0 += molsize
        if (curr == tar1):
            wordList[2] = str(changeTo1)
            tar1 += molsize
        if (curr == tar2):
            wordList[2] = str(changeTo2)
            tar2 += molsize
        if (curr == tar3):
            wordList[2] = str(changeTo3)
            tar3 += molsize
        currMol +=1
    else:
        newfArray.append(wordList)
    curr+=1
    
#print(newfArray)
for line in newfArray:
    newFile.write(str(line[0]) + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(line[3]) + "\t" + str(line[4]) + " " + str(line[5]) + " " + str(line[6]) + " ")
    newFile.write("\n")
newFile.close()