# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 08:52:49 2019

@author: mattcohe
"""
i = 0
j =0
k = 0
m = 0
newfile = open("rest-ALLdata.txt",'a')
newfile.write("\nC\n ID: " + str(i+1) + "X: "+ str(coordinates[i*3]) + " Y: " str(coordinates[i*3+1])+ " Z: " str(coordinates[i*3+2]))
