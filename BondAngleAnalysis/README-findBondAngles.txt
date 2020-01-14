#README:
To use findBondAngles.py:

put this file in a directory with the listed files. 

	put vertex, and end types in line 20, 21, 22. these should correspond with 
	the number types in the files. 

	Put current step, final step, and increment size on lines 27, 28, 29.

	finally, put file names inside open() function on line 50, 51.  dmp should correspond to 
	dump final, and bonds should correspond to the bond file. 
		the format should look like: dump = open("dump_final.lammps","r")
							  ^modify file name here, place inside quotes.