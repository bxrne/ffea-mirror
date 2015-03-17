#!/usr/bin/env python

import os, sys
import math

if len(sys.argv) != 3 and len(sys.argv) != 5 :
	sys.exit("Usage: python translate_surf [INPUT TRANSLATING .SURF FILE] [INPUT TEMPLATE .SURF FILE]\nor\tpython translate_surf [INPUT .SURF FILE] [X TRANSLATION] [Y TRANSLATION] [Z TRANSLATION]")

infile1 = open(sys.argv[1], "r")
outfile = open(sys.argv[1].split(".")[0] + "_translated.surf", "w")

# Explicit traslation
if len(sys.argv) == 5:
	x_trans = float(sys.argv[2])
	y_trans = float(sys.argv[3])
	z_trans = float(sys.argv[4])

# Translating Surface to fit on to of another surf file
if len(sys.argv) == 3:
	x_trans = 0.0
	y_trans = 0.0
	z_trans = 0.0
	centroid1_x = 0.0
	centroid1_y = 0.0
	centroid1_z = 0.0
	centroid2_x = 0.0
	centroid2_y = 0.0
	centroid2_z = 0.0

	infile2 = open(sys.argv[2], "r")

	#calculating centroids	
	if(infile1.readline()[0:11] != "surfacemesh"):
		sys.exit("Error: Input file 1 is not a netgen .surf file.")
	num_nodes1 = int(infile1.readline())
	for i in range(0, num_nodes1):
		line = infile1.readline().split()
		centroid1_x = centroid1_x + float(line[0])
		centroid1_y = centroid1_y + float(line[1])
		centroid1_z = centroid1_z + float(line[2])

	centroid1_x = centroid1_x * 1.0 / num_nodes1
	centroid1_y = centroid1_y * 1.0 / num_nodes1
	centroid1_z = centroid1_z * 1.0 / num_nodes1
	print centroid1_x, centroid1_y, centroid1_z

	if(infile2.readline()[0:11] != "surfacemesh"):
		sys.exit("Error: Input file 2 is not a netgen .surf file.")
	num_nodes2 = int(infile2.readline())
	for i in range(0, num_nodes2):
		line = infile2.readline().split()
		centroid2_x = centroid2_x + float(line[0])
		centroid2_y = centroid2_y + float(line[1])
		centroid2_z = centroid2_z + float(line[2])

	centroid2_x = centroid2_x * 1.0 / num_nodes2
	centroid2_y = centroid2_y * 1.0 / num_nodes2
	centroid2_z = centroid2_z * 1.0 / num_nodes2

	x_trans = centroid2_x - centroid1_x
	y_trans = centroid2_y - centroid1_y
	z_trans = centroid2_z - centroid1_z	

# Outputting translated file
infile3 = open(sys.argv[1], "r")
if(infile3.readline()[0:11] != "surfacemesh"):
	sys.exit("Error: Input file 1 is not a netgen .surf file.")

outfile.write("surfacemesh\n")
num_nodes = int(infile3.readline())
outfile.write(str(num_nodes) + "\n")
for i in range(0, num_nodes):
	line = infile3.readline().split()
	outfile.write(str(float(line[0]) + x_trans) + " " + str(float(line[1]) + y_trans) + " " + str(float(line[2]) + z_trans) + "\n")
num_faces = int(infile3.readline())
outfile.write(str(num_faces) + "\n")
for i in range(0, num_faces):
	outfile.write(infile3.readline())	
