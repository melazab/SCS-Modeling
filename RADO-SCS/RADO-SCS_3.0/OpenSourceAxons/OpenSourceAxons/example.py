''' 
This example file will load the X,Y and Z center coordinates for an axon in the RADO-SCS model, as well as the 
diameter of the fiber.

This file was written by Hans Zander, a PhD candidate at The University of Michigan
email: hzander@umich.edu

Updated by Mohamed Elazab (07/10/2026)
email: mohamed.elazab@case.edu
'''


### Import the correct packages
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from pathlib import Path

# Define a function to load both the xyz center coordinates of the axon, as well as the diameter
def loadAxon(axonNumber,path2Files):
	'''
	Inputs:
		axonNumber	-	int (between 0 and 2039)		Which axon (out of 2039) should be loaded
		path2Files	-	str								The path to the folder containing both the axon files, 
														as well as the diameter file

	Outputs:
		fiberD 		-	float		[um]				The diameter of the fiber in um (5.7, 7.3, 8.7, 10.0 or 11.5)
		axonXYZ 	-	Nx3 array	[m]					The x (column 1), y (col2) and z (col3) center locations
														of each compartment in the axon in meters. 
														The locations are sequential, and have the following structure:

														row0 	Node 0
														row1 	MYSA 0
														row2 	FLUT 0
														row3-8	STIN 0-5
														row9	FLUT 1
														row10	MYSA 1
														row11	Node 1
														row12	MYSA 2
														...
														rowN	Node N

														The files will both begin and end on a node.

														McIntyre, Cameron C., Andrew G. Richardson, and Warren M. Grill. 
														"Modeling the excitability of mammalian nerve fibers: influence 
														of afterpotentials on the recovery cycle." Journal of 
														neurophysiology 87.2 (2002): 995-1006.
	'''

	### Load the axon files
	AxonDiameters = np.load(path2Files / 'AxonDiameters.npy')
	axonXYZ = np.load(path2Files / ('axon_' + str(axonNumber) + '.npy'))

	fiberD = AxonDiameters[axonNumber]

	return fiberD, axonXYZ


### Set the path to the file

path2Files = Path.home() / 'Box'/'Electrical Nerve Block Institute'/ 'Admin'/'Grant Applications' / 'NIH Grants' / '2026 R01 SCS MoA' / 'Modeling' / 'RADO-SCS' / 'RADO-SCS_3.0' / 'OpenSourceAxons'/ 'OpenSourceAxons'

# -------------------------- Example 1 ---------------------------------
# Determine the diameter and plot the Y,Z trajectory of axon 2038
fiberD, axonXYZ = loadAxon(2038,path2Files)

print('The diameter of fiber 2038 is ' + str(fiberD) + ' um')

fig,ax = plt.subplots()
ax.plot(axonXYZ[:,1],axonXYZ[:,2])

ax.set(xlabel = 'Y location (m)', ylabel = 'Z location (m)',
	title = 'The Y and Z trajectory of axon 2038')


# -------------------------- Example 2 ---------------------------------
# Make a 3d scatter plot of only the nodes of axon 1800
fiberD, axonXYZ = loadAxon(1800,path2Files)

# Define the number of compartments between nodes
compBtwnNodes = 11

fig2 = plt.figure()
ax2 = plt.axes(projection="3d")

ax2.scatter3D(axonXYZ[::compBtwnNodes,0], axonXYZ[::compBtwnNodes,1], axonXYZ[::compBtwnNodes,2], c=axonXYZ[::compBtwnNodes,2], cmap='hsv');
ax2.set(xlabel = 'X location (m)', ylabel = 'Y location (m)', zlabel = 'Z location (m)',
	title = 'The XYZ trajectory of the nodes of axon 1800')


# -------------------------- Example 3 ---------------------------------
# Make a 3d plot of the first 25 axons

fig3 = plt.figure()
ax3 = plt.axes(projection="3d")

for axonNum in range(0,25):
	fiberD, axonXYZ = loadAxon(axonNum,path2Files)
	ax3.plot(axonXYZ[:,0],axonXYZ[:,1],axonXYZ[:,2])

ax3.set(xlabel = 'X location (m)', ylabel = 'Y location (m)', zlabel = 'Z location (m)',
	title = 'The XYZ trajectory of the first 25 axons')

plt.show()