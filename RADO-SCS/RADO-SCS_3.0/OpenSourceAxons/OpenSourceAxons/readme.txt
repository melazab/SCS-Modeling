RADO-SCS axon trajectories and diameters


Author: Hans Zander, PhD candidate
Institution: University of Michigan, Ann Arbor
email: hzander@umich.edu



This folder should contain the following:


axonDiameters.npy	-	A file containing the diameters of all 2039 axons for the fibers used in the RADO-SCS paper. All diameters are in um. This is an .npy file, and can be loaded in python with numpy.load(filePath).


axon_0.npy
axon_1.npy
...
axon_2038.npy		-	These files contain the x (column 1), y (col2) and z (col3) center locations of each compartment of each axon in meters. The diameter of the axon is in the axonDiameters.npy file (axon_0.npy -> axonDiameter[0] in the axonDiameters file).
	The center coordinate locations are sequential, and have the following structure:

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

	McIntyre, Cameron C., Andrew G. Richardson, and Warren M. Grill. "Modeling the excitability of mammalian nerve fibers: influence of afterpotentials on the recovery cycle." Journal of neurophysiology 87.2 (2002): 995-1006.
	

example.py 			- A python file with several examples of getting the diameter and the trajectories of various axons, as well as plotting them.