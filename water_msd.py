#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
##!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

import numpy as np
from numpy.linalg import *
import MDAnalysis
from MDAnalysis.analysis.align import *
import sys
import os
from distance_functions import *
from pocket_residues import *

# ----------------------------------------
# VARIABLE DECLARATION

pdb_file = sys.argv[1]
prmtop_file = sys.argv[2]
traj_file = sys.argv[3]
system = sys.argv[4]

zeros = np.zeros
dot_prod = np.dot
sqrt = np.sqrt
flush = sys.stdout.flush

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
	print '%s' %(string)
	flush()

def summary(nSteps):
	sum_file = open('%s.diffusion.summary' %(system),'w')
	sum_file.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
	sum_file.write('To recreate this analysis, run this line:\n')
	for i in range(len(sys.argv)):
		sum_file.write('%s ' %(sys.argv[i]))
	sum_file.write('\n\n')
	sum_file.write('Reprinting the pocket_residues.py script for recreating this analysis:\n')
	sum_file.write("wat_resname = '%s'\n pocket_sel = '%s'\n radius = %f" %(wat_resname,pocket_sel,radius)) 
	sum_file.write('\n\n')
	sum_file.write('Output files are:\n')
	sum_file.write('	%s.nRes.dat 		--> holds the number of residues within the defined pocket\n' %(system))
	sum_file.write('	%s.waters.dat 		--> holds the atom numbers of the water oxygens that are within the defined pocket\n' %(system))
	sum_file.write('	%s.cog.xyz		--> A XYZ trajectory of the Center Of Geometry (COG) of the defined pocket; waters within radius of this point are considered within the defined pocket; Note you should always check that the radius is accurately describing all timesteps in your trajectory...\n' %(system))
	sum_file.write('	%s.msd.pocket.dat 	--> holds the dt, counts, <MSD>, and <MSD^2> data of waters within the defined pocket\n' %(system))
	sum_file.write('	%s.long_lived.dat	--> outputs a vmd atomselection to use to show the long lived waters in the un-truncated trajectories\n' %(system))
	sum_file.close()

# ----------------------------------------
# MAIN:

ffprint('Loading Reference Structure')
ref = MDAnalysis.Universe(prmtop_file,pdb_file)
ref_all = ref.select_atoms('all')
ref_pocket = ref.select_atoms(pocket_sel)
ref_all.translate(-ref_pocket.center_of_geometry())		
ref0 = ref_pocket.positions

ffprint('Loading Analysis Universe')
u = MDAnalysis.Universe(prmtop_file,traj_file)
u_all = u.select_atoms('all')
wat = u.select_atoms(wat_resname)
u_pocket = u.select_atoms(pocket_sel)

nSteps = len(u.trajectory)		# number of steps
nWats = wat.n_residues			# number of water residues
nAtoms = wat.n_atoms			# number of atoms in water selection...
nAtoms0 = wat.atoms[0].index		# atoms.index is 0 indexed, so leave as is...

if nWats*3 != nAtoms:
	ffprint('Something is fucked up. Selection issues. nWats*3 != nAtoms...')
	sys.exit()

# Memory Declaration
allCoord = zeros((nSteps,nAtoms,3),dtype=np.float64)	# array holding all xyz data for all atoms in the system
msd = zeros((nSteps,3),dtype=np.float64)		# array holding the msd data

# File Declaration
nRes_file = open('%s.nRes.dat' %(system),'w')	# File to print out the number of residues within the binding pocket
res_num = open('%s.waters.dat' %(system),'w')	# File to print out the atom numbers of the oxygen atom of TIP3 residues within the binding pocket selection
COG_file = open('%s.cog.xyz' %(system),'w')	# File to print a vmd xyz trajectory showing the motion of the COG of the binding pocket selection

ffprint('Beginning trajectory analysis')
# Loop through trajectory
for ts in u.trajectory:
	t = u_pocket.center_of_geometry()
	#Writing an xyz trajectory of the center of geometry of the binding pocket; the COG particle is labeled as a dummy atom X
	COG_file.write('1\n  generated by MDAnalysis and RBD\n X         %10.4f         %10.4f         %10.4f\n' %(t[0], t[1], t[2]))	
	# Align to reference (moves COG of the pocket to origin)
	u_all.translate(-t)

	# obtain dimension values to be used for wrapping atoms
	dims = u.dimensions[:3]
	dims2 = dims/2.0
	# Fix the wrapping issues
	for i in range(0,nAtoms,3):
		temp = wat.atoms[i].position
		t = wrapping(temp,dims,dims2)
		# translate the atoms of the residue using the translational matrix
		wat.atoms[i:i+3].translate(t)

	# Calculate the rotational matrix to align u to the ref, using the pocket selection as the reference selection
	R, rmsd = rotation_matrix(u_pocket.positions,ref0)
	u_all.rotate(R)

	x,y,z = u_pocket.center_of_geometry()
	pocket_waters = u.select_atoms('%s and byres point %s %s %s %d' %(wat_resname,x,y,z,radius)) # Atom selection for the waters within radius angstroms of the COG of the pocket	### CHECK THAT THIS ATOM SELECTION STILL WORKS
	
	nRes = pocket_waters.n_residues
	num_atoms = nRes*3
	nRes_file.write('%d\n' %(nRes))		# Outputting the number of water residues at timestep ts
	time = ts.frame-1
	for i in range(num_atoms):
		atom = pocket_waters.atoms[i]
		atom_num = atom.index
		zeroed_atom_num = atom_num-nAtoms0
		allCoord[time,zeroed_atom_num,:] = atom.position	# Saving xyz coordinates of pocket waters to the allCoord array   ### NEED TO CHECK THAT THIS INDEXING WORKS
		if i%3 == 0:
			res_num.write('%d   ' %(atom_num))
	res_num.write('\n')

nRes_file.close()
res_num.close()
COG_file.close()
ffprint('Done with saving coordinates of waters within the pocket, writing COG traj, etc...')

# Analyze Binary to determine which residues/timesteps are to be analyzed and perform dist2 analysis
ffprint('Beginning msd calculations (binary analysis and dist2 running average calcs)')
long_lived = set()
for i in range(nWats):
	temp = 3*i
	for j in range(nSteps):
		if allCoord[j,temp,0] != 0:
			dt=1
			pos0 = zeros((3,3),dtype=float)			### IS THIS LINE NEEDED ACTUALLY??? I NEED TO FIGURE THIS OUT ONCE AND FOR ALL...
			pos0 = allCoord[j,temp:temp+3,:]

			while (j+dt)<nSteps and allCoord[j+dt,temp,0] != 0:
				if dt > 200:
					long_lived.add(temp+nAtoms0)
				pos1 = zeros((3,3),dtype=float)		### IS THIS LINE NEEDED ACTUALLY??? I NEED TO FIGURE THIS OUT ONCE AND FOR ALL...
				pos1 = allCoord[j+dt,temp:temp+3,:]
				dist2=0.
				dist2 = MSD(pos0,pos1,3)
				msd[dt,0]+=1
				msd[dt,1]+= dist2
				msd[dt,2]+= dist2**2
				dt+=1

ffprint('Finished with dist2 calculations. Beginning to average and print out msd values')
msd_file = open('%s.msd.pocket.dat' %(system),'w')

for i in range(1,nSteps):
	if msd[i,0]>1.0:
		msd[i,1] /= msd[i,0]		# AVERAGE OF MSD
		msd[i,2] /= msd[i,0]		# AVERAGE OF MSD**2
	msd_file.write('%10.d    %10.d    %10.6f    %10.6f \n' %(i, msd[i,0], msd[i,1], msd[i,2]))
msd_file.close()

ffprint('Writing the unique TIP3 oxygen numbers that are found to be within the binding pocket for longer than 200 frames...')
nf = open('%s.long_lived.dat' %(system), 'w')
ll_list = list(long_lived)
nf.write('%s and same residue as index ' %(wat_resname))
for i in range(len(ll_list)):
	nf.write('%d ' %(ll_list[i]))
nf.write('\n')
nf.close()

summary(nSteps)

