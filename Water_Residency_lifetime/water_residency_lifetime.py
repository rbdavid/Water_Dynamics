#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
##!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

import sys
import numpy as np
from numpy.linalg import *
import MDAnalysis
from MDAnalysis.analysis.align import *
from distance_functions import *

zeros = np.zeros
dot_prod = np.dot
sqrt = np.sqrt
flush = sys.stdout.flush

# ----------------------------------------
# VARIABLE DECLARATION

config_file = sys.argv[1]

necessary_parameters = ['pdb_file','prmtop_file','traj_file','pocket_selection','wat_O_name','wat_H_name']
all_parameters = ['pdb_file','prmtop_file','traj_file','pocket_selection','wat_resname','pocket_radius','number_of_wats_filename','wat_res_nums_filename','center_of_geometry_filename','correlation_filename','long_lived_wat_filename','Wrapped','water_OH_bond_dist','summary_bool','summary_filename','exclude_waters_bool','exclude_waters_selection','VMD_vis_bool','VMD_vis_filename','VMD_step']

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
	print '%s' %(string)
	flush()

def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''

	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
	parameters['wat_resname'] = 'WAT'
	parameters['pocket_radius'] = 6.0
	parameters['number_of_wats_filename'] = 'num_wats_pocket.dat'
	parameters['wat_res_nums_filename'] = 'res_nums_wats.dat'
	parameters['center_of_geometry_filename'] = 'COG_pocket.xyz'
	parameters['correlation_filename'] = 'autocorrelation.dat'
	parameters['long_lived_wat_filename'] = 'long_lived_wats.vmd'
	parameters['Wrapped'] = True
	parameters['water_OH_bond_dist'] = 0.9572 		# assumes the use of TIP3P water model
	parameters['summary_bool'] = True
	parameters['summary_filename'] = 'water_diffusion_analysis.summary'
	parameters['exclude_waters_bool'] = False
	parameters['exclude_waters_selection'] = None
	parameters['VMD_vis_bool'] = False
	parameters['VMD_vis_filename'] = 'COG_vis_state.vmd'
	parameters['VMD_step'] = 100

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
			sys.exit()

def summary(filename):
	with open(filename,'w') as W:
		W.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
		W.write('To recreate this analysis, run this line:\n')
		for i in range(len(sys.argv)):
			W.write('%s ' %(sys.argv[i]))
		W.write('\n\nParameters used:\n')
		for i in all_parameters:
			W.write('%s = %s \n' %(i,parameters[i]))
		W.write('\n\n')

# ----------------------------------------
# MAIN:
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

# ----------------------------------------
# LOADING IN ANALYSIS UNIVERSE AND CREATING THE NECESSARY ATOM SELECTIONS
ffprint('Loading Analysis Universe.')
u = MDAnalysis.Universe(parameters['prmtop_file'],parameters['pdb_file'])
u_all = u.select_atoms('all')
wat = u.select_atoms('resname %s' %(parameters['wat_resname']))
u_pocket = u.select_atoms(parameters['pocket_selection'])

ffprint('Grabbing the Positions of the Binding Pocket to use as reference.')
u_all.translate(-u_pocket.center_of_geometry())
pocket_ref = u_pocket.positions

ffprint('Loading in the trajectory.')
u.load_new(parameters['traj_file'])

nSteps = len(u.trajectory)		# number of steps
nWats = wat.n_residues			# number of water residues
nRes0 = wat.residues[0].resid		# residue index is 0 indexed 

if nWats*3 != wat.n_atoms:
	ffprint('nWats*3 != wat.n_atoms. Unexpected number of water atoms. Possible selection error with water residue name.')
	sys.exit()

# ----------------------------------------
# MEMORY DECLARATION
oxygen_Coord = np.full((nSteps,nWats,3),np.nan)		# array holding all xyz data for water oxygens within the binding pocket; need to initialize this array for all water residues since any residue could find its way into the pocket; filling this array with nan values for subsequent easy replacement and testing.
OH_vector = np.full((nSteps,nWats,3),np.nan)		# array holding the O-H bond vector for waters within the binding pocket; need to initialize this array for all water residues since any residue could find its way into the pocket; filling this array with nan values for subsequent easy replacement and testing.
correlation_data = zeros((nSteps,4),dtype=np.float64)		# array holding the msd/bond autocorrelation data

# ----------------------------------------
# TRAJECTORY ANALYSIS 
if parameters['exclude_waters_bool']:
	with open(parameters['number_of_wats_filename'],'w') as X, open(parameters['wat_res_nums_filename'],'w') as Y, open(parameters['center_of_geometry_filename'],'w') as Z:
		ffprint('Beginning trajectory analysis')
		# Loop through trajectory
		for ts in u.trajectory:
			t = u_pocket.center_of_geometry()
			Z.write('1\n  generated by MDAnalysis and RBD\n X         %10.4f         %10.4f         %10.4f\n' %(t[0], t[1], t[2]))	#Writing an xyz trajectory of the center of geometry of the binding pocket; the COG particle is labeled as a dummy atom X
			
			u_all.translate(-t)	# Align to reference (moves COG of the pocket to origin)
	
			if not parameters['Wrapped']:
				dims = u.dimensions[:3]	# obtain dimension values to be used for wrapping atoms
				dims2 = dims/2.0
				for i in range(nWats):
					temp = wat.reidues[i].atom[0].position
					t = wrapping(temp,dims,dims2)
					wat.residues[i].translate(t)
	
			R, rmsd = rotation_matrix(u_pocket.positions,pocket_ref)	# Calculate the rotational matrix to align u to the ref, using the pocket selection as the reference selection
			u_all.rotate(R)
		
			pocket_waters = wat.select_atoms('byres point 0 0 0 %d' %(parameters['pocket_radius'])) # Atom selection for the waters within radius angstroms of the COG of the pocket; Assumes that the COG of the pocket is at 0,0,0 xyz coordinates (which it should be bc the translational motion of the pocket is removed...
			exclude_waters = pocket_waters.select_atoms('not %s'%(parameters['exclude_waters_selection']))

			nRes = pocket_waters.n_residues		# Calculate the number of waters within the pocket volume
			exclude_nRes = exclude_waters.n_residues	# Calculate the number of waters within the pocket volume that aren't excluded
			X.write('%d   %d\n' %(nRes,exclude_nRes))		# Outputting the number of water residues at timestep ts

			timestep = ts.frame
			for i in range(exclude_nRes):
				res = pocket_waters.residues[i]		#
				res_index = res.resid-nRes0			# grabbing the resid of the residue; needs to be zero-indexed for appropriate array assignment
				ox_pos = res.select_atoms('name %s' %(parameters['wat_O_name'])).positions[0]
				hy_pos = res.select_atoms('name %s' %(parameters['wat_H_name'])).positions[0]
				oxygen_Coord[timestep,res_index,:] = ox_pos
				OH_vector[timestep,res_index,:] = (ox_pos - hy_pos)/parameters['water_OH_bond_dist'] # NO NEED TO CALC THE MAGNITUDE OF THIS VECTOR BECAUSE I KNOW IF FROM THE PARAMETERS OF TIP3 (OR OTHER WATER MODEL) 
				Y.write('%d   ' %(res_index+nRes0))	# atom_num is zero-indexed; for vmd, need one-indexed values...
			Y.write('\n')
	ffprint('Done with saving coordinates of waters within the pocket, O-H bond vectors, writing COG traj, etc.\n Beginning msd calculations.')

else:
	with open(parameters['number_of_wats_filename'],'w') as X, open(parameters['wat_res_nums_filename'],'w') as Y, open(parameters['center_of_geometry_filename'],'w') as Z:
		ffprint('Beginning trajectory analysis')
		# Loop through trajectory
		for ts in u.trajectory:
			t = u_pocket.center_of_geometry()
			Z.write('1\n  generated by MDAnalysis and RBD\n X         %10.4f         %10.4f         %10.4f\n' %(t[0], t[1], t[2]))	#Writing an xyz trajectory of the center of geometry of the binding pocket; the COG particle is labeled as a dummy atom X
			
			u_all.translate(-t)	# Align to reference (moves COG of the pocket to origin)
	
			if not parameters['Wrapped']:
				dims = u.dimensions[:3]	# obtain dimension values to be used for wrapping atoms
				dims2 = dims/2.0
				for i in range(nWats):
					temp = wat.reidues[i].atom[0].position
					t = wrapping(temp,dims,dims2)
					wat.residues[i].translate(t)
	
			R, rmsd = rotation_matrix(u_pocket.positions,pocket_ref)	# Calculate the rotational matrix to align u to the ref, using the pocket selection as the reference selection
			u_all.rotate(R)
		
			pocket_waters = wat.select_atoms('byres point 0 0 0 %d' %(parameters['pocket_radius'])) # Atom selection for the waters within radius angstroms of the COG of the pocket; Assumes that the COG of the pocket is at 0,0,0 xyz coordinates (which it should be bc the translational motion of the pocket is removed...
		
			nRes = pocket_waters.n_residues		# Calculate the number of wates within the pocket volume
			X.write('%d\n' %(nRes))		# Outputting the number of water residues at timestep ts
	
			timestep = ts.frame
			for i in range(nRes):
				res = pocket_waters.residues[i]		#
				res_index = res.resid-nRes0			# grabbing the resid of the residue; needs to be zero-indexed for appropriate array assignment
				ox_pos = res.select_atoms('name %s' %(parameters['wat_O_name'])).positions[0]
				hy_pos = res.select_atoms('name %s' %(parameters['wat_H_name'])).positions[0]
				oxygen_Coord[timestep,res_index,:] = ox_pos
				OH_vector[timestep,res_index,:] = (ox_pos - hy_pos)/parameters['water_OH_bond_dist'] # NO NEED TO CALC THE MAGNITUDE OF THIS VECTOR BECAUSE I KNOW IF FROM THE PARAMETERS OF TIP3 (OR OTHER WATER MODEL) 
				Y.write('%d   ' %(res_index+nRes0))	# atom_num is zero-indexed; for vmd, need one-indexed values...
			Y.write('\n')
	
	ffprint('Done with saving coordinates of waters within the pocket, O-H bond vectors, writing COG traj, etc.\n Beginning msd calculations.')

# ------------------------------------------
# ANALYSIS OF TRAJECTORY DATA - MSD AND O-H BOND AUTOCORRELATION AND LONG-LIVED WATER RESIDUES

long_lived = set()
for i in range(nWats):		# Looping through all water residues.
	for j in range(nSteps):	# Looping through all timesteps for a single water.
		if oxygen_Coord[j,i,0] == oxygen_Coord[j,i,0]:	# boolean test to see if array object has a nan value or not; nan values will not equate and produce a FALSE;
			dt=1
			pos0 = oxygen_Coord[j,i,:]
			vec0 = OH_vector[j,i,:]
			while (j+dt)<nSteps and oxygen_Coord[j+dt,i,0] == oxygen_Coord[j+dt,i,0]:	# 
				if dt == 200 and i+nRes0+1 not in long_lived:	# the water molecule has resided in the pocket for 200 frames (or more) AND has not been added to the set already;
					long_lived.add(i+nRes0+1)	# saving the one-indexed residue index for long-lived water molecules 
				
				pos1 = oxygen_Coord[j+dt,i,:]
				vec1 = OH_vector[j+dt,i,:]

				dist, dist2 = euclid_dist(pos0,pos1)	# Calculates the MSD of the oxygen atoms in the water molecule
				scalar_product = dot_prod(vec0,vec1)

				correlation_data[dt,0]+=1		# count array element
				correlation_data[dt,1]+= dist2	# sum of MSD values
				correlation_data[dt,2]+= dist2**2	# sum of MSD^2 values
				correlation_data[dt,3]+= scalar_product
				dt+=1			# increment the dt value

ffprint('Finished with dist2 calculations. Beginning to average and print out msd values')

# ----------------------------------------
# OUTPUTTING DATA TO FILE

with open(parameters['correlation_filename'],'w') as W:
	for i in range(1,nSteps):
		if correlation_data[i,0]>1.0:
			correlation_data[i,1]/=correlation_data[i,0]	# Finish the average of the MSD value for the dt
			correlation_data[i,2]/=correlation_data[i,0]	# Finish the average of the MSD^2 value for the dt
			correlation_data[i,3]/=correlation_data[i,0]	# Finish the average of the Velocity*Velocity autocorrelation 
		W.write('%10.d   %10.d   %10.6f   %10.6f   %10.6f \n' %(i,correlation_data[i,0],correlation_data[i,1],correlation_data[i,2],correlation_data[i,3]))

ffprint('Writing the unique TIP3 oxygen numbers that are found to be within the binding pocket for longer than 200 frames.')
with open(parameters['long_lived_wat_filename'],'w') as W:
	ll_list = list(long_lived)
	W.write('resid ')
	for i in range(len(ll_list)):
		W. write('%d '%(ll_list[i]))
	W.write('\n')

# ----------------------------------------
# OUTPUTTING VMD VIS STATE FILE

if parameters['VMD_vis_bool']:
	with open(parameters['VMD_vis_filename'],'w') as W:
		W.write('# VMD Visualization State - COG of protein pocket within which water dynamics are analyzed; Written by RBD\n# Load into VMD using "vmd -e %s \n\nmol new %s\n mol addfile %s step %d\n mol modselect 0 0 "protein"\n mol modstyle 0 0 "NewCartoon"\n\n mol new %s step %d\n set sel [atomselect top "name X"]\n $sel set radius %d\nmol modselect 0 1 "name X"\n mol modstyle 0 1 "VDW"\n mol modmaterial 0 1 "Transparent"\n mol modcolor 0 1 "ColorID 23"\n\n'%(parameters['VMD_vis_filename'],parameters['prmtop_file'],parameters['traj_file'],parameters['VMD_step'],parameters['center_of_geometry_filename'],parameters['VMD_step'],parameters['pocket_radius']))

# ----------------------------------------
# OUTPUTTING SUMMARY OF THE ANALYSIS
if parameters['summary_bool']:
	summary(parameters['summary_filename'])

