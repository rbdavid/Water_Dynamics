#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:
#./nucleophilic_water_plotting.py config_file

# ----------------------------------------
# PREAMBLE:

import numpy as np
import numpy.ma as ma
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

sqrt = np.sqrt
zeros = np.zeros
nullfmt = NullFormatter()
prob_density_cmap = plt.cm.get_cmap('Blues')
prob_density_cmap.set_under('w',alpha=0)
#my_cmap.set_over('w',alpha=0)
free_energy_cmap = plt.cm.get_cmap('jet')

config_file = sys.argv[1]

k = 0.001987 # Kcal K^-1 mol^-1
T = 300. # K
kT = k*T
four_pi =4.0* np.pi

# ----------------------------------------
# FUNCTIONS

necessary_parameters = ['nucleophile_water_datafile','average_water_datafile','distance_cutoff','angle_cutoff','start','stop','step','system_descriptor']
all_parameters = ['nucleophile_water_datafile','average_water_datafile','distance_cutoff','angle_cutoff','distance_bin_width','angle_bin_width','delta_write','homemade_time_range_bool','homemade_time_range_list','plotting_bool','Free_energy_bool','nucleophile_water_postanalysis_output','average_water_postanalysis_output']

def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''

	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
	parameters['distance_bin_width'] = float(0.05)		# Units: Angstroms
	parameters['angle_bin_width'] = float(1.50)		# Units: Degrees
	parameters['delta_write'] = float(0.002)		# Units: ns
	parameters['homemade_time_range_bool'] = False
	parameters['homemade_time_range_list'] = None
	parameters['plotting_bool'] = False
	parameters['Free_energy_bool'] = False
	parameters['nucleophile_water_postanalysis_output'] = 'nucleophilic_waters.dat'
	parameters['average_water_postanalysis_output'] = 'average_waters.dat'
	parameters['plot_x_range_tuple'] = (2.50,8.00)
	parameters['plot_y_range_tuple'] = (25.0,180.0)

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)

def Nucleophilic_probability(dist_cutoff,dist_array,ang_cutoff,ang_array,timestep_array):
	wat_count = 0
	frame_count = 0
	old_frame_number = np.nan 
	nWats = len(dist_array)
	if nWats != len(ang_array):
		print 'Number of elements in the distance and angle arrays read into the Nucleophilic_probability function do not match; This should not happen.'
		sys.exit()
	else:
		for i in range(nWats):
			current_frame_number = timestep_array[i]
			if dist_array[i]<dist_cutoff and ang_array[i]>ang_cutoff:	# DETERMINE IF THE WATER IS IN THE NUCLEOPHILIC POSITION
				wat_count += 1						# IF IT IS, ADD TO THE COUNTER
				if current_frame_number == old_frame_number:		# A NUCLEOPHILIC WATER HAS ALREADY BEEN FOUND AT THIS FRAME NUMBER 
					print 'multiple nucl. waters in frame:', timestep_array[i]
				else:							# A NUCLEOPHILIC WATER HAS NOT BEEN FOUND AT THIS FRAME NUMBER
					frame_count += 1
					old_frame_number = current_frame_number

	return wat_count, frame_count	# RETURNING TWO VALUES; THE TOTAL NUMBER OF NUCLEOPHILIC WATERS FOUND IN THE POCKET; THE TOTAL NUMBER OF FRAMES WITH A NUCLEOPHILIC WATER 

# ----------------------------------------
# MAIN
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

distance_cutoff = float(parameters['distance_cutoff'])
angle_cutoff = float(parameters['angle_cutoff'])
start = int(parameters['start'])
stop = int(parameters['stop'])
step = int(parameters['step'])
delta_write = float(parameters['delta_write'])
distance_bin_width = float(parameters['distance_bin_width'])
angle_bin_width = float(parameters['angle_bin_width'])

# PREP TIME RANGES TO ANALYZE THE NUCLEOPHILIC/AVERAGE WATER DATA
frame_ranges = []
for i in range(start,stop,step):	# i has units of ns
	j = i + step 			# j has units of ns
	first_frame = int(i/delta_write)
	last_frame = int(j/delta_write)-1
	nSteps = len(range(first_frame,last_frame))+1
	frame_ranges.append([nSteps,i,j,first_frame,last_frame])

if parameters['homemade_time_range_bool']:
	for i in range(len(parameters['homemade_time_range_list'])): # assumes the time_range_list is a list of lists....
		first_frame = int(parameters['homemade_time_range_list'][i][0]/delta_write)
		last_frame = int(parameters['homemade_time_range_list'][i][1]/delta_write)-1
		nSteps = len(range(first_frame,last_frame))+1
		frame_ranges.append([nSteps,parameters['homemade_time_range_list'][i][0],parameters['homemade_time_range_list'][i][1],first_frame,last_frame])

# ----------------------------------------
# LOAD THE NUCLEOPHILIC WATER DATAFILE INTO MEMORY
nucleophilic_water_data = np.loadtxt(parameters['nucleophile_water_datafile'])
# DETERMINE THE NUMBER OF BINS TO BE USED IN PLOTTING THE NUCLEOPHILIC DATA
distance_minimum = np.min(nucleophilic_water_data[:,2])
distance_maximum = np.max(nucleophilic_water_data[:,2])
angle_minimum = np.min(nucleophilic_water_data[:,3])
angle_maximum = np.max(nucleophilic_water_data[:,3])

num_bins = [int((distance_maximum-distance_minimum)/distance_bin_width)+1,int((angle_maximum-angle_minimum)/angle_bin_width)+1]

delta_dist = (distance_maximum - distance_minimum)/num_bins[0]
delta_angle = (angle_maximum - angle_minimum)/num_bins[1]

dist_edges = zeros(num_bins[0]+1,dtype=np.float64)
angle_edges = zeros(num_bins[1]+1,dtype=np.float64)
for i in range(num_bins[0]+1):
	dist_edges[i] = distance_minimum + i*delta_dist
for i in range(num_bins[1]+1):
	angle_edges[i] = angle_minimum + i*delta_angle

# ----------------------------------------
# BEGIN DATA ANALYSIS OF THE NUCLEOPHILIC WATER DATA
with open(parameters['nucleophile_water_postanalysis_output'],'w') as W:
	W.write('Distance Cutoff: %f \nAngle Cutoff: %f\n\nTotal Dataset\n	Smallest P_gamma-O_water distance (nucleophilic): %f\n	Largest P_gamma-O_water distance: %f\n	Smallest O_water-P_gamma-O_beta,gamma angle: %f\n	Largest O_water-P_gamma-O_beta,gamma angle (nucleophilic): %f\n\nNumber of bins in the distance dimension: %d\nNumber of bins in the angle dimension: %d\n Distance bin dimension: %f; Angle bin dimension: %f\n\n' %(distance_cutoff,angle_cutoff,distance_minimum,distance_maximum,angle_minimum,angle_maximum,num_bins[0],num_bins[1],delta_dist,delta_angle))

	for i in frame_ranges:
		first_index = np.argwhere(nucleophilic_water_data[:,0]==i[3])[0][0]	# SINCE THE DATA IS NOT ORGANIZED BY TIMESTEP, I NEED TO DETERMINE THE CORRECT INDEX FOR THE FIRST FRAME OF THE TIME RANGE
		last_index = np.argwhere(nucleophilic_water_data[:,0]==i[4])[-1][0]	# SIMILAR TO ABOVE; DETERMINE THE INDEX OF THE LAST FRAME IN THE TIME RANGE
		W.write('Analyzing %d to %d ns; corresponds to frames %d to %d, which have indices %d to %d\nCorresponds to %d steps\n'%(i[1],i[2],i[3],i[4],first_index,last_index,i[0]))

		# ----------------------------------------
		# PLOT THE PROB DENSITY OF THE NUCLEOPHILIC WATER DATA; NOTE: HIST2D FUNCTION CALCULATES THE PROB DENSITY AS NORMALIZED BY THE TOTAL NUMBER OF DATA POINTS (INSTEAD OF NSTEPS) * DX * DY
		if parameters['plotting_bool']:
			counts,xedges,yedges,image = plt.hist2d(nucleophilic_water_data[first_index:last_index,2],nucleophilic_water_data[first_index:last_index,3],bins=num_bins,cmap=prob_density_cmap,normed=True,vmin=0.001,vmax=0.020)
			cb1 = plt.colorbar(extend='max')		# extend='max'
			cb1.set_label('Probability Density')
			plt.xlabel(r'P$_{\gamma}$ - O$_{water}$ Distance ($\AA$)',size=20)
			plt.ylabel(r'O$_{\beta,\gamma}$-P$_{\gamma}$-O$_{water}$ Angle (Deg)',size=20)
			plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
			plt.xlim(parameters['plot_x_range_tuple'])
			plt.ylim(parameters['plot_y_range_tuple'])
			plt.tight_layout()
			plt.savefig('%04d.%04d.nucleophilic_wat.hist2d.png'%(i[1],i[2]),transparent=True,dpi=300)
			plt.close()

		# ----------------------------------------
		# CALCULATE,WRITE TO OUTPUT, AND PLOT THE RELATIVE FREE ENERGY OF THE NUCLEOPHILIC WATER DATA; NORMALIZED FOR SPHERICAL VOLUME
		if parameters['Free_energy_bool']:
			# DEFINE CONSTANTS
			nValues = i[0]
			fe_divisor = nValues*delta_dist*delta_angle*four_pi
			# ARRAY DECLARATION
			counts = np.zeros(tuple(num_bins),dtype=np.float64)
			fe_counts = np.zeros(tuple(num_bins),dtype=np.float64)
			# LOOP THROUGH AND HISTOGRAM DATA
			for a in nucleophilic_water_data[first_index:last_index]:
				distance_index = int((a[2] - distance_minimum)/delta_dist)
				angle_index = int((a[3] - angle_minimum)/delta_angle)
				
				if distance_index not in range(num_bins[0]+1) or angle_index not in range(num_bins[1]+1):
					print 'A specific water is trying to be binned in a bad bin. Here is the data that is bad:'
					print i
					sys.exit()
				if distance_index == num_bins[0]: 
					distance_index = -1 
				if angle_index == num_bins[1]: 
					angle_index = -1 
				else:
					counts[distance_index][angle_index] += 1
					fe_counts[distance_index][angle_index] += 1/(fe_divisor*a[2]**2)	# nSteps*delta_dist*delta_angle normalizes the probability density; 4.0*pi*a[2]**2 volume corrects the spherical volume these waters were found in. 
			# FROM THE HISTOGRAMMED ARRAYS, CALCULATE THE FREE ENERGY SURFACE
			for j in range(num_bins[0]):
				for k in range(num_bins[1]):
						fe_counts[j][k] = -kT*np.log(fe_counts[j][k])
			fe_counts -= np.ndarray.min(fe_counts)		# Make lowest energy bin the zero point reference
			
			# PCOLORMESH ERRORS WHEN THE C MATRIX HAS NANS OR INFS; CREATE A MASKED ARRAY SO THAT PCOLORMESH WILL IGNORE INFS IN THE FE_COUNTS ARRAY
			masked_fe_counts = ma.masked_where(np.isinf(fe_counts),fe_counts)
			# PLOT FREE ENERGY SURFACE
			free_plot = plt.pcolormesh(dist_edges,angle_edges,masked_fe_counts.T,cmap=free_energy_cmap)
			cb1 = plt.colorbar(extend='max')		# extend='max'
			cb1.set_label(r'Relative Free Energy (kcal mol$^{-1}$)')
			plt.xlabel(r'P$_{\gamma}$ - O$_{water}$ Distance ($\AA$)',size=20)
			plt.ylabel(r'O$_{\beta,\gamma}$-P$_{\gamma}$-O$_{water}$ Angle (Deg)',size=20)
			plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
			plt.xlim(parameters['plot_x_range_tuple'])
			plt.ylim(parameters['plot_y_range_tuple'])
			plt.tight_layout()
			plt.savefig('%04d.%04d.nucleophilic_wat.free_energy.png'%(i[1],i[2]),transparent=True,dpi=300)
			plt.close()

			with open('%04d.%04d.nucleophilic_wat.free_energy.dat'%(i[1],i[2]),'w') as Y:
				np.savetxt(Y,fe_counts)
		
		# ----------------------------------------
		# CALCULATE THE PROBABILITY OF OBSERVING A NUCLEOPHILIC WATER IN THE TIME RANGE; OUTPUT TO FILE
		wat_count, frame_count = Nucleophilic_probability(distance_cutoff,nucleophilic_water_data[first_index:last_index,2],angle_cutoff,nucleophilic_water_data[first_index:last_index,3],nucleophilic_water_data[first_index:last_index,0])
		wat_prob = wat_count/float(i[0])
		frame_prob = frame_count/float(i[0])
		wat_prob_error = sqrt(wat_count)/float(i[0]) 
		frame_prob_error = sqrt(frame_count)/float(i[0]) 

		W.write('Total number of nucleophilic waters in this time range: %d\nProbability of waters being nucleophilic (includes frames with multiple nucleophilic waters): %f\nProbability error arrived at by sqrt rule of counting experiments: %f\nTotal number of frames with nucleophilic waters: %d\nProbability of a frame having a nucleophilc water: %f\nProbability error arrived at by sqrt rule of counting experiments: %f\n\n' %(wat_count,wat_prob,wat_prob_error,frame_count,frame_prob,frame_prob_error))

# ----------------------------------------
# CALCULATE THE AVERAGE (AND STD DEV) NUMBER OF WATERS WITHIN THE BINDING POCKET
avg_water_data = np.loadtxt(parameters['average_water_datafile'])
with open(parameters['average_water_postanalysis_output'],'w') as W:
	for i in frame_ranges:
		first_index = i[3]
		last_index = i[4]
		avg_waters = np.sum(avg_water_data[first_index:last_index])
		avg_waters /= i[0]
		std_waters = np.std(avg_water_data[first_index:last_index])
		W.write('Time range: %d to %d ns; average number of waters in the ATP binding pocket: avg %f st dev %f \n' %(i[0],i[1],avg_waters,std_waters))
	
