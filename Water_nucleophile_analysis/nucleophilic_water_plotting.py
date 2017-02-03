#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:
#./nucleophilic_water_plotting.py config_file

# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

sqrt = np.sqrt
nullfmt = NullFormatter()
my_cmap = plt.cm.get_cmap('jet')
my_cmap.set_under('w',alpha=0)

config_file = sys.argv[1]

# ----------------------------------------
# FUNCTIONS

necessary_parameters = ['nucleophile_water_datafile','average_water_datafile','distance_cutoff','angle_cutoff','start','stop','step','system_descriptor']
all_parameters = ['nucleophile_water_datafile','average_water_datafile','distance_cutoff','angle_cutoff','distance_bin_width','angle_bin_width','delta_write','homemade_time_range_bool','homemade_time_range_list','plotting_bool','nucleophile_water_postanalysis_output','average_water_postanalysis_output']

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
nucleophilic_water_data = np.loadtxt(parameters['nucleophile_water_datafile'])

distance_minimum = np.min(nucleophilic_water_data[:,2])
distance_maximum = np.max(nucleophilic_water_data[:,2])
angle_minimum = np.min(nucleophilic_water_data[:,3])
angle_maximum = np.max(nucleophilic_water_data[:,3])
num_bins = [int((distance_maximum-distance_minimum)/distance_bin_width)+1,int((angle_maximum-angle_minimum)/angle_bin_width)+1]

with open(parameters['nucleophile_water_postanalysis_output'],'w') as W:
	W.write('Distance Cutoff: %f \nAngle Cutoff: %f\n\nTotal Dataset\n	Smallest P_gamma-O_water distance (nucleophilic): %f\n	Largest P_gamma-O_water distance: %f\n	Smallest O_water-P_gamma-O_beta angle: %f\n	Largest O_water-P_gamma-O_beta angle (nucleophilic): %f\n\nNumber of bins in the distance dimension: %d\nNumber of bins in the angle dimension: %d\n\n' %(distance_cutoff,angle_cutoff,distance_minimum,distance_maximum,angle_minimum,angle_maximum,num_bins[0],num_bins[1]))

	for i in range(len(frame_ranges)):
		first_index = np.argwhere(nucleophilic_water_data[:,0]==frame_ranges[i][3])[0]
		last_index = np.argwhere(nucleophilic_water_data[:,0]==frame_ranges[i][4])[-1]
		W.write('Analyzing %d to %d ns; corresponds to frames %d to %d, which have indices %d to %d\nCorresponds to %d steps\n'%(frame_ranges[i][1],frame_ranges[i][2],frame_ranges[i][3],frame_ranges[i][4],first_index,last_index,frame_ranges[i][0]))
		if parameters['plotting_bool']:
			counts,xedges,yedges,image = plt.hist2d(nucleophilic_water_data[first_index:last_index,2],nucleophilic_water_data[first_index:last_index,3],bins=num_bins,cmap=my_cmap,normed=True,vmin=0.00001,vmax=0.06)
			cb1 = plt.colorbar()
			cb1.set_label('Probability Density')
			plt.xlabel(r'P$_{\gamma}$ - O$_{water}$ Distance ($\AA$)',size=20)
			plt.ylabel(r'O$_{\beta}$-P$_{\gamma}$-O$_{water}$ Angle (Deg)',size=20)
			plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
			plt.title('Water Position in relation to the Nucleophilic Attack\n%d to %d ns - %s' %(frame_ranges[i][1],frame_ranges[i][2],parameters['system_descriptor']))
			plt.xlim(parameters['plot_x_range_tuple'])
			plt.ylim(parameters['plot_y_range_tuple'])
			plt.tight_layout()
			plt.savefig('%04d.%04d.nucleophilic_wat.hist2d.png'%(frame_ranges[i][1],frame_ranges[i][2]),transparent=True,dpi=160)
			plt.close()
	
		wat_count, frame_count = Nucleophilic_probability(distance_cutoff,nucleophilic_water_data[first_index:last_index,2],angle_cutoff,nucleophilic_water_data[first_index:last_index,3],nucleophilic_water_data[first_index:last_index,0])
		wat_prob = wat_count/float(frame_ranges[i][0])
		frame_prob = frame_count/float(frame_ranges[i][0])
		wat_prob_error = sqrt(wat_count)/float(frame_ranges[i][0]) 
		frame_prob_error = sqrt(frame_count)/float(frame_ranges[i][0]) 

		W.write('Total number of nucleophilic waters in this time range: %d\nProbability of waters being nucleophilic (includes frames with multiple nucleophilic waters): %f\nProbability error arrived at by sqrt rule of counting experiments: %f\nTotal number of frames with nucleophilic waters: %d\nProbability of a frame having a nucleophilc water: %f\nProbability error arrived at by sqrt rule of counting experiments: %f\n\n' %(wat_count,wat_prob,wat_prob_error,frame_count,frame_prob,frame_prob_error))

# ----------------------------------------
avg_water_data = np.loadtxt(parameters['average_water_datafile'])
with open(parameters['average_water_postanalysis_output'],'w') as W:
	for i in range(len(frame_ranges)):
		first_index = frame_ranges[i][3]
		last_index = frame_ranges[i][4]
		avg_waters = np.sum(avg_water_data[first_index:last_index])
		avg_waters /= frame_ranges[i][0]
		std_waters = np.std(avg_water_data[first_index:last_index])
		W.write('Time range: %d to %d ns; average number of waters in the ATP binding pocket: avg %f st dev %f \n' %(frame_ranges[i][0],frame_ranges[i][1],avg_waters,std_waters))
	
