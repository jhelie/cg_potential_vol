#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.4"
parser = argparse.ArgumentParser(prog = 'cg_potential_vol', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
***********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@gmail.com)
git: https://github.com/jhelie/cg_potential_vol
DOI: 
***********************************************

DESCRIPTION
 
 This script calculates the electrosatic potential along z in CG systems - i.e. taking
 into account the force shift and cutoff.

 The 3D electrostatic potential is saved in an OpenDx file which can then be processed
 by the 'dx_plot' utility to produce the desired 2D or 1D graphs.

 A file containing the charged particles can also be supplied to calculate the density
 of charges.
 
 The script follows a 2-step process:
  1. calculate charge density along z for each slice
  2. sum the corresponding shifted potentials over all the slices

 charges identification
 ----------------------
   You can specify which particles to take into account for the calculation of the total
   charge density by supplying a file via the --charges option. Each line of this file
   should follow the format (without quotation marks):
    -> 'label,value,MDAnalysis selection string'

   Note that the MDAnalysis selection string should not contain any commas.

   The absolute charge for each group will be plotted on the charge density profile. The
   group colour must be specified for each charge.
  
   By default charges are defined using the Martini 2.1 settings:
    -> Na+,1,name NA+
    -> CL-,-1,name CL-
    -> PO4,-1,name PO4
    -> NC3,1,name NC3
    -> NH3,1,name NH3
    -> Lys,1,resname LYS and name SC2 (position: 3/3)
    -> Arg,1,resname ARG and name SC2 (position: 3/3)
    -> Asp,-1,resname ASP and name SC1 (position: 2/2)
    -> Glu,-1,resname GLU and name SC1 (position: 2/2)
  
   Another default set of charges can be used by specifiying --charges 2.2P :
    -> Na+,1,name NA+
    -> CL-,-1,name CL-
    -> WP,0.46,name WP
    -> WM,-0.46,name WM
    -> PO4,-1,name PO4
    -> NC3,1,name NC3
    -> NH3,1,name NH3
    -> Lys,1,resname LYS and name SC2 (position: 4/4)
    -> Arg,1,resname ARG and name SC2 (position: 4/4)
    -> Asp,-1,resname ASP and name SC1 (position: 3/3)
    -> Glu,-1,resname GLU and name SC1 (position: 3/3)
    -> Asn_p,0.46,resname ASN and name SCP (position: 3/4)
    -> Asn_n,-0.46,resname ASN and name SCN (position: 4/4)
    -> Gln_p,0.46,resname GLN and name SCP (position: 3/4)
    -> Gln_n,-0.46,resname GLN and name SCN (position: 4/4)
    -> Thr_p,0.31,resname THR and name SCP (position: 3/4)
    -> Thr_n,-0.31,resname THR and name SCN (position: 4/4)
    -> Ser_p,0.4,resname SER and name SCP (position: 3/4)
    -> Ser_n,-0.4,resname SER and name SCN (position: 4/4)

   By default proteins termini are considered uncapped and +1 and -1 are added to the
   first and last backbone beads ("name BB") respectively. If this is not what you want
   just use the option --capped.

   Note that for systems containing different protein species automatic addition of
   charges at the termini has not yet been implemented: you should use --capped and
   provide a file via --charges specifying all charged particles in the system, including
   termini beads.

 electrostatic potential
 -----------------------
   The shifted electrostatic potential is calculated as per Gromacs manual for each point charge.
   In particular when rs = 0, the potential created by a charge q at a distance r is:
    -> phi(r) = q/(4*pi*e0*er)*(1/r - 5/(3*rc) + 5*r^3/(3*rc^4) - r^4/rc^5)


REQUIREMENTS

 The following python modules are needed :
  - MDAnalysis
  - matplotlib
  - numpy
  - scipy
  - networkX (if option --algorithm is set to 'min' or 'cog')
  - sklearn (if option --algorithm is set to 'density')


NOTES

 1. The density is calculated with respect to the z axis, not the bilayer normal. So the
    more your system deforms the noiser and the less meaningful the 1D and 2D profiles get.

 
USAGE

Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro] (required)
-x			: trajectory file [.xtc]
-o			: name of output folder
-b			: beginning time (ns) (the bilayer must exist by then!)
-e			: ending time (ns)	
-t 		1	: process every t-frames
--charges	2.1	: definition of charged particles, see 'DESCRIPTION' 
--capped		: assumes protein termini are capped
--positions		: use bead positions in residues instead of bead names to select charged protein beads
--protein		: only take into account charges due to protein(s)
 
Grid options
-----------------------------------------------------
--sz		200 	: number of slizes along z
--sx		200 	: number of slizes along x and y

Electrostatic potential
-----------------------------------------------------
--er		15	: dielectric constant
--rs		0 	: distance from charge where the electrostatic starts to be shifted (Angstrom) [TO DO]
--rc		12	: distance from charge where the electrostatic should reach 0 (Angstrom)

Other options
-----------------------------------------------------
--debug		: for development
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-x', nargs=1, dest='xtcfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-b', nargs=1, dest='t_start', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-e', nargs=1, dest='t_end', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-t', nargs=1, dest='frames_dt', default=[1], type=int, help=argparse.SUPPRESS)
parser.add_argument('--charges', nargs=1, dest='chargesfilename', default=['2.1'], help=argparse.SUPPRESS)
parser.add_argument('--capped', dest='capped', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('--positions', dest='positions', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('--protein', dest='protein', action='store_true', help=argparse.SUPPRESS)

#grid options
parser.add_argument('--sx', nargs=1, dest='sx', default=[200], type=int, help=argparse.SUPPRESS)
parser.add_argument('--sz', nargs=1, dest='sz', default=[200], type=int, help=argparse.SUPPRESS)

#electrostatic potential
parser.add_argument('--er', nargs=1, dest='er', default=[15], type=float, help=argparse.SUPPRESS)
parser.add_argument('--rs', nargs=1, dest='rs', default=[0], type=float, help=argparse.SUPPRESS)
parser.add_argument('--rc', nargs=1, dest='rc', default=[12], type=float, help=argparse.SUPPRESS)

#lipids identification options
parser.add_argument('--beads', nargs=1, dest='beadsfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--flipflops', nargs=1, dest='selection_file_ff', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--leaflets', nargs=1, dest='cutoff_leaflet', default=['optimise'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--debug', dest='debug', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
#data options
args.grofilename = args.grofilename[0]
args.xtcfilename = args.xtcfilename[0]
args.output_folder = args.output_folder[0]
args.t_start = args.t_start[0]
args.t_end = args.t_end[0]
args.frames_dt = args.frames_dt[0]
#density profile options
args.chargesfilename = args.chargesfilename[0]
args.sx = args.sx[0]
args.sz = args.sz[0]
##electrostatic potential
args.er = args.er[0]
args.rs = args.rs[0]/float(10)				#conversion to nm for unit consistency
args.rc = args.rc[0]/float(10)				#conversion to nm for unit consistency

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================
#generic science modules
try:
	import math
except:
	print "Error: you need to install the maths module."
	sys.exit(1)
try:
	import numpy as np
except:
	print "Error: you need to install the numpy module."
	sys.exit(1)
try:
	import scipy as sp
	import scipy.stats
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#MDAnalysis module
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.analysis
	import MDAnalysis.analysis.leaflet
	import MDAnalysis.analysis.distances
	#set MDAnalysis to use periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

#=========================================================================================
# sanity check
#=========================================================================================
if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)
if args.chargesfilename not in ["2.1","2.2P"] and not os.path.isfile(args.chargesfilename):
	print "Error: file " + str(args.chargesfilename) + " not found."
	sys.exit(1)
if args.chargesfilename not in ["2.1","2.2P"] and args.positions:
	print "Error: --position cannot be used with user specified --charges"
	sys.exit(1)
if args.t_end != -1 and args.t_end < args.t_start:
	print "Error: the starting time (" + str(args.t_start) + "ns) for analysis is later than the ending time (" + str(args.t_end) + "ns)."
	sys.exit(1)

if args.xtcfilename == "no":
	if '-t' in sys.argv:
		print "Error: -t option specified but no xtc file specified."
		sys.exit(1)
	elif '-b' in sys.argv:
		print "Error: -b option specified but no xtc file specified."
		sys.exit(1)
	elif '-e' in sys.argv:
		print "Error: -e option specified but no xtc file specified."
		sys.exit(1)
elif not os.path.isfile(args.xtcfilename):
	print "Error: file " + str(args.xtcfilename) + " not found."
	sys.exit(1)

#=========================================================================================
# create folders and log file
#=========================================================================================
if args.output_folder == "no":
	if args.xtcfilename == "no":
		args.output_folder = "cg_potential_vol_" + args.grofilename[:-4]
	else:
		args.output_folder = "cg_potential_vol_" + args.xtcfilename[:-4]
if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + " already exists, choose a different output name via -o."
	sys.exit(1)
else:
	os.mkdir(args.output_folder)
	filename_log = os.getcwd() + '/' + str(args.output_folder) + '/cg_potential_vol.log'
	output_log = open(filename_log, 'w')		
	output_log.write("[cg_potential_vol v" + str(version_nb) + "]\n")
	output_log.write("\nThis folder and its content were created using the following command:\n\n")
	tmp_log="python cg_potential_vol.py"
	for c in sys.argv[1:]:
		tmp_log += " " + c
	output_log.write(tmp_log + "\n")
	output_log.close()
	
	#copy input files
	if args.chargesfilename not in ["no","2.1","2.2P"]:
		shutil.copy2(args.chargesfilename,args.output_folder + "/")

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

global nb_x
global nb_y
global nb_z	
global delta_x
global delta_y
global delta_z	
global coords_x
global coords_z
global z_middle
global upper_avg
global lower_avg
global f_factor
global pot_conv
global potential
global potential_1D
global potential_2D
global charge_density
global charge_density_1D
global charge_density_2D

f_factor = 138.935485/float(args.er)
pot_conv = 0.010364272			#to convert from kJ.mol-1.e-1 to V
potential = np.zeros((args.sx, args.sx, args.sz))
potential_1D = np.zeros(args.sz)
potential_2D = np.zeros((args.sx,args.sz))
charge_density = np.zeros((args.sx, args.sx, args.sz))
charge_density_1D = np.zeros(args.sz)
charge_density_2D = np.zeros((args.sx,args.sz))

#=========================================================================================
# data loading
#=========================================================================================

def set_charges():
			
	global charges_groups
	global charged_residues
	charges_groups = {}
	
	#use martini 2.1
	#---------------
	if args.chargesfilename == "2.1":		
		if not args.protein:
			#solvent
			charges_groups["Na+"] = {}
			charges_groups["Na+"]["value"] = 1
			charges_groups["Na+"]["sele_string"] = "name NA+"
			charges_groups["CL-"] = {}
			charges_groups["CL-"]["value"] = -1
			charges_groups["CL-"]["sele_string"] = "name CL-"
	
			#lipids
			charges_groups["PO4"] = {}
			charges_groups["PO4"]["value"] = -1
			charges_groups["PO4"]["sele_string"] = "name PO4"
			charges_groups["NH3"] = {}
			charges_groups["NH3"]["value"] = 1
			charges_groups["NH3"]["sele_string"] = "name NH3"
			charges_groups["NC3"] = {}
			charges_groups["NC3"]["value"] = 1
			charges_groups["NC3"]["sele_string"] = "name NC3"
	
		#protein: use bead names
		charged_residues = ["LYS","ARG","ASP","GLU"]
		if not args.positions:
			charges_groups["LYS"] = {}
			charges_groups["LYS"]["value"] = 1
			charges_groups["LYS"]["sele_string"] = "resname LYS and name SC2"
			charges_groups["ARG"] = {}
			charges_groups["ARG"]["value"] = 1
			charges_groups["ARG"]["sele_string"] = "resname ARG and name SC2"
			charges_groups["ASP"] = {}
			charges_groups["ASP"]["value"] = -1
			charges_groups["ASP"]["sele_string"] = "resname ASP and name SC1"
			charges_groups["GLU"] = {}
			charges_groups["GLU"]["value"] = -1
			charges_groups["GLU"]["sele_string"] = "resname GLU and name SC1"
		#protein: use bead position
		else:
			charges_groups["LYS"] = {}
			charges_groups["LYS"]["value"] = 1
			charges_groups["LYS"]["sele_string"] = "resname LYS"
			charges_groups["LYS"]["position"] = 2
			charges_groups["ARG"] = {}
			charges_groups["ARG"]["value"] = 1
			charges_groups["ARG"]["sele_string"] = "resname ARG"
			charges_groups["ARG"]["position"] = 2
			charges_groups["ASP"] = {}
			charges_groups["ASP"]["value"] = -1
			charges_groups["ASP"]["sele_string"] = "resname ASP"
			charges_groups["ASP"]["position"] = 1
			charges_groups["GLU"] = {}
			charges_groups["GLU"]["value"] = -1
			charges_groups["GLU"]["sele_string"] = "resname GLU"
			charges_groups["GLU"]["position"] = 1
			
 
	#use martini 2.2P
	#----------------
	elif args.chargesfilename == "2.2P":
		if not args.protein:
			#solvent
			charges_groups["Na+"] = {}
			charges_groups["Na+"]["value"] = 1
			charges_groups["Na+"]["sele_string"] = "name NA+"
			charges_groups["CL-"] = {}
			charges_groups["CL-"]["value"] = -1
			charges_groups["CL-"]["sele_string"] = "name CL-"
			charges_groups["WP"] = {}
			charges_groups["WP"]["value"] = 0.46
			charges_groups["WP"]["sele_string"] = "name WP"
			charges_groups["WM"] = {}
			charges_groups["WM"]["value"] = -0.46
			charges_groups["WM"]["sele_string"] = "name WM"
	
			#lipids
			charges_groups["PO4"] = {}
			charges_groups["PO4"]["value"] = -1
			charges_groups["PO4"]["sele_string"] = "name PO4"
			charges_groups["NH3"] = {}
			charges_groups["NH3"]["value"] = 1
			charges_groups["NH3"]["sele_string"] = "name NH3"
			charges_groups["NC3"] = {}
			charges_groups["NC3"]["value"] = 1
			charges_groups["NC3"]["sele_string"] = "name NC3"

		#protein
		charged_residues = ["LYS","ARG","ASP","GLU","ASN_p","GLN_p","THR_p","SER_p","ASN_n","GLN_n","THR_n","SER_n"]
		if not args.positions:
			charges_groups["LYS"] = {}
			charges_groups["LYS"]["value"] = 1
			charges_groups["LYS"]["sele_string"] = "resname LYS and name SCP"
			charges_groups["ARG"] = {}
			charges_groups["ARG"]["value"] = 1
			charges_groups["ARG"]["sele_string"] = "resname ARG and name SCP"
			charges_groups["ASP"] = {}
			charges_groups["ASP"]["value"] = -1
			charges_groups["ASP"]["sele_string"] = "resname ASP and name SCN"
			charges_groups["GLU"] = {}
			charges_groups["GLU"]["value"] = -1
			charges_groups["GLU"]["sele_string"] = "resname GLU and name SCN"
			charges_groups["ASN_p"] = {}
			charges_groups["ASN_p"]["value"] = 0.46
			charges_groups["ASN_p"]["sele_string"] = "resname ASN and name SCP"
			charges_groups["GLN_p"] = {}
			charges_groups["GLN_p"]["value"] = 0.46
			charges_groups["GLN_p"]["sele_string"] = "resname GLN and name SCP"
			charges_groups["THR_p"] = {}
			charges_groups["THR_p"]["value"] = 0.31
			charges_groups["THR_p"]["sele_string"] = "resname THR and name SCP"
			charges_groups["SER_p"] = {}
			charges_groups["SER_p"]["value"] = 0.4
			charges_groups["SER_p"]["sele_string"] = "resname SER and name SCP"
			charges_groups["ASN_n"] = {}
			charges_groups["ASN_n"]["value"] = -0.46
			charges_groups["ASN_n"]["sele_string"] = "resname ASN and name SCN"
			charges_groups["GLN_n"] = {}
			charges_groups["GLN_n"]["value"] = -0.46
			charges_groups["GLN_n"]["sele_string"] = "resname GLN and name SCN"
			charges_groups["THR_n"] = {}
			charges_groups["THR_n"]["value"] = -0.31
			charges_groups["THR_n"]["sele_string"] = "resname THR and name SCN"
			charges_groups["SER_n"] = {}
			charges_groups["SER_n"]["value"] = -0.4
			charges_groups["SER_n"]["sele_string"] = "resname SER and name SCN"
		else:
			charges_groups["LYS"] = {}
			charges_groups["LYS"]["value"] = 1
			charges_groups["LYS"]["sele_string"] = "resname LYS"
			charges_groups["LYS"]["position"] = 3
			charges_groups["ARG"] = {}
			charges_groups["ARG"]["value"] = 1
			charges_groups["ARG"]["sele_string"] = "resname ARG"
			charges_groups["ARG"]["position"] = 3
			charges_groups["ASP"] = {}
			charges_groups["ASP"]["value"] = -1
			charges_groups["ASP"]["sele_string"] = "resname ASP"
			charges_groups["ASP"]["position"] = 2
			charges_groups["GLU"] = {}
			charges_groups["GLU"]["value"] = -1
			charges_groups["GLU"]["sele_string"] = "resname GLU"
			charges_groups["GLU"]["position"] = 2
			
			charges_groups["ASN_p"] = {}
			charges_groups["ASN_p"]["value"] = 0.46
			charges_groups["ASN_p"]["sele_string"] = "resname ASN"
			charges_groups["ASN_p"]["position"] = 2
			charges_groups["GLN_p"] = {}
			charges_groups["GLN_p"]["value"] = 0.46
			charges_groups["GLN_p"]["sele_string"] = "resname GLN"
			charges_groups["GLN_p"]["position"] = 2
			charges_groups["THR_p"] = {}
			charges_groups["THR_p"]["value"] = 0.31
			charges_groups["THR_p"]["sele_string"] = "resname THR"
			charges_groups["THR_p"]["position"] = 2
			charges_groups["SER_p"] = {}
			charges_groups["SER_p"]["value"] = 0.4
			charges_groups["SER_p"]["sele_string"] = "resname SER"
			charges_groups["SER_p"]["position"] = 2
			
			charges_groups["ASN_n"] = {}
			charges_groups["ASN_n"]["value"] = -0.46
			charges_groups["ASN_n"]["sele_string"] = "resname ASN"
			charges_groups["ASN_n"]["position"] = 3
			charges_groups["GLN_n"] = {}
			charges_groups["GLN_n"]["value"] = -0.46
			charges_groups["GLN_n"]["sele_string"] = "resname GLN"
			charges_groups["GLN_n"]["position"] = 3
			charges_groups["THR_n"] = {}
			charges_groups["THR_n"]["value"] = -0.31
			charges_groups["THR_n"]["sele_string"] = "resname THR"
			charges_groups["THR_n"]["position"] = 3
			charges_groups["SER_n"] = {}
			charges_groups["SER_n"]["value"] = -0.4
			charges_groups["SER_n"]["sele_string"] = "resname SER"
			charges_groups["SER_n"]["position"] = 3
			

	#use user supplied
	#-----------------
	else:	
		with open(args.chargesfilename) as f:
			lines = f.readlines()
		for l_index in range(0, len(lines)):
			line = lines[l_index]
			if line[-1] == "\n":
				line = line[:-1]
			l_content = line.split(',')
			if len(l_content) != 3:
				print "Error: wrong format on line " + str(l_index + 1) + " of " + str(args.chargesfilename) + ", see DESCRIPTION in TM_density --help."
				sys.exit(1)
			tmp_label = l_content[0]
			tmp_value = float(l_content[1])
			tmp_sele = l_content[2]
			if tmp_label not in charges_groups.keys():
				charges_groups[tmp_group] = {}
				charges_groups[tmp_group]["value"] = tmp_value
				charges_groups[tmp_group]["sele_string"] = tmp_sele

	return
def load_MDA_universe():
	
	global U
	global all_atoms
	global nb_atoms
	global nb_frames_xtc
	global frames_to_process
	global nb_frames_to_process
	global f_start
	global f_end
	global residues_list
	global water_pres
	global water_sele
	global leaflet_sele
	global upper_sele
	global lower_sele
	global charges_groups_pres
	f_start = 0
	charges_groups_pres = {}
		
	#load universe
	#-------------
	if args.xtcfilename == "no":
		print "\nLoading file..."
		U = Universe(args.grofilename)
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = 1
		frames_to_process = [0]
		nb_frames_to_process = 1
	else:
		print "\nLoading trajectory..."
		U = Universe(args.grofilename, args.xtcfilename)
		U_timestep = U.trajectory.dt
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = U.trajectory.numframes

		U.trajectory.rewind()
		#sanity check
		if U.trajectory[nb_frames_xtc-1].time/float(1000) < args.t_start:
			print "Error: the trajectory duration (" + str(U.trajectory.time/float(1000)) + "ns) is shorted than the starting stime specified (" + str(args.t_start) + "ns)."
			sys.exit(1)
		if U.trajectory.numframes < args.frames_dt:
			print "Warning: the trajectory contains fewer frames (" + str(nb_frames_xtc) + ") than the frame step specified (" + str(args.frames_dt) + ")."

		#create list of index of frames to process
		if args.t_end != -1:
			f_end = int((args.t_end*1000 - U.trajectory[0].time) / float(U_timestep))
			if f_end < 0:
				print "Error: the starting time specified is before the beginning of the xtc."
				sys.exit(1)
		else:
			f_end = nb_frames_xtc - 1
		if args.t_start != -1:
			f_start = int((args.t_start*1000 - U.trajectory[0].time) / float(U_timestep))
			if f_start > f_end:
				print "Error: the starting time specified is after the end of the xtc."
				sys.exit(1)
		if (f_end - f_start)%args.frames_dt == 0:
			tmp_offset = 0
		else:
			tmp_offset = 1
		frames_to_process = map(lambda f:f_start + args.frames_dt*f, range(0,(f_end - f_start)//args.frames_dt+tmp_offset))
		nb_frames_to_process = len(frames_to_process)

	#identify leaflets the lazy way (we assume a plane bilayer as we plot against z so no need to be fancy for now)
	#------------------------------
	leaflet_sele = U.selectAtoms("name PO4")
	tmp_lipids_avg_z = leaflet_sele.centerOfGeometry()[2]
	upper_sele = leaflet_sele.selectAtoms("prop z > " + str(tmp_lipids_avg_z))
	lower_sele = leaflet_sele.selectAtoms("prop z < " + str(tmp_lipids_avg_z))

	#add protein termini charge if present
	#-------------------------------------
	print "Identifying charges particles..."
	if not args.capped and not args.debug:
		prot = U.selectAtoms("protein")
		if prot.numberOfAtoms() > 0:
			prot.resnums()
			charges_groups["Nter"] = {}
			charges_groups["Cter"] = {}
			charges_groups["Nter"]["value"] = 1
			charges_groups["Cter"]["value"] = -1
			charged_residues.append("Nter")
			charged_residues.append("Cter")
			if not args.positions:
				charges_groups["Nter"]["sele_string"] = "resnum 1 and name BB"
				charges_groups["Cter"]["sele_string"] = "resnum " + str(int(prot.resnums()[-1])) + " and name BB"
			else:
				charges_groups["Nter"]["sele_string"] = "resnum 1"
				charges_groups["Nter"]["position"] = 0
				charges_groups["Cter"]["sele_string"] = "resnum " + str(int(prot.resnums()[-1]))
				charges_groups["Cter"]["position"] = 0
			print " added charge +1 to backbone bead on residue 1 of proteins"
			print " added charge -1 to backbone bead on residue " + str(int(prot.numberOfResidues())) + " of proteins"
			charges_groups_pres["Nter"] = True
			charges_groups_pres["Cter"] = True
		elif args.protein:
			print "Error: --protein specified but no proteins found."
			sys.exit(1)
	
	#create charged particles selections
	#-----------------------------------
	charge_pres = False
		
	for q in charges_groups.keys():
		charges_groups[q]["sele"] = U.selectAtoms(charges_groups[q]["sele_string"])
		if args.positions and q in charged_residues:
			charges_groups[q]["sele"] = charges_groups[q]["sele"].atoms[int(charges_groups[q]["position"]):int(charges_groups[q]["sele"].numberOfAtoms()):int(charges_groups[q]["sele"].residues[0].numberOfAtoms())]
		
		#debug mode
		if args.debug:
			print "atom names for charge group: ", q
			print charges_groups[q]["sele"].names()
		
		charges_groups[q]["nb"] = int(charges_groups[q]["sele"].numberOfAtoms())
		if charges_groups[q]["nb"] == 0:
			print " warning: charge selection string '" + str(charges_groups[q]["sele_string"]) + "' returned 0 atoms."
			charges_groups_pres[q] = False
		else:
			charge_pres = True
			charges_groups_pres[q] = True

	#debug mode
	if args.debug:
		charges_groups["BCQa"] = {}
		charges_groups["BCQa"]["sele"] = U.selectAtoms("name BCQa or name BHQa or name BSQa")
		charges_groups["BCQa"]["value"] = -1
		charges_groups_pres["BCQa"] = True
		charges_groups["BCQa"]["nb"] = int(charges_groups["BCQa"]["sele"].numberOfAtoms())
		charges_groups["BCQd"] = {}
		charges_groups["BCQd"]["sele"] = U.selectAtoms("name BCQd or name BHQd or name BSQd")
		charges_groups["BCQd"]["value"] = 1
		charges_groups_pres["BCQd"] = True
		charges_groups["BCQd"]["nb"] = int(charges_groups["BCQd"]["sele"].numberOfAtoms())
	
	if not charge_pres:
		print "Error: no charged particles found, try using --charges to supply correct charges definition."
		sys.exit(1)
		
	return

#=========================================================================================
# core functions
#=========================================================================================

def calculate_density(f_index, box_dim):
	
	global nb_x
	global nb_y
	global nb_z	
	global delta_x
	global delta_y
	global delta_z	
	global upper
	global lower

	#define bins
	tmp_bins_x = np.linspace(0,box_dim[0],args.sx+1)
	tmp_bins_y = np.linspace(0,box_dim[1],args.sx+1)
	tmp_bins_z = np.linspace(0,box_dim[2],args.sz+1)
	delta_x = tmp_bins_x[1]-tmp_bins_x[0]				#distance 2 cells along x (in Angstrom)
	delta_y = tmp_bins_y[1]-tmp_bins_y[0]	
	delta_z = tmp_bins_z[1]-tmp_bins_z[0]
	nb_x = int(np.floor(args.rc*10/float(delta_x))+1)	#nb of cells along x within the cutoff distance from another cell
	nb_y = int(np.floor(args.rc*10/float(delta_y))+1)
	nb_z = int(np.floor(args.rc*10/float(delta_z))+1)

	#store leaflets z position
	upper[f_index] = upper_sele.centerOfGeometry()[2]
	lower[f_index] = lower_sele.centerOfGeometry()[2]

	#calculate charge density in each bin
	#debug mode
	if args.debug:
		print "\ncharge groups used for stats:"
	for q in charges_groups.keys():
		if charges_groups_pres[q]:
			#debug mode
			if args.debug:
				print q
			
			#get coordinates
			tmp_coord = charges_groups[q]["sele"].coordinates()

			#transform coordinates into index of bins (the fanciness is to handle coords smaller/bigger than 0/U.dim)
			tmp_coord[:,0] = np.minimum(np.floor(tmp_coord[:,0]/float(delta_x)), np.remainder(np.floor(tmp_coord[:,0]/float(delta_x)),args.sx))
			tmp_coord[:,1] = np.minimum(np.floor(tmp_coord[:,1]/float(delta_y)), np.remainder(np.floor(tmp_coord[:,1]/float(delta_y)),args.sx))
			tmp_coord[:,2] = np.minimum(np.floor(tmp_coord[:,2]/float(delta_z)), np.remainder(np.floor(tmp_coord[:,2]/float(delta_z)),args.sz))
			tmp_coord = tmp_coord.astype(int)
			for qi in range(0,charges_groups[q]["nb"]):
				charge_density[tmp_coord[qi,0],tmp_coord[qi,1],tmp_coord[qi,2]] += charges_groups[q]["value"]
	return
def calculate_stats():
	global coords_x
	global coords_z
	global z_middle
	global upper_avg
	global lower_avg
	global potential
	global potential_1D
	global potential_2D
	global charge_density
	global charge_density_1D
	global charge_density_2D
	
	#calculate coords and leaflets positions
	#---------------------------------------
	upper_avg = np.average(upper)
	lower_avg = np.average(lower)
	z_middle = (upper_avg+lower_avg) / float(2)
	upper_avg -= z_middle
	lower_avg -= z_middle	
	tmp_coords_x = np.linspace(0,U.dimensions[0],args.sx+1)
	tmp_coords_z = np.linspace(0,U.dimensions[2],args.sz+1)
	coords_x = tmp_coords_x[0:args.sx] + (tmp_coords_x[1]-tmp_coords_x[0])/float(2) - U.dimensions[0]/float(2)
	coords_z = tmp_coords_z[0:args.sz] + (tmp_coords_z[1]-tmp_coords_z[0])/float(2) - z_middle

	#calculate average charge in each voxel
	#--------------------------------------
	charge_density = charge_density / float(nb_frames_to_process)

 	#calculate potential (in V)
 	#-------------------
	#calculate "distance" matrix: this matrix represents a cube, extending +/- nb_x (resp nb_y, nb_z) in each direction,
	#the value of each cells representing the distance from this cell to the center of the cube (the only trick is that
	#the indexing has to be done using positive integers only)
	#Note: the "distance" is everything in the potential expression which is not the charge (i.e. if we multiply the value
	#of a cell in the distance matrix by a charge we obtain the value of the potential in that cell created by this charge
	#if positioned in the center of the cube)
	tmp_dist = np.zeros((2*nb_x,2*nb_y,2*nb_z))			
	for nxx in range(0,2*nb_x+1):
		for nyy in range(0,2*nb_y+1):
			for nzz in range(0,2*nb_z+1):
				r = np.sqrt(((nxx-nb_x)*delta_x)**2 + ((nyy-nb_y)*delta_y)**2 + ((nzz-nb_z)*delta_z)**2) / float(10)	
				if r > 0 and r < args.rc:
					tmp_dist[nxx,nyy,nzz] = pot_conv * f_factor * (1/float(r) - 5/float(3*args.rc) + 5 * r**3 / float(3 * args.rc**4) - r**4 / float(args.rc**5))

	#open files
	if args.xtcfilename == "no":
		tmp_filename_dx = os.getcwd() + '/' + args.output_folder + '/' + str(args.grofilename[:-4]) + '_potential_3D_tmp.dx'
		tmp_filename_dx_charge = os.getcwd() + '/' + args.output_folder + '/' + str(args.grofilename[:-4]) + '_charge_3D_tmp.dx'
	else:
		tmp_filename_dx = os.getcwd() + '/' + args.output_folder + '/' + str(args.xtcfilename[:-4]) + '_potential_3D_tmp.dx'
		tmp_filename_dx_charge = os.getcwd() + '/' + args.output_folder + '/' + str(args.xtcfilename[:-4]) + '_charge_3D_tmp.dx'
	output_dx = open(tmp_filename_dx, 'w')
	output_dx_charge = open(tmp_filename_dx_charge, 'w')
	
	dx_counter = 0
	#browse each occupied voxel and add its potential contribution to its neighbours within rc
	for nx in range(0, args.sx):
		for ny in range(0, args.sx):
			for nz in range(0, args.sz):
				#display progress
				progress = '\r -processing voxel (' + str(nx) + ',' + str(ny) + ',' + str(nz) + ')       '  
				sys.stdout.flush()
				sys.stdout.write(progress)							

				#retrieve neighbours charge matrix				
				tmp_q = charge_density[np.remainder(np.arange(nx-nb_x,nx+nb_x),args.sx),:,:][:,np.remainder(np.arange(ny-nb_y,ny+nb_y),args.sx),:][:,:,np.remainder(np.arange(nz-nb_z,nz+nb_z),args.sz)]

				#multiply elements with distance matrix
				tmp_pot = np.multiply(tmp_q,tmp_dist)

				#add sum of contributions to current voxel
				#basically what we are doing is that we're multiplying each neighbouring charge by its "distance" from
				#the current point so as to calculate the potential it contributes to the current cell, and then we
				#just add all of these up
				potential[nx,ny,nz] += np.sum(tmp_pot)

				#append this to voxel
				output_dx.write(str(round(potential[nx,ny,nz],4)))
				output_dx_charge.write(str(round(charge_density[nx,ny,nz],4)))
				
				if dx_counter == 2:
					output_dx.write("\n")
					output_dx_charge.write("\n")
					dx_counter = 0
				else:
					output_dx.write(" ")
					output_dx_charge.write("\n")
					dx_counter += 1
	
	#close temporary dx files
	output_dx.write("attribute \"dep\" string \"positions\"\n")
	output_dx.write("object \"Electrostatic potential (V)\" class field\n")
	output_dx.write("component \"positions\" value 1\n")
	output_dx.write("component \"connections\" value 2\n")
	output_dx.write("component \"data\" value 3\n")
	output_dx.close()
	output_dx_charge.write("attribute \"dep\" string \"positions\"\n")
	output_dx_charge.write("object \"Charge density (e.A-3)\" class field\n")
	output_dx_charge.write("component \"positions\" value 1\n")
	output_dx_charge.write("component \"connections\" value 2\n")
	output_dx_charge.write("component \"data\" value 3\n")
	output_dx_charge.close()

	#calculate average potential
	#---------------------------
	for nz in range(0, args.sz):
		potential_1D[nz] = np.average(potential[:,:,nz])
		for nx in range(0, args.sx):
			potential_2D[nx,nz] = np.average(potential[nx,:,nz])

	#calculate charge density
	#------------------------
	charge_density = charge_density / float(delta_x*delta_y*delta_z)
	for nz in range(0, args.sz):
		charge_density_1D[nz] = np.average(charge_density[:,:,nz])
		for nx in range(0, args.sx):
			charge_density_2D[nx,nz] = np.average(charge_density[nx,:,nz])

	return

#=========================================================================================
# outputs
#=========================================================================================

#xvg files
def write_xvg_charges():
	
	#open files
	filename_xvg = os.getcwd() + '/' + args.output_folder + '/' + str(args.xtcfilename[:-4]) + '_charge_1D.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [electrostatic potential profile - written by cg_potential_vol v" + str(version_nb) + "]\n")
	output_xvg.write("#  -> nb of slices x and y: " + str(args.sx) + "\n")
	output_xvg.write("#  -> nb of slices z: " + str(args.sz) + "\n")
	output_xvg.write("#  -> slices volume: " + str(round(delta_x*delta_y*delta_z,2)) + ") (Angstrom3)\n")
	output_xvg.write("# nb of frames which contributed to this profile:\n")
	output_xvg.write("# -> weight = " + str(nb_frames_to_process) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Charge densisty along z\"\n")
	output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
	output_xvg.write("@ yaxis label \"charge density ($e^{-1}.\AA^{3}$)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 1\n")
	output_xvg.write("@ s0 legend \"total\"\n")
	
	#data
	for nz in range(0,args.sz):
		results = str(coords_z[nz]) + "	" + "{:.6e}".format(charge_density_1D[nz])
		output_xvg.write(results + "\n")	
	output_xvg.close()
		
	return
def write_xvg_potential():
	
	#open files
	filename_xvg = os.getcwd() + '/' + args.output_folder + '/' + str(args.xtcfilename[:-4]) + '_potential_1D.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [electrostatic potential profile - written by cg_potential_vol v" + str(version_nb) + "]\n")
	output_xvg.write("#  -> nb of slices x and y: " + str(args.sx) + "\n")
	output_xvg.write("#  -> nb of slices z: " + str(args.sz) + "\n")
	output_xvg.write("#  -> slices volume: " + str(round(delta_x*delta_y*delta_z,2)) + " (Angstrom3)\n")
	output_xvg.write("# nb of frames which contributed to this profile:\n")
	output_xvg.write("# -> weight = " + str(nb_frames_to_process) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Electrostatic potential along z\"\n")
	output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
	output_xvg.write("@ yaxis label \"potential (V)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 1\n")
	output_xvg.write("@ s0 legend \"total\"\n")
	
	#data
	for nz in range(0,args.sz):
		results = str(coords_z[nz]) + "	" + "{:.6e}".format(potential_1D[nz])
		output_xvg.write(results + "\n")	
	output_xvg.close()

		
	return

#graphs
def graph_charges():
		
	#filenames
	filename_svg = os.getcwd() + '/' + args.output_folder + '/' + str(args.xtcfilename[:-4]) + '_charges_1D.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Charge density profile along z")

	#plot data
	ax = fig.add_subplot(111)
	plt.plot(coords_z, charge_density_1D, color = 'k', linewidth = 2)
	plt.vlines(lower_avg, min(charge_density_1D), max(charge_density_1D), linestyles = 'dashed')
	plt.vlines(upper_avg, min(charge_density_1D), max(charge_density_1D), linestyles = 'dashed')
	plt.vlines(0, min(charge_density_1D), max(charge_density_1D), linestyles = 'dashdot')
	plt.hlines(0, min(coords_z), max(coords_z))
	plt.xlabel('z distance to bilayer center ($\AA$)')
	plt.ylabel('average charge density ($e.\AA^{-3}$)')
	
	#save figure
	ax.set_xlim(min(coords_z), max(coords_z))
	#ax.set_ylim(min_density_charges, max_density_charges)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
	ax.xaxis.labelpad = 10
	ax.yaxis.labelpad = 10
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.85)
	fig.savefig(filename_svg)
	plt.close()

	#2D profile
	#----------
	
	#filenames
	filename_svg = os.getcwd() + '/' + args.output_folder + '/' + str(args.xtcfilename[:-4]) + '_charges_2D.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Electrostatic profile slice")

	#plot data
	ax = fig.add_subplot(111)

	#rotate data so that the x is horizontal and z vertical after imshow plotting
	charge_density_2D_oriented = np.zeros((args.sz,args.sx))
	for nx in range(0, args.sx):
		for nz in range(0, args.sz):
			charge_density_2D_oriented[nz,nx] = charge_density_2D[nx,args.sz-1-nz]

	#plot data
	ax = fig.add_subplot(111)
	im = plt.imshow(charge_density_2D_oriented, extent = [min(coords_x),max(coords_x),min(coords_z),max(coords_z)], cmap = matplotlib.cm.jet_r)
	plt.vlines(lower_avg, min(charge_density_2D[:,0]), max(charge_density_2D[:,0]), linestyles = 'dashed')
	plt.vlines(upper_avg, min(charge_density_2D[:,0]), max(charge_density_2D[:,0]), linestyles = 'dashed')
	plt.vlines(0, min(charge_density_2D[:,0]), max(charge_density_2D[:,0]), linestyles = 'dashdot')
	plt.xlabel('x ($\AA$)')
	plt.ylabel('z distance to bilayer center ($\AA$)')
	
	#color bar
	cax = fig.add_axes([0.85, 0.26, 0.025, 0.48])
	cbar = fig.colorbar(im, orientation='vertical', cax=cax)
	cbar.ax.tick_params(axis='y', direction='out')
	cbar.set_label(r'charge density ($e.\AA^{-3}$)')
		
	#save figure
	ax.set_xlim(min(coords_x), max(coords_x))
	ax.set_ylim(min(coords_z), max(coords_z))
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
	ax.xaxis.labelpad = 10
	ax.yaxis.labelpad = 10
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.1, right = 0.8)
	fig.savefig(filename_svg)
	plt.close()

	return
def graph_potential():
			
	#1D profile
	#----------
	#filenames
	filename_svg = os.getcwd() + '/' + args.output_folder + '/' + str(args.xtcfilename[:-4]) + '_potential_1D.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Electrostatic profile along z")

	#plot data
	ax = fig.add_subplot(111)
	plt.plot(coords_z, potential_1D, color = 'k', linewidth = 2)
	plt.vlines(lower_avg, min(potential_1D), max(potential_1D), linestyles = 'dashed')
	plt.vlines(upper_avg, min(potential_1D), max(potential_1D), linestyles = 'dashed')
	plt.vlines(0, min(potential_1D), max(potential_1D), linestyles = 'dashdot')
	plt.hlines(0, min(coords_z), max(coords_z))
	plt.hlines(0, min(coords_z), max(coords_z))
	plt.xlabel('z distance to bilayer center ($\AA$)')
	plt.ylabel('electrostatic potential (V)')
	
	#save figure
	#ax.set_xlim(min(coords_z), max(coords_z))
	ax.set_xlim(-60, 60)
	#ax.set_ylim(min_density_charges, max_density_charges)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=9))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
	ax.xaxis.labelpad = 10
	ax.yaxis.labelpad = 10
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.85)
	fig.savefig(filename_svg)
	plt.close()


	#2D profile
	#----------
	#filenames
	filename_svg = os.getcwd() + '/' + args.output_folder + '/' + str(args.xtcfilename[:-4]) + '_potential_2D.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Electrostatic profile slice")

	#rotate data so that the x is horizontal and z vertical after imshow plotting
	potential_2D_oriented = np.zeros((args.sz,args.sx))
	for nx in range(0, args.sx):
		for nz in range(0, args.sz):
			potential_2D_oriented[nz,nx] = potential_2D[nx,args.sz-1-nz]

	#plot data
	ax = fig.add_subplot(111)
	#im = plt.imshow(potential_2D_oriented, extent = [min(coords_x),max(coords_x),min(coords_z),max(coords_z)], cmap = matplotlib.cm.jet_r, vmin = -0.08, vmax = 0.04)
	im = plt.imshow(potential_2D_oriented[19:140,:], extent = [min(coords_x),max(coords_x),-60,60], cmap = matplotlib.cm.jet_r, vmin = -0.06, vmax = 0.03)
	plt.vlines(lower_avg, min(potential_2D[:,0]), max(potential_2D[:,0]), linestyles = 'dashed')
	plt.vlines(upper_avg, min(potential_2D[:,0]), max(potential_2D[:,0]), linestyles = 'dashed')
	plt.vlines(0, min(potential_2D[:,0]), max(potential_2D[:,0]), linestyles = 'dashdot')
	plt.xlabel('x ($\AA$)')
	plt.ylabel('z distance to bilayer center ($\AA$)')
	
	#color bar
	cax = fig.add_axes([0.85, 0.26, 0.025, 0.48])
	cbar = fig.colorbar(im, orientation='vertical', cax=cax)
	cbar.ax.tick_params(axis='y', direction='out')
	cbar.set_label(r'potential (V)')
		
	#save figure
	ax.set_xlim(min(coords_x), max(coords_x))
	#ax.set_ylim(min(coords_z), max(coords_z))
	ax.set_ylim(-60,60)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=9))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
	ax.xaxis.labelpad = 10
	ax.yaxis.labelpad = 10
	ax.tick_params(axis='x', direction='out')
	ax.tick_params(axis='y', direction='out')
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.1, right = 0.8)
	fig.savefig(filename_svg)
	plt.close()

	return

#Opendx file
def write_dx_potential():

	#open files and read in data
	if args.xtcfilename == "no":
		filename_dx = os.getcwd() + '/' + args.output_folder + '/' + str(args.grofilename[:-4]) + '_potential_3D'
	else:
		filename_dx = os.getcwd() + '/' + args.output_folder + '/' + str(args.xtcfilename[:-4]) + '_potential_3D'
	with file(filename_dx + '_tmp.dx', 'r') as f:
		dx_data = f.read()
	output_dx = open(filename_dx + '.dx', 'w')
	
	#general header
	output_dx.write("#electrostatic potential - written by cg_potential_vol v" + str(version_nb) + "]\n")
	output_dx.write("# -> nb frames processed: " + str(nb_frames_to_process) + "\n")
	output_dx.write("# -> epsilon_r: " + str(args.er) + "\n")
	output_dx.write("# -> r_switch: " + str(args.rs) + "\n")
	output_dx.write("# -> r_cutoff: " + str(args.rc) + "\n")

	#array info
	output_dx.write("object 1 class gridpositions counts " + str(args.sx) + " " + str(args.sx) + " " + str(args.sz) + "\n")
	output_dx.write("origin 0 0 0\n")
	output_dx.write("delta " + str(round(delta_x,4)) + " 0 0\n")
	output_dx.write("delta 0 " + str(round(delta_y,4)) + " 0\n")
	output_dx.write("delta 0 0 " + str(round(delta_z,4)) + "\n")
	output_dx.write("object 2 class gridconnections counts " + str(args.sx) + " " + str(args.sx) + " " + str(args.sz) + "\n")
	output_dx.write("object 3 class array type double rank 0 items " + str(int(args.sx*args.sx*args.sz)) + " data follows\n")
	
	#write data
	output_dx.write(dx_data)
	output_dx.close()
	os.remove(filename_dx + '_tmp.dx')
	
	return

def write_dx_charge():
	#open files and read in data
	if args.xtcfilename == "no":
		filename_dx = os.getcwd() + '/' + args.output_folder + '/' + str(args.grofilename[:-4]) + '_charge_3D'
	else:
		filename_dx = os.getcwd() + '/' + args.output_folder + '/' + str(args.xtcfilename[:-4]) + '_charge_3D'
	with file(filename_dx + '_tmp.dx', 'r') as f:
		dx_data = f.read()
	output_dx = open(filename_dx + '.dx', 'w')
	
	#general header
	output_dx.write("#electrostatic potential - written by cg_potential_vol v" + str(version_nb) + "]\n")
	output_dx.write("# -> nb frames processed: " + str(nb_frames_to_process) + "\n")
	output_dx.write("# -> epsilon_r: " + str(args.er) + "\n")
	output_dx.write("# -> r_switch: " + str(args.rs) + "\n")
	output_dx.write("# -> r_cutoff: " + str(args.rc) + "\n")

	#array info
	output_dx.write("object 1 class gridpositions counts " + str(args.sx) + " " + str(args.sx) + " " + str(args.sz) + "\n")
	output_dx.write("origin 0 0 0\n")
	output_dx.write("delta " + str(round(delta_x,4)) + " 0 0\n")
	output_dx.write("delta 0 " + str(round(delta_y,4)) + " 0\n")
	output_dx.write("delta 0 0 " + str(round(delta_z,4)) + "\n")
	output_dx.write("object 2 class gridconnections counts " + str(args.sx) + " " + str(args.sx) + " " + str(args.sz) + "\n")
	output_dx.write("object 3 class array type double rank 0 items " + str(int(args.sx*args.sx*args.sz)) + " data follows\n")
	
	#write data
	output_dx.write(dx_data)
	output_dx.close()
	os.remove(filename_dx + '_tmp.dx')

	return

##########################################################################################
# ALGORITHM
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================
set_charges()
load_MDA_universe()
upper = np.zeros(nb_frames_to_process)
lower = np.zeros(nb_frames_to_process)

#=========================================================================================
# generate data
#=========================================================================================
print "\nCalculating charge density..."
#case: structure only
#--------------------
if args.xtcfilename=="no":
	calculate_density(0,U.trajectory.ts.dimensions)
#case: browse xtc frames
#-----------------------
else:
	for f_index in range(0,nb_frames_to_process):
		ts = U.trajectory[frames_to_process[f_index]]
		progress = '\r -processing frame ' + str(f_index+1) + '/' + str(nb_frames_to_process) + ' (every ' + str(args.frames_dt) + ' frame(s) from frame ' + str(f_start) + ' to frame ' + str(f_end) + ' out of ' + str(nb_frames_xtc) + ')      '  
		sys.stdout.flush()
		sys.stdout.write(progress)							
		calculate_density(f_index, U.trajectory.ts.dimensions)
	print ""

#=========================================================================================
# process data
#=========================================================================================
print "\nCalculating electrostatic potential..."
calculate_stats()

#=========================================================================================
# produce outputs
#=========================================================================================
print "\n\nWriting outputs..."
write_xvg_charges()
write_xvg_potential()
write_dx_potential()
write_dx_charge()
graph_charges()
graph_potential()	
	
#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check output in ./" + args.output_folder + "/"
print ""
sys.exit(0)
