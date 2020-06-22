#!/usr/bin/env python

### ____________________________________________________________________ ###
###                                                                      ###
### OptiFind (optifind.py) - Script for catalogue-based source finding   ###
###                          with SoFiA 2                                ###
### Copyright (C) 2020 Tobias Westmeier                                  ###
### ____________________________________________________________________ ###
###                                                                      ###
### Address:  Tobias Westmeier                                           ###
###           ICRAR M468                                                 ###
###           The University of Western Australia                        ###
###           35 Stirling Highway                                        ###
###           Crawley WA 6009                                            ###
###           Australia                                                  ###
###                                                                      ###
### E-mail:   tobias.westmeier [at] uwa.edu.au                           ###
### ____________________________________________________________________ ###
###                                                                      ###
### This program is free software: you can redistribute it and/or modify ###
### it under the terms of the GNU General Public License as published by ###
### the Free Software Foundation, either version 3 of the License, or    ###
### (at your option) any later version.                                  ###
###                                                                      ###
### This program is distributed in the hope that it will be useful,      ###
### but WITHOUT ANY WARRANTY; without even the implied warranty of       ###
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         ###
### GNU General Public License for more details.                         ###
###                                                                      ###
### You should have received a copy of the GNU General Public License    ###
### along with this program. If not, see http://www.gnu.org/licenses/.   ###
### ____________________________________________________________________ ###
###                                                                      ###


import os
import sys
import math
from astropy.io import fits
from astropy.wcs import WCS



# Error handling
def exit_error(message):
	sys.stderr.write("ERROR: {0}\n".format(message));
	sys.exit(1);



# Print usage instructions
def print_instructions():
	sys.stdout.write("\n");
	sys.stdout.write(" Usage:\n   optifind.py <par_file> <source_list> <r_spat> <r_spec> [<sofia_exe>]\n");
	sys.stdout.write("\n");
	sys.stdout.write(" Arguments:\n");
	sys.stdout.write("   <par_file>     Name of the SoFiA 2 control parameter file to be used.\n");
	sys.stdout.write("   <source_list>  Name of the input source catalogue file.\n");
	sys.stdout.write("   <r_spat>       Spatial radius of the sub-region in pixels.\n");
	sys.stdout.write("   <r_spec>       Spectral radius of the sub-region in channels.\n");
	sys.stdout.write("   <sofia_exe>    Optional name of the SoFiA 2 executable. Default: sofia.\n");
	sys.stdout.write("\n");
	sys.stdout.write(" OptiFind serves as a wrapper script around the SoFiA 2 source finding pipeline\n");
	sys.stdout.write(" to allow source finding  on multiple sub-regions of a data cube centred on the\n");
	sys.stdout.write(" positions from a user-supplied source catalogue. The user will need to provide\n");
	sys.stdout.write(" a template SoFiA 2 parameter file that will be used in each source finding run\n");
	sys.stdout.write(" spawned by OptiFind.\n");
	sys.stdout.write("\n");
	sys.stdout.write(" The template file must specify the input data cube to be searched and can also\n");
	sys.stdout.write(" define an output file name which will be used as the base name for all output.\n");
	sys.stdout.write(" In addition,  the user must specify a source catalogue containing the position\n");
	sys.stdout.write(" of each source to be searched. The catalogue must specify the world coordinate\n");
	sys.stdout.write(" position of each source on a separate line. The following, comma-separated pa-\n");
	sys.stdout.write(" rameters must be supplied with each source:\n");
	sys.stdout.write("\n");
	sys.stdout.write("   id, coord_1, coord_2, coord_3, ...\n");
	sys.stdout.write("\n");
	sys.stdout.write(" Here, id is a unique source ID  that will be used as an identifier  for output\n");
	sys.stdout.write(" products,  while coord_n denotes the coordinates in each dimension of the data\n");
	sys.stdout.write(" cube.  The coordinate must be specified in the raw units  of the data cube  as\n");
	sys.stdout.write(" specified in the FITS header.  A coordinate value must be given  for each axis\n");
	sys.stdout.write(" of the cube in the correct order.\n");
	sys.stdout.write("\n");
	sys.stdout.write(" For example, if a cube has four axes  (right ascension, declination, frequency\n");
	sys.stdout.write(" and Stokes),  then four coordinate values need to be provided for each source,\n");
	sys.stdout.write(" and they must be given in the native units  specified in the header,  e.g. de-\n");
	sys.stdout.write(" grees for right ascension and Hz for frequency.\n");
	sys.stdout.write("\n");
	sys.stdout.write(" Separate output catalogues and products will be created  for each SoFiA 2 run.\n");
	sys.stdout.write(" Their base name will be either \"optifind\" + suffix or output.filename + suffix\n");
	sys.stdout.write(" depending on whether or not an output file name  was defined  in the parameter\n");
	sys.stdout.write(" file.  The suffix will be an underscore  followed by the source ID provided in\n");
	sys.stdout.write(" the catalogue file.\n");
	sys.stdout.write("\n");
	sys.stdout.write(" In addition to the individual  output catalogues from each run,  OptiFind will\n");
	sys.stdout.write(" also create a single, merged catalogue called  \"optifind_merged_catalogue.txt\"\n");
	sys.stdout.write(" in the same output directory.  Note that this feature is currently only avail-\n");
	sys.stdout.write(" able for plain-text ASCII catalogues, but not for XML or SQL catalogues.\n");
	sys.stdout.write("\n");
	return;



# Check for command-line arguments
if(len(sys.argv) < 5):
	print_instructions();
	sys.exit(0);
radius_spat = int(sys.argv[3]);
radius_spec = int(sys.argv[4]);
if(len(sys.argv) > 5): sofia_2_executable = sys.argv[5];
else: sofia_2_executable = "sofia";


# Read parameter file
pars = {};
try:
	with open(sys.argv[1]) as fp:
		for line in fp:
			line = line.strip();
			if(line and line[0] != "#" and "=" in line):
				(key, value) = line.split("=", 1);
				key = key.strip();
				value = value.strip();
				if(value):
					comment = value.find("#");
					if(comment < 0): pars[key] = value;
					else: pars[key] = value.split("#", 1)[0].strip();
				else: pars[key] = "";
except:
	exit_error("Failed to read parameter file: " + sys.argv[1]);

if(len(pars) == 0): exit_error("No valid parameter settings found.");


# Read WCS from input cube header
try:
	wcs = WCS(pars["input.data"]);
	axis_size = wcs.array_shape[::-1];
	# NOTE: Need to reverse order here, as otherwise inconsistent with order of wcs.axis_type_names!!!
	naxes = len(axis_size);
	axis_type = wcs.axis_type_names;
except:
	exit_error("Failed to read WCS from input data cube.");


# Work out spatial and spectral axes
axes = [-1, -1, -1];
for i in range(naxes):
	if(axis_type[i] == "RA" or axis_type[i] == "GLON"): axes[0] = i;
	elif(axis_type[i] == "DEC" or axis_type[i] == "GLAT"): axes[1] = i;
	elif(axis_type[i] == "FREQ" or axis_type[i] == "VELO" or axis_type[i] == "VRAD" or axis_type[i] == "VOPT" or axis_type[i] == "FELO"): axes[2] = i;
if(axes[0] == -1 or axes[1] == -1 or axes[2] == -1): exit_error("Failed to identify spatial and/or spectral axis of data cube.");


# Read input catalogue
sources = [];
n_cols = 0;
try:
	with open(sys.argv[2]) as fp:
		for line in fp:
			line = line.strip();
			if(line and line[0] != "#"):
				cols = line.split(",");
				if(n_cols != len(cols)):
					if(n_cols > 0): exit_error("Variable number of catalogue columns encountered.");
					else: n_cols = len(cols);
				cols = [item.strip() for item in cols];
				sources.append(cols);
except:
	exit_error("Failed to read input catalogue: " + sys.argv[2]);

if(len(sources) == 0): exit_error("No sources found in input catalogue.");
if(n_cols != naxes + 1): exit_error("Data cube is {:d}D, but {:d} coordinate values given in catalogue.".format(naxes, n_cols - 1));


# Loop over all sources
cat_names = [];
for src in sources:
	# Convert source position into pixel coordinates
	coord_wld = [[float(item) for item in src[1:]]]
	coord_pix = wcs.wcs_world2pix(coord_wld, 0)[0];
	
	# Work out region
	x_min = max(0, int(math.floor(coord_pix[axes[0]] - radius_spat)));
	x_max = min(axis_size[axes[0]] - 1, int(math.ceil(coord_pix[axes[0]] + radius_spat)));
	y_min = max(0, int(math.floor(coord_pix[axes[1]] - radius_spat)));
	y_max = min(axis_size[axes[1]] - 1, int(math.ceil(coord_pix[axes[1]] + radius_spat)));
	z_min = max(0, int(math.floor(coord_pix[axes[2]] - radius_spec)));
	z_max = min(axis_size[axes[2]] - 1, int(math.ceil(coord_pix[axes[2]] + radius_spec)));
	
	# Check if region within bounds
	if(x_max - x_min >= radius_spat and y_max - y_min >= radius_spat and z_max - z_min >= radius_spec):
		# Copy and adjust parameter settings
		pars_copy = pars.copy();
		
		pars_copy["input.region"] = "{:d}, {:d}, {:d}, {:d}, {:d}, {:d}".format(x_min, x_max, y_min, y_max, z_min, z_max);
		print("Processing source \"{}\"".format(src[0]));
		print("  Region: {}".format(str(pars_copy["input.region"])));
		
		suffix = "_" + src[0].replace(" ", "_").replace("/", "-");
		
		if(pars_copy["output.filename"]):
			pars_copy["output.filename"] += suffix;
		else:
			pars_copy["output.filename"] = "optifind" + suffix;
		cat_names.append(pars_copy["output.filename"]);
		
		pars_copy["parameter.offset"] = "true";
		
		# Write temporary parameter file
		try:
			with open("sofia_optifind_parameter_file.tmp", "w") as fp:
				for key in pars_copy:
					fp.write("{}\t=\t{}\n".format(key, pars_copy[key]));
		except:
			exit_error("Failed to write temporary SoFiA 2 parameter file. Is the\n       current working directory write-protected?");
		
		# Run SoFiA and delete temporary parameter file again
		os.system("{} sofia_optifind_parameter_file.tmp".format(sofia_2_executable));
		os.system("rm -f sofia_optifind_parameter_file.tmp");
	else:
		print("Skipping source \"{}\": outside bounds.".format(src[0]));


# Merge output catalogues
n = len(cat_names);
if(n):
	# Figure out where the catalogue files are
	output_dir = pars["output.directory"];
	fmt_plain  = pars["output.writeCatASCII"];
	
	if(output_dir and output_dir[-1] != "/"): output_dir += "/";
	
	if(fmt_plain):
		sys.stdout.write("Creating merged ASCII catalogue:\n  " + output_dir + "optifind_merged_catalogue.txt\n");
		header  = "";
		content = "";
		need_header = True;
		
		for item in cat_names:
			# Read catalogue
			try:
				with open(output_dir + item + "_cat.txt") as fp:
					for line in fp:
						if(need_header and line and line[0] == "#"): header += line;
						if(line and not line.isspace() and line[0] != "#"): content += line;
			except:
				exit_error("Failed to read SoFiA 2 catalogue: " + output_dir + item + "_cat.txt");
			
			if(need_header and header): need_header = False;
		
		# Create output catalogue
		try:
			with open(output_dir + "optifind_merged_catalogue.txt", "w") as fp:
				fp.write(header + "\n");
				fp.write(content);
		except:
			exit_error("Failed to write merged output catalogue. Is the\n       current working directory write-protected?");
