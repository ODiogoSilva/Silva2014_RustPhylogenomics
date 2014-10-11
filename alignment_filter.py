#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Alignment_filter.py
#  
#  Copyright 2012 Unknown <diogo@arch>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  

import argparse
import subprocess
import sys
import pickle
import Alignment
import time
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Filters alignment files")

parser.add_argument("-in", dest="infile", nargs="+", required=True, help="Provide the input file name. If multiple "
					"files are provided, please separated the names with spaces")
parser.add_argument("-if", dest="InputFormat", default="fasta", choices=["fasta", "nexus", "phylip"], help="Format of "
					"the input file(s) (default is '%(default)s')")
parser.add_argument("-g", dest="gap", default="-", help="Symbol for gap (default is '%(default)s')")
parser.add_argument("-m", dest="missing", default="X", help="Symbol for missing data (default is '%(default)s')")
parser.add_argument("-tg", dest="threshold_gap", nargs=1, help="Threshold for the maximum proportion of gaps allowed "
					"in each column")
parser.add_argument("-tm", dest="threshold_missing", nargs=1, help="Threshold for the maximum proportion of missing "
					"data allowed in each column")
parser.add_argument("-o", dest="outfile", help="Name of the output file")
parser.add_argument("-plot", dest="plot", action="store_const", const="True", help="Create a plot of the missing and "
					"gap data for the alignment")

arg = parser.parse_args()


def loading(current_state, size, prefix, width):
	""" Function that prints the loading progress of the script """
	percentage = int(((current_state + 1) / size) * 100)
	complete = int(width * percentage * 0.01)
	if percentage == 100:
		sys.stdout.write("\r%s [%s%s] %s%% -- Done!\n" % (prefix, "#" * complete, "." * (width - complete), percentage))
	else:
		sys.stdout.write("\r%s [%s%s] %s%%" % (prefix, "#" * complete, "." * (width-complete), percentage))
	sys.stdout.flush()


def dataset_creator(infile_list):
	""" Deals with input files. Creates one instance in case there is only one input alignment, or another instance
	in case there are multiple input files"""

	if len(infile_list) == 1:
		alignment_object = Alignment.Alignment("".join(infile_list))

	else:
		alignments = Alignment.AlignmentList(infile_list)
		alignment_object = alignments.concatenate()

	alignment_dic = alignment_object.alignment
	alignment_order = alignment_object.iter_sequences()

	return alignment_dic, alignment_order
	
def data_breaker(storage, alignment_length, break_range=20000):
	""" Function that will break the storage into smaller parts (of 50 000 positions each) and each one will be processed separately and then merged in the end. This is only necessary for very large alignments, in which case the partitioning prevents the machine from running out of memory, and it's considerably faster """
	subpart_list = []
	previous_value = 0
	for i in range(break_range,alignment_length,break_range):
		temp_storage = dict((taxa,seq[previous_value:i]) for taxa, seq in storage.items())
		pickle.dump(temp_storage, open ("alignment_subpart_"+str(i), "wb"))
		subpart_list.append("alignment_subpart_"+str(i))
		previous_value = i
	else:
		temp_storage = dict((taxa,seq[previous_value:alignment_length]) for taxa, seq in storage.items())
		pickle.dump(temp_storage, open ("alignment_subpart_"+str(alignment_length), "wb"))
		subpart_list.append("alignment_subpart_"+str(alignment_length))
	return subpart_list

def extremitiesFiller (sequence, missing_symbol, gap_symbol,extremities="both"):
	""" Function used to substitute the gap filled extremities of aligned sequences with the missing data symbol """
	nucl,nucl_rev = 0,-1 # These are the counter to iterate over the sequence. The nucl counter starts at the begining, going forward. The nucl_rev starts at the end of the sequence, going backwards
	sequence_list = list(sequence)
	if len(set(sequence_list)) == 1: # Checks if there are more than 1 type of character. If a sequence is composed only with, let's say, gaps, then the next line replaces it all for the missing data symbol
		sequence_final = missing_symbol*len(sequence) #
		return sequence_final
	if extremities == "both" or extremities == "begin":
		while sequence_list[nucl] == gap_symbol:
			sequence_list[nucl] = missing_symbol
			nucl += 1
	if extremities == "both" or extremities == "end":
		while sequence_list[nucl_rev] == gap_symbol:
			sequence_list[nucl_rev] = missing_symbol
			nucl_rev -= 1
		sequence_final = "".join(sequence_list)	
	elif extremities == "no":
		sequence_final = sequence
	return sequence_final
	
def columnFilter (sequence_dict, missing_threshold, gap_threshold, missing_symbol, gap_symbol, multi="no"):
	""" This function checks the columns of the alignment and either removes or masks columns that are not up to the defined quality standards """
	
	mod_sequence_list = dict((header,list(sequence)) for header,sequence in sequence_dict.items())
	sequence_length = len(list(sequence_dict.values())[0])
	taxa_number = len(sequence_dict.values())
	
	excluded_missing, excluded_gaps = [], []
	total_missing, total_gap = [],[]
	
	current_position = 0
	
	# Creating the column variable in which quality checks will be made
	for position in range(sequence_length-1,-1,-1):
		if multi == "no": # When there are multiple files being processed, do not show the progress bar for each alignment
			loading(current_position,sequence_length,"Processing alignment", 50)
		current_column = []
		for char in [sequence[position] for sequence in sequence_dict.values()]:
			current_column.append(char)
			
		# Making quality checks
		proportion_missing_data = (float(current_column.count(missing_symbol))/float(taxa_number))*float(100)
		proportion_gaps = (float(current_column.count(gap_symbol))/float(taxa_number))*float(100)
		total_gap_proportion = proportion_missing_data+proportion_gaps
		total_missing.append((position,proportion_missing_data))
		total_gap.append((position, proportion_gaps))
		
		# Comparing current values to thresholds
		if proportion_missing_data > float(missing_threshold):
			excluded_missing.append(position)
			[sequence.pop(position) for sequence in mod_sequence_list.values()]
		elif total_gap_proportion > float(gap_threshold):
			excluded_gaps.append(position)
			[sequence.pop(position) for sequence in mod_sequence_list.values()]
		current_position += 1
	mod_sequence_final = dict((header,"".join(sequence)) for header,sequence in mod_sequence_list.items())
	return mod_sequence_final, excluded_missing[::-1], excluded_gaps[::-1], sequence_length, total_missing, total_gap
	

def file_filter (sequences_dict, missing_symbol, gap_symbol, missing_threshold, gap_threshold, multi="no", extremities_2="both"):
	""" Main function that treats gaps and applies thresholds on the current alignment """
	new_storage = {}
	# Adds missing data to extremities
	for header, sequence in sequences_dict.items():
		new_storage[header] = extremitiesFiller (sequence, missing_symbol, gap_symbol, extremities=extremities_2)
	# Filters sequences according to the given thresholds
	filter_sequences, excluded_missing, excluded_gaps, sequence_length, total_missing, total_gap = columnFilter (new_storage, missing_threshold, gap_threshold, missing_symbol,gap_symbol, multi)
	return filter_sequences, excluded_missing, excluded_gaps, sequence_length, total_missing, total_gap
	
def list_fusion (inlist):
	""" Function that fuses sequential number in a list (e.g. [1,2,3,4,6,8,9] into [1-4,5,8,9]) """
	prev_number = -2
	final_list = []
	temp_list = []
	flag = 0
	for number in inlist:
			num = number - 1
			if num != prev_number and flag == 0:
				final_list.append(number)
			elif num == prev_number and flag == 0:
				final_list.pop(-1)
				temp_list.append(prev_number)
				temp_list.append(number)
				flag = 1
			elif num == prev_number and flag == 1 and number == inlist[-1]:
				final_list.append(str(temp_list[0])+"-"+str(number))
			elif num == prev_number and flag == 1:
				temp_list.append(number)
			elif num != prev_number and flag == 1:
				final_list.append(str(temp_list[0])+"-"+str(temp_list[-1]))
				final_list.append(number)
				temp_list = []
				flag = 0
			prev_number = number
	return final_list
				
def logger (infile, excluded_missing, excluded_gaps, missing_threshold, gap_threshold):
	""" Collects information from the missing and gap data filtered from a single file and saves on a log file """
	log_code = infile.split(".")[0]+"_filterLog.txt"
	log_handle = open(log_code,"w")
	excluded_missing_fused = list_fusion(excluded_missing)
	excluded_gap_fused = list_fusion(excluded_gaps)
	log_handle.write("This is Alignment_filter.py v0.1 run at %s\n\nProcessing file %s\n" % (time.asctime(time.localtime()),infile))
	log_handle.write("\nExcluded positions due to missing data (threshold was %s):\nTotal: %s\n" % (missing_threshold,len(excluded_missing)))
	[log_handle.write("%s; " % position) for position in excluded_missing_fused]
	log_handle.write("\nExcluded positions due to excessive gap data (threshold was %s):\nTotal: %s\n" % (gap_threshold,len(excluded_gaps)))
	[log_handle.write("%s; " % position) for position in excluded_gap_fused]
	log_handle.close()
	
def master_log (infile, missing_number, gap_number, missing_threshold, gap_threshold):
	""" Logs the missing and gap data filtered for all files that were processed """
	log_handle = open("Master_log.txt","w")
	log_handle.write("This is the master log for Alignment_filter.py v0.1 run at %s\n\nProcessing file %s\n\nThe number of excluded positions due to missing data (threshold was %s):%s\n\nExcluded positions due to excessive gap data (threshold was %s): %s" % (time.asctime(time.localtime()),infile,missing_threshold,missing_number,gap_threshold,gap_number))

def plot_populate (dictionary, tuple_list):
	for i in tuple_list:
		dictionary[i[0]] = i[1]
	return dictionary
	
def missing_data_ploter (infile, excluded_missing, excluded_gaps,sequence_length):
	# Preparing dataset
	position = [x[0] for x in excluded_missing]
	missing_list = [x[1] for x in excluded_missing]
	gap_list = [x[1] for x in excluded_gaps]
	# Plotting
	plt.plot(position, missing_list, "r-", position, gap_list, "b-")
	plot_code = infile.split(".")[0]+".svg"
	plt.savefig(plot_code)
	plt.close()
def dump_fasta (storage, storage_order, infile):
	outfile_code = infile.split(".")[0]+"_missingFilter"
	out_file = open(outfile_code+".fas","w")
	for key in storage_order:
		out_file.write(">"+key.replace(" ","_")+"\n"+storage[key]+"\n")
	out_file.close()

def main ():
	input_file_list = arg.infile
	missing = "".join(arg.threshold_missing)
	gap = "".join(arg.threshold_gap)
	missing_data_number, gaps_number = 0,0
	# Run the filters for each alignment individually
	for infile in input_file_list:
		# The dataset_creator was designed to handle lists of files, hence the list() function
		storage,taxa_order = dataset_creator([infile])
		alignment_length = len(list(storage.values())[0])
		master_storage = dict((sp,"") for sp in taxa_order)
		if alignment_length > 20000:
			subpart_list = data_breaker(storage, alignment_length)		
			for subpart in subpart_list:
				loading(subpart_list.index(subpart),len(subpart_list),"Processing alignment", 50)
				sub_storage = pickle.load( open(subpart, "rb"))
				if subpart == subpart_list[0]:
					filter_sequences, excluded_missing, excluded_gaps, sequence_length, total_missing, total_gap = file_filter (sub_storage, arg.missing, arg.gap, missing, gap, multi="yes", extremities_2="begin")
				elif subpart == subpart_list[-1]:
					filter_sequences, excluded_missing, excluded_gaps, sequence_length, total_missing, total_gap = file_filter (sub_storage, arg.missing, arg.gap, missing, gap, multi="yes", extremities_2="end")
				else:
					filter_sequences, excluded_missing, excluded_gaps, sequence_length, total_missing, total_gap = file_filter (sub_storage, arg.missing, arg.gap, missing, gap, multi="yes", extremities_2="no")
				logger (subpart,excluded_missing,excluded_gaps, missing, gap)
				missing_data_ploter(subpart, total_missing, total_gap, sequence_length)
				for header, sequence in filter_sequences.items():
					master_storage[header] += sequence
			dump_fasta (master_storage, taxa_order, infile)
			[subprocess.Popen(["rm %s" % x],shell=True).wait() for x in subpart_list]
		else:
			if len(input_file_list) == 1:
				filter_sequences, excluded_missing, excluded_gaps, sequence_length, total_missing, total_gap = file_filter (storage, arg.missing, arg.gap, missing, gap)
				missing_data_number += len(excluded_missing)
				gaps_number += len(excluded_gaps)
				logger (infile,excluded_missing,excluded_gaps, missing, gap)
				dump_fasta (filter_sequences, taxa_order, infile)
				missing_data_ploter(infile, total_missing, total_gap, sequence_length)
			elif len(input_file_list) > 1:
				loading(input_file_list.index(infile)+1, len(input_file_list), "Processing alignments (%s)" % infile, 50)
				filter_sequences, excluded_missing, excluded_gaps, sequence_length, total_missing, total_gap = file_filter (storage, arg.missing, arg.gap, missing, gap, multi="yes")
				missing_data_number += len(excluded_missing)
				gaps_number += len(excluded_gaps)
				logger (infile,excluded_missing,excluded_gaps, missing, gap)
				dump_fasta (filter_sequences, taxa_order, infile)
				if arg.plot != None:
					missing_data_ploter(infile, total_missing, total_gap, sequence_length)
	master_log (infile, missing_data_number, gaps_number, missing, gap)
	
main()
