#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#  codeml_parser.py
#  
#  Copyright 2014 Diogo N. Silva
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
#  

import os
import locale
import vincent
import argparse
import subprocess
from collections import OrderedDict
import statsmodels.sandbox.stats.multicomp as multi_correction


parser = argparse.ArgumentParser(description="codeml_parser is a python 2 script developed to parse the output of "
											"multiple branch site model results from PAML, performs chi-square test "
											"of significance and produces informative plots and tables")

parser.add_argument("-in", dest="folder_list", nargs="+", required=True, help="The directories containing the output "
																			"of the branch site model")
# TODO: Provide choices in the clades option
parser.add_argument("--clade", "-c", dest="clade", nargs="*", help="For clade specific statistics, provide the taxa "
																"members"
															"of such clade separated by whitespace")
parser.add_argument("-o", dest="output_file", help="Provide the name for the output file")
parser.add_argument("-hist", dest="histogram", nargs="*", choices=["all_selected"], help="Construct histograms from the"
																						"data and according to a"
																						"run mode of your choice")
parser.add_argument("-p", dest="plot_options", nargs="*", choices=["a", "b"], help="Various plot options: 1 - Donut "
																				"plot"
																				"with the distribution of sites"
																				"across several selection classes")
parser.add_argument("-w", dest="write_alignments", action="store_const", const=True, help="Use this option if you want"
																						"to write the gap-free "
																						"alignments used by codeml")

arg = parser.parse_args()


class PamlPair ():
	""" This class contains the information for a given pair of Null/Alternative analyses from the PAML branch site
	model. It contains all the relevant information about that single pair. If the gene's pair was identified under
	positive selection by LRT, then several attributes related to positive selection findings will be available """

	def __init__(self, folder):
		""" To instantiate a PamlPair object, a folder containing the Null and Alternative sub folders """

		# Initialize p-value attribute
		self.pvalue = None
		self.fdr_value = None

		# Initialize likelihood attributes
		self.alternative_lnL = None
		self.null_lnL = None

		# A status attribute to escape possible errors
		self.status = True

		# Initialize aa filter attributes
		self.conserved_aa = None
		self.unique_aa = None
		self.diverse_aa = None
		self.all_clade_unique = None
		self.mostly_unique = None
		self.mostly_diverse = None
		self.mostly_conserved = None

		# Set folder variable
		self.folder = folder

		self.__parse_main_alternative__("y_selection")
		self.__parse_main_null__("y_selection")

	def set_fdr(self, fdr):
		""" Update the pvalue """

		self.fdr_value = fdr

	def __parse_main_alternative__(self, file_suffix):
		""" This function parses the main output file of the alternative model, setting a number of object attributes.
		For this to work, the output files of the alternative hypothesis must be inside a folder name "Alternative"

		A list of attributes follows:
		- self.gene_length : (int) The length of the gene (gene length is 0 is alignment is non-existent);
		- self.alternative_lnL : (str) The likelihood of the alternative model;
		- self.selected_aa: (list) A list of tuples, with each element containing the position of the aa and the PP
		associated with it;
		- self.conserved_prop : (tuple) A tuple containing the proportion, background w value, foreground w value for
		the conserved w class;
		- self.neutral_prop : (tuple) A tuple containing the proportion, background w value, foreground w value for
		the neutral w class;
		- self.high_selection_prop : (tuple) A tuple containing the proportion, background w value, foreground w
		value for the w class with positive selection on the foreground branch and negative selection on the background
		branch;
		- self.selection_prop : (tuple) A tuple containing the proportion, background w value, foreground w value for
		the w class with positive selection on the foreground branch and neutral selection on the background branch;
		"""

		# Assuming that the output files of the alternative hypothesis are inside a folder named "Alternative",
		# this will find the codeml output file based on a specific substring (file_suffix)
		folder_contents = os.listdir(self.folder + "/Alternative")
		file_path = self.folder + "/Alternative/" + "".join([x for x in folder_contents if file_suffix in x])

		# If the output file does not exist, set the status for this Pair object to false and return
		if file_path == self.folder + "/Alternative/":
			self.status = False
			return 1

		# Opening output file
		try:
			file_handle = open(file_path)
		except:
			print("Cannot open file %s" % file_path)
			raise SystemExit

		# Creating variable that will store the alignment
		self.alignment = OrderedDict()
		counter, alignment_counter, fall_back_alignment_counter = 0, 0, 0
		self.selected_aa = []

		for line in file_handle:
			if line.strip().startswith("seed used = "):
				fall_back_alignment_counter = 1
				# Getting preliminary gene length. If the alignment is not reduced due to gaps, this gene_length
				# variable will be used. Otherwise it will be replaced by the gap free alignment length
				gene_length_str = next(file_handle)

				if gene_length_str.strip() != "":
					self.gene_length = gene_length_str.strip().split()[1]

			if line.strip().startswith("Before deleting alignment gaps"):
				fall_back_alignment_counter = 0

			# Getting the gene length
			if line.strip().startswith("After deleting gaps."):
				self.gene_length = int("".join([x for x in line if x.isdigit()]))

			# Getting the likelihood of the model
			if line.strip().startswith("lnL"):
				self.alternative_lnL = float(line.split(":")[-1].split()[0])

			# Getting the proportion of sites for each class
			if line.strip().startswith("proportion"):
				proportion_vals = line.split()[1:]
				proportion_vals = [float(x) for x in proportion_vals]

			# Getting the background w values
			if line.strip().startswith("background w"):
				background_w = line.split()[2:]
				background_w = [float(x) for x in background_w]

			# Getting the foreground w values
			if line.strip().startswith("foreground w"):
				foreground_w = line.split()[2:]
				foreground_w = [float(x) for x in foreground_w]

			if line.strip().startswith("Bayes Empirical Bayes (BEB)"):
				counter = 1
				next(file_handle)

			# Reset alignment counter
			if counter == 1 and line.strip().startswith("The grid"):
				counter = 0

			# Reset fall_back_alignment counter
			if fall_back_alignment_counter == 1 and line.strip().startswith("Printing out site pattern counts"):
				fall_back_alignment_counter = 0

			if counter == 1 and line.strip() != "":
				aa = line.split()
				if "*" in aa[-1]:
					self.selected_aa.append((aa[0], aa[-1]))  # aa[0]: position of the aa; aa[-1] the PP value of
					# positive selection

			if line.strip().startswith("After deleting gaps"):

				# In case the cleaned alignment is empty
				if " 0 sites" not in line:
					alignment_counter = 1
				else:
					self.gene_length = 0

			if line.strip().startswith("Printing out site pattern counts"):
				alignment_counter = 0

			if len(line.split("_")) != 1 and alignment_counter == 1:

				fields = line.strip().split()
				species = fields[0]
				sequence = fields[1:]

				self.alignment[species] = sequence

			if len(line.split("_")) != 1 and fall_back_alignment_counter == 1 and self.alignment != {}:

				fields = line.strip().split()
				species = fields[0]
				sequence = fields[1:]

				self.alignment[species] = sequence

		else:

			try:
				# Assigning proportions and w values to object attributes
				self.conserved_prop = (proportion_vals[0], background_w[0], foreground_w[0])  # The list contains
				# [proportion, background_w, foreground_w] for the w class 0
				self.neutral_prop = (proportion_vals[1], background_w[1], foreground_w[1])
				self.high_selection_prop = (proportion_vals[2], background_w[2], foreground_w[2])
				self.selection_prop = (proportion_vals[3], background_w[3], foreground_w[3])
			except NameError:
				pass

		file_handle.close()

	def __parse_main_null__(self, file_suffix):
		""" This function parses the main output file of the null model, setting a number of object attributes. A list
		of attributes follows:
		self.null_lnL : (str) The likelihood of the alternative model;
		"""

		# Assuming that the output files of the null hypothesis are inside a folder named "Null",
		# this will find the codeml output file based on a specific substring (file_suffix)
		folder_contents = os.listdir(self.folder + "/Null")
		file_path = self.folder + "/Null/" + "".join([x for x in folder_contents if file_suffix in x])

		# If the output file does not exist, set the status for this Pair object to false and return
		if file_path == self.folder + "/Null/":
			self.status = False
			return 1

		try:
			file_handle = open(file_path)
		except:
			print("Cannot open file %s" % file_path)
			raise SystemExit

		for line in file_handle:

			if line.strip().startswith("lnL"):
				self.null_lnL = float(line.split(":")[-1].split()[0])

		file_handle.close()

	def likelihood_ratio_test(self):
		""" Conducts a likelihood ratio test of the self.alternative_lnL and self.null_lnL values and returns a float
		with the p-value. It returns the chi2 p-value and also sets it as an object attribute: self.lrt_pvalue """

		# In case one of the likelihood values does not exist in one of the files, set p-value as 1. This indicates a
		#  problem somewhere in the analysis of this pair and therefore I attributed a non significant p-value
		if self.alternative_lnL is None or self.null_lnL is None:
			self.pvalue = 1.0
			return

		lrt = 2 * (self.alternative_lnL - self.null_lnL)

		# If the LRT value is below 0, there is something wrong since a null hypothesis with less parameters should
		# not have a higher likelihood than the alternative hypothesis. This usually means that there was some
		# convergence problems in the alternative test.
		if lrt < 0:
			self.pvalue = 1.0
			return

		# Calculating chi-square test
		proc = subprocess.Popen(["chi2 1 %s" % lrt], shell=True, stdout=subprocess.PIPE)
		chi2_output = proc.stdout.read()

		encoding = locale.getdefaultlocale()[1]
		chi2_output = chi2_output.decode(encoding)

		self.pvalue = float(chi2_output.split("=")[2])

	def filter_aa(self, clade, set_aa_columns=None):
		""" This function returns a number of selected amino acid filters, such as conserved, unique or diversifying
		amino acids. A clade of species must be provided and the number of unique and diversifying selected amino
		acids to that clade will be returned. It always returns a tuple. If there are positively selected sites,
		the tuple has three elements, otherwise it contains only an "NA" string
		The set_aa_columns option can be set to True, so that the function also sets, for each site class, a list
		containing the entire alignment column. This allows nucleotide/codon trends in the data set but decreases
		speed."""

		def detect_unique_aa(alignment, taxa_list):
			""" Returns the number of unique and exclusive amino acids of a given clade, irrespective of the presence
			of selection """
			aa_count = 0

			# Alignment sanity check
			if len(alignment) > 0 and len(list(alignment.values())[0]) > 0:

				for i in range(len(list(alignment.values())[0])):
					taxa_specific_aa = [codon_table[char[i]] for sp, char in alignment.items() if sp in taxa_list]
					other_taxa_aa = [codon_table[char[i]] for sp, char in alignment.items() if sp not in taxa_list]

					if len(set(taxa_specific_aa)) == 1 and taxa_specific_aa[0] not in set(other_taxa_aa):
						aa_count += 1

			return aa_count

		def detect_conserved_aa(alignment, taxa_list):
			""" Returns a list of sets containing all conserved and mostly conserved codons in the clade group,
			excluding the selected positions """

			all_conserved = []
			all_mostly_conserved = []

			# Alignment sanity check
			if len(alignment) > 0 and len(list(alignment.values())[0]) > 0:

				for i in range(len(list(alignment.values())[0])):

					if i not in selected_positions:
						column = [codon_table[char[i]] for sp, char in alignment.items()]
						most_common_aa = [x for x in set(column) if all([column.count(x) >= column.count(y) for y in set(column)])]
						most_common_aa_frequency = float(column.count(most_common_aa[0])) / float(len(column))

						if len(set(column)) == 1:

							clade_codons = [char[i] for sp, char in alignment.items() if sp in taxa_list]
							other_codons = [char[i] for sp, char in alignment.items() if sp not in taxa_list]
							all_conserved.append((clade_codons, other_codons))

						if most_common_aa_frequency > 0.70:

							clade_codons = [char[i] for sp, char in alignment.items() if sp in taxa_list]
							other_codons = [char[i] for sp, char in alignment.items() if sp not in taxa_list]
							all_mostly_conserved.append((clade_codons, other_codons))

			return all_conserved, all_mostly_conserved

		# List of preset clades
		preset_dic = {"pucciniales_genome": ["Puccinia_triticina", "Melampsora_laricis-populina",
											"Puccinia_graminis"], "pucciniales": ["Puccinia_triticina",
											"Melampsora_laricis-populina", "Puccinia_graminis", "Hemileia_vastatrix"]}

		# Inspecting if a preset clade is to be used
		for preset in preset_dic:
			if [preset] == clade:
				clade = preset_dic[preset]

		codon_table = {
					'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
					'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
					'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
					'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
					'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
					'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
					'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
					'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
					'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
					'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
					'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
					'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
					'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
					'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
					'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
					'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}

		selected_positions = [int(pos[0]) - 1 for pos in self.selected_aa]

		# Check if there are any selected aa in this pair
		if self.selected_aa is not []:

			# Simple numerical attributes of the object
			self.conserved_aa, self.unique_aa, self.diverse_aa, self.all_clade_unique, self.mostly_conserved, \
			self.mostly_unique, self.mostly_diverse = 0, 0, 0, 0, 0, 0, 0

			if set_aa_columns is True:
				# More complex attributes containing lists of tuples
				self.conserved_aa_list = []
				self.unique_aa_list = []
				self.diverse_aa_list = []
				self.mostly_conserved_aa_list = []

			for aminoacid in self.selected_aa:
				position = int(aminoacid[0]) - 1  # The position of the aa in the alignment
				aa_column = [codon_table[char[position]] for char in self.alignment.values()]
				#codon_column = [char[position] for char in self.alignment.values()]
				unique_aa_colum = set(aa_column)
				self.most_common_aa = [x for x in unique_aa_colum if all([aa_column.count(x) >= aa_column.count(y) for y in unique_aa_colum])]

				# Check if there is only one variant
				if len(unique_aa_colum) == 1:
					self.conserved_aa += 1

					if set_aa_columns is True and clade is not None:
						clade_codon_list = [char[position] for sp, char in self.alignment.items() if sp in clade]
						other_codon_list = [char[position] for sp, char in self.alignment.items() if sp not in clade]
						self.conserved_aa_list.append((clade_codon_list, other_codon_list))

					continue

				if clade is not None:

					clade_specific_aa = [codon_table[char[position]] for sp, char in self.alignment.items() if sp in clade]
					other_aa = [codon_table[char[position]] for sp, char in self.alignment.items() if sp not in clade]

					# Counts the number of positively selected sites exclusive and homogeneous to a given clade
					if len(set(clade_specific_aa)) == 1 and clade_specific_aa[0] not in other_aa:
						self.unique_aa += 1
						continue

					# Counts the number of positively selected sites exclusive but with variation within a given clade
					elif [x for x in clade_specific_aa if x in other_aa] is []:
						self.diverse_aa += 1
						continue

					elif self.most_common_aa is not []:
						if len(set(clade_specific_aa)) == 1 and clade_specific_aa[0] not in self.most_common_aa:
							self.mostly_unique += 1
							continue
						if [x for x in clade_specific_aa if x in self.most_common_aa] is []:
							self.mostly_diverse += 1
							continue

				# Check if the most common aa is also mostly prevalent. This will add to the mostly conserved class,
				# which relaxes the conserved_aa variant by allowing a small percentage of sites to mutate. Threshold
				# is at 0.85
				if self.most_common_aa:
					frequency_most_common_aa = float(aa_column.count(self.most_common_aa[0])) / float(len(aa_column))

					if frequency_most_common_aa > 0.70:
						self.mostly_conserved += 1

						if set_aa_columns is True:
							clade_codon_list = [char[position] for sp, char in self.alignment.items() if sp in clade]
							other_codon_list = [char[position] for sp, char in self.alignment.items() if sp not in clade]
							self.mostly_conserved_aa_list.append((clade_codon_list, other_codon_list))

							#clade_aa_list = [codon_table[char[position]] for sp, char in self.alignment.items() if
							# sp in clade]
							#other_aa_list = [codon_table[char[position]] for sp, char in self.alignment.items() if
							# sp not in clade]

						continue

			self.all_clade_unique = detect_unique_aa(self.alignment, clade)

			if set_aa_columns is True:

				self.all_conserved, self.all_mostly_conserved = detect_conserved_aa(self.alignment, clade)

			return self.conserved_aa, self.unique_aa, self.diverse_aa, self.mostly_conserved, self.mostly_unique, \
				self.mostly_diverse

		else:
			return "NA"


class PamlPairSet ():
	""" This class will contain a variable number of PamlPair objects and will provide a number of methods for their
	analyses """

	def __init__(self, folder_list):
		""" The object is initialized with a folder list, each of which contains both the Null and Alternative model
		folders, which will be used to create the PamlPair object. The PamlPair objects will be stored in a dictionary
		with their corresponding folder (gene name) as a key """

		self.paml_pairs = OrderedDict()

		for folder in folder_list:
			print("\rProcessing folder %s out of %s (%s)" % (folder_list.index(folder), len(folder_list), folder)),

			paml_object = PamlPair(folder)

			if paml_object.status is True:

				self.paml_pairs[folder] = PamlPair(folder)

	def pair_list(self):
		""" Returns a list with the PamlPair objects """

		return [val for val in self.paml_pairs.values()]

	def test_selection_suite(self):
		""" wrapper for the basic selection test and FDR correction """

		self.test_selection()
		self.fdr_correction()

	def write_alignments(self):
		""" Writes the actual alignments used in codeml, i.e., without gaps, in Fasta format """

		for gene, pair in self.paml_pairs.items():

			output_handle = open(gene + ".fas", "w")

			for sp, seq in pair.alignment.items():

				output_handle.write(">%s\n%s\n" % (sp, "".join(seq)))

			output_handle.close()

	def get_number_aa(self):
		""" Prints the number of mostly conserved sites for the R,S and L amino acids """

		self.R, self.L, self.S = 0, 0, 0

		for gene, pair in self.paml_pairs.items():

			if pair.fdr_value < 0.05:

				if pair.mostly_conserved:
					try:
						most_common = "".join(pair.most_common_aa).strip()
						if most_common == "S":
							self.S += 1
						if most_common == "R":
							self.R += 1
						if most_common == "L":
							self.L += 1
						if most_common == "W":
							print(gene)
					except:
						continue

		print("S: %s, L: %s, R: %s" % (self.S, self.L, self.R))

	def test_selection(self):
		""" For each PamlPair object in the self.paml_pairs conduct a likelihood ratio test for positive selection """

		for pair in self.paml_pairs.values():
			print("\rProcessing selection tests on file %s of %s (%s)" % (list(self.paml_pairs.values()).index(pair),
																		len(self.paml_pairs.values()), pair.folder)),

			pair.likelihood_ratio_test()

	def get_class_proportion(self):
		""" sets the number of genes that contain a given site class as new attributes """

		self.class_proportions = {"conserved": 0, "mostly_conserved": 0, "unique": 0, "diversifying": 0,
								"mostly_unique": 0, "mostly_diverse": 0}
		selected_genes = 0

		for pair in self.paml_pairs.values():

			if pair.conserved_aa is not None and pair.conserved_aa > 0 and pair.fdr_value < 0.05:

				self.class_proportions["conserved"] += 1

			if pair.mostly_conserved is not None and pair.mostly_conserved > 0 and pair.fdr_value < 0.05:

				self.class_proportions["mostly_conserved"] += 1

			if pair.unique_aa is not None and pair.unique_aa > 0 and pair.fdr_value < 0.05:

				self.class_proportions["unique"] += 1

			if pair.diverse_aa is not None and pair.diverse_aa > 0 and pair.fdr_value < 0.05:

				self.class_proportions["diversifying"] += 1

			if pair.mostly_unique is not None and pair.mostly_unique > 0 and pair.fdr_value < 0.05:

				self.class_proportions["mostly_unique"] += 1

			if pair.mostly_diverse is not None and pair.mostly_diverse > 0 and pair.fdr_value < 0.05:

				self.class_proportions["mostly_diverse"] += 1

			if pair.fdr_value < 0.05:

				selected_genes += 1

		for key, val in self.class_proportions.items():

			self.class_proportions[key] = float(val) / float(selected_genes)

	def get_gene_class_proportions(self):
		""" For each gene, get the propotion of sites for each site class """

		def gene_pie(storage):
			""" Creates an histogram with the number of genes with the most prevalent class for each site class """

			conserved_count, mostly_conserved_count, unique_count, diversifying_count = 0, 0, 0, 0

			for gene, vals in storage.items():

				maximum_val = max(vals)

				most_prevalent = [i for i, j in enumerate(vals) if j == maximum_val]

				for pos in most_prevalent:

					if pos == 0:

						conserved_count += 1

					if pos == 1:

						unique_count += 1

					if pos == 2:

						diversifying_count += 1

			pie_data = {"Conserved": conserved_count, "Unique": unique_count, "Diversifying": diversifying_count}

			pie_chart = vincent.Pie(pie_data)
			pie_chart.legend("Site classes")
			pie_chart.to_json("Gene_site_class_distribution.json")

		gene_storage = OrderedDict()  # Order of the list elements [conserved, mostly conserved, unique, diversifying]

		for gene, pair in self.paml_pairs.items():

			if pair.fdr_value < 0.05:

				if pair.selected_aa:

					number_selected_aa = float(len(pair.selected_aa))

					conserved_proportion = float(pair.conserved_aa) / number_selected_aa
					mostly_conserved_proportion = float(pair.mostly_conserved) / number_selected_aa
					unique_proportion = float(pair.unique_aa) / number_selected_aa
					diversifying_proportion = float(pair.diverse_aa) / number_selected_aa

					gene_storage[gene] = [conserved_proportion + mostly_conserved_proportion, unique_proportion,
										diversifying_proportion]

		output_file = open("Gene_class_proportion.csv", "w")

		output_file.write("Gene; Conserved; Unique; Diversifying\n")

		for gene, vals in gene_storage.items():

			output_file.write("%s; %s; %s; %s\n" % (gene, vals[0], vals[1], vals[2]))

		output_file.close()

		gene_pie(gene_storage)

	def fdr_correction(self, alpha=0.05):
		""" Applies a False Discovery Rate correction to the p-values of the PamlPair objects """

		pvalue_dict = OrderedDict()

		for gene, pair in self.paml_pairs.items():

			if pair. pvalue is not None:

				pvalue_dict[gene] = pair.pvalue

		pvalue_list = [pval for pval in pvalue_dict.values()]

		fdr_bool_list, fdr_pvalue_list, alpha_S, alpha_B = multi_correction.multipletests(pvalue_list, alpha=alpha,
																						method="fdr_bh")

		# Updating PamlPairs with corrected p-value
		for gene, fdr_val in zip(pvalue_dict, fdr_pvalue_list):

			self.paml_pairs[gene].set_fdr(fdr_val)

	def filter_aa(self, clade, set_aa_columns=None):
		""" A wrapper that applies the filter_aa method of the PamlPair to every pair """

		for gene, pair in self.paml_pairs.items():

			pair.filter_aa(clade, set_aa_columns=set_aa_columns)

	def site_histogram(self, runmode):
		""" This method will produce a variety of histograms depending on the run modes.The supported run modes follow:
			all_selected - Histogram of the distribution of selected sites per gene """

		if "all_selected" in runmode:

			# Retrieving a list containing the selected sites for each PamlPair with evidence of selection (fdr < 0.05)
			selected_site_list = []

			for gene, pair in self.paml_pairs.items():
				if pair.fdr_value < 0.05:
					selected_site_list.append(len(pair.selected_aa))

			# Setting Bar plot list
			bar_list = []
			for i in range(max(selected_site_list) + 1):
				bar_list.append(selected_site_list.count(i))

			# Constructing plot
			site_hist = vincent.Bar(bar_list)
			site_hist.axis_titles(x="Number of selected sites (FDR < 0.05)", y="Frequency")
			site_hist.to_json("Site_histogram.json")

	def donut_selected_classes (self):
		""" Creating of a donut-style plot with the number of positively selected sites for each of the following
		classes: conserved, unique and diversifying """

		selected_aa, conserved_aa, unique_aa, diverse_aa, mostly_conserved, mostly_unique, mostly_diverse = 0, 0, 0, \
																											0, 0, 0, 0

		for gene, pair in self.paml_pairs.items():

			if pair.fdr_value < 0.05:

				selected_aa += len(pair.selected_aa)

				if pair.conserved_aa is not None:
					conserved_aa += pair.conserved_aa

				if pair.unique_aa is not None:
					unique_aa += pair.unique_aa

				if pair.diverse_aa is not None:
					diverse_aa += pair.diverse_aa

				if pair.mostly_conserved is not None:
					mostly_conserved += pair.mostly_conserved

				try:
					mostly_unique += pair.mostly_unique
				except:
					pass

				try:
					mostly_diverse += pair.mostly_diverse
				except:
					pass

		no_annotation_sites = selected_aa - (conserved_aa + unique_aa + diverse_aa + mostly_conserved + mostly_unique +
											mostly_diverse)

		data_series = OrderedDict([("Conserved sites (%s)" % conserved_aa, conserved_aa), ("Mostly conversed sites "
					"(""%s)" % mostly_conserved, mostly_conserved), ("Unique sites (%s)" % unique_aa, unique_aa),
					("Diverse sites (%s)" % diverse_aa, diverse_aa), ("No annotation (%s)" % no_annotation_sites,
					no_annotation_sites), ("Mostly unique (%s)" % mostly_unique, mostly_unique), ("Mostly diverse ("
					"%s)" % mostly_diverse, mostly_diverse)])

		#{"Conserved sites (%s)" % (conserved_aa): conserved_aa, "Unique sites (%s)" % (unique_aa): unique_aa, "Diverse
		#  sites (%s)" %(diverse_aa):diverse_aa, "No annotation (%s)" % (no_annotation_sites): no_annotation_sites}

		class_chart = vincent.Pie(data_series, inner_radius=150)
		class_chart.colors(brew="Set2")
		class_chart.legend("Selection classes")

		class_chart.to_json("Class_distribution.json")

	def check_trend_conserve(self):
		""" This method can only be applied after the filter_aa method. It parses all codon columns of the conserved
		and mostly conserved sites to check for a trend in codon/nucleotide usage """

		import pandas as pd

		codon_table = {
					'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
					'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
					'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
					'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
					'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
					'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
					'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
					'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
					'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
					'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
					'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
					'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
					'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
					'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
					'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
					'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}

		def check_codons(set_list):
			""" Given a list of sets, with the first element of a set containing the clade codons and the second
			element the other codons, this determines codon usage bias. For a given amino acid, it counts the usage of
			the codons for the clade and the other taxa """

			clade_count = {}
			other_count = {}

			for group in set_list:

				if group:
					for codon in group[0][0]:
						if codon_table[codon] in clade_count:
							clade_count[codon_table[codon]].append(codon)
						else:
							clade_count[codon_table[codon]] = [codon]

					for codon in group[0][1]:
						if codon_table[codon] in other_count:
							other_count[codon_table[codon]].append(codon)
						else:
							other_count[codon_table[codon]] = [codon]

			common_aa_list = [aa for aa in clade_count if aa in other_count]

			clade_storage = []
			other_storage = []

			for aa in common_aa_list:

				clade_temp_dic = {}
				other_temp_dic = {}

				clade_list = clade_count[aa]
				other_list = other_count[aa]

				total_list = clade_list + other_list

				for codon in set(total_list):

					clade_temp_dic[codon] = float(clade_list.count(codon)) / float(len(clade_list))
					other_temp_dic[codon] = float(other_list.count(codon)) / float(len(other_list))

				clade_storage.append((aa, clade_temp_dic))
				other_storage.append((aa, other_temp_dic))

			return other_storage, clade_storage

		def check_nucleotides(set_list):
			""" Given a list of sets, with the first element of a set containing the clade codons and the second
			element the other codons, it returns a list with two tuples [(clade nucleotide count (1,2,3,4)),
			(other nucleotide count (1,2,3,4)))] """

			clade_count = {"nucA": 0, "nucT": 0, "nucG": 0, "nucC": 0}
			other_count = {"nucA": 0, "nucT": 0, "nucG": 0, "nucC": 0}

			for group in set_list:
				for codon in group[0][0]:
					for nuc in nucleotides:
						clade_count["nuc%s" % nuc] += codon.count(nuc)

				# count nucleotides for other
				for codon in group[0][1]:
					for nuc in nucleotides:
						other_count["nuc%s" % nuc] += codon.count(nuc)

			# Get nucleotide proportions instead of absolute numbers
			clade_nuc_total = sum(clade_count.values())
			other_nuc_total = sum(other_count.values())

			clade_prop = dict((x, float(y) / float(clade_nuc_total)) for x, y in clade_count.items())
			other_prop = dict((x, float(y) / float(other_nuc_total)) for x, y in other_count.items())

			return clade_prop, other_prop

		def write_plot(clade_storage, other_storage, output_file_name):
			""" wrapper that creates a plot from a storage variable """

			index = ["Clade", "Other"]
			data = [clade_storage, other_storage]
			df = pd.DataFrame(data, index=index)

			grouped_bar = vincent.GroupedBar(df)
			grouped_bar.axis_titles(x='Groups', y='Proportion')
			grouped_bar.legend(title="Legend")
			grouped_bar.to_json(output_file_name)

		nucleotides = ["A", "T", "G", "C"]
		conserved_storage = []
		mostly_conserved_storage = []
		all_conserved_storage = []
		all_mostly_conserved_storage = []

		for gene,pair in self.paml_pairs.items():

			if pair.conserved_aa is not None and pair.conserved_aa > 0:

				conserved_storage.append(pair.conserved_aa_list)

			if pair.mostly_conserved is not None and pair.mostly_conserved > 0:

				mostly_conserved_storage.append(pair.mostly_conserved_aa_list)

			try:
				all_conserved_storage.append(pair.all_conserved)
			except:
				pass

			try:
				all_mostly_conserved_storage.append(pair.all_mostly_conserved)
			except:
				pass

		conserved_counts_clade, conserved_counts_other = check_nucleotides(conserved_storage)
		mostly_conserved_counts_clade, mostly_conserved_counts_other = check_nucleotides(mostly_conserved_storage)

		# Creating plot
		write_plot(conserved_counts_clade, conserved_counts_other, "Conserved_nucleotide_trend.json")
		write_plot(mostly_conserved_counts_clade, mostly_conserved_counts_other, "Mostly_conserved_nucleotide_trend.json")

		conserved_clade_codon, conserve_other_codon = check_codons(conserved_storage)
		mostly_conserved_clade_codon, mostly_conserved_other_codon = check_codons(mostly_conserved_storage)

		all_conserved_clade_codon, all_conserved_other_codon = check_codons(all_conserved_storage)
		all_mostly_conserved_clade_codon, all_mostly_conserved_other_codon = check_codons(all_mostly_conserved_storage)

		for clade, other in zip(conserved_clade_codon, conserve_other_codon):

			write_plot(clade[1], other[1], "Conserved_codon_trend_%s.json" % (clade[0]))

		for clade, other in zip(mostly_conserved_clade_codon, mostly_conserved_other_codon):

			write_plot(clade[1], other[1], "Mostly_conserved_codon_trend_%s.json" % (clade[0]))

		for clade, other in zip(all_conserved_clade_codon, all_conserved_other_codon):

			write_plot(clade[1], other[1], "All_conserved_codon_trend_%s.json" % (clade[0]))

		for clade, other in zip(all_mostly_conserved_clade_codon, all_mostly_conserved_other_codon):

			write_plot(clade[1], other[1], "all_mostly_conserved_codon_trend_%s.json" % (clade[0]))

	def write_table(self, output_file):
		""" Writes the information on the PamlPair objects into a csv table """

		output_handle = open(output_file, "w")
		output_handle.write("Gene;lnL Alternative;lnL Null;p-value;FDR correction;N sites;w class 0;w class 1;w class "
							"2;w class 3;Selected sites; Sites position; Conserved sites;Mostly conserved sites;Unique "
							"sites;Mostly_unique;Diversifying sites;Mostly diverse; All unique sites \n")

		for gene, pair in self.paml_pairs.items():
			print("\rProcessing selection tests on file %s of %s (%s)" % (list(self.paml_pairs.values()).index(pair),
																		len(self.paml_pairs.values()), pair.folder)),

			try:
				gene_length = pair.gene_length
			except:
				gene_length = None

			if pair.fdr_value < 0.05:
				output_handle.write("%s; %s; %s; %s; %s; %s; %s; %s;%s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s\n" % (
										gene,
										pair.alternative_lnL,
										pair.null_lnL,
										pair.pvalue,
										pair.fdr_value,
										gene_length,
										pair.conserved_prop,
										pair.neutral_prop,
										pair.high_selection_prop,
										pair.selection_prop,
										len(pair.selected_aa),
										pair.selected_aa,
										pair.conserved_aa,
										pair.mostly_conserved,
										pair.unique_aa,
										pair.mostly_unique,
										pair.diverse_aa,
										pair.mostly_diverse,
										pair.all_clade_unique))
			else:
				output_handle.write("%s; %s; %s; %s; %s; %s\n" % (gene,
																			pair.alternative_lnL,
																			pair.null_lnL,
																			pair.pvalue,
																			"",
																			gene_length))

		output_handle.close()


if __name__ == "__main__":

	def main():
		# Arguments
		folder_list = arg.folder_list
		clade_list = arg.clade
		output_file = arg.output_file
		plot_options = arg.plot_options
		write_alignments = arg.write_alignments

		# Plotting arguments
		histogram_runmode = arg.histogram

		paml_output = PamlPairSet(folder_list)
		paml_output.test_selection_suite()

		if clade_list is not None:
			if plot_options is not None and "b" in plot_options:
				paml_output.filter_aa(clade_list, set_aa_columns=True)
			else:
				paml_output.filter_aa(clade_list)

		if histogram_runmode is not None and histogram_runmode != []:
			for runmode in histogram_runmode:
				paml_output.site_histogram(runmode)

		if plot_options is not None:

			if "a" in plot_options:
				paml_output.donut_selected_classes()

			if "b" in plot_options:
				paml_output.check_trend_conserve()

		if write_alignments:

			paml_output.write_alignments()

		#paml_output.get_class_proportion()
		paml_output.get_gene_class_proportions()
		paml_output.get_number_aa()
		paml_output.write_table(output_file)

	main()