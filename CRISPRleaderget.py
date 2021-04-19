#!/usr/bin/env python2

from __future__ import print_function
import os
import time
import re
import subprocess
from subprocess import call
from os import path
import sys
import argparse
import distutils.spawn
import textwrap

#######################################################################################################################################################################
############ This script is developed to isolate CRISPR leaders from CRISPRDetect output file in order to identify the motifs within the leader sequences #############
############ Released Version 1 (beta_Version 5) - order is: fetch specified subtype leader sequences, modifying the fasta file to include the CRISPR detehct info, MAFFT alignment #############
############ MAX align, MAFFT Re-align, FastTree, TreeCluster, Calling clade sequences from cluster output												   #############
#######################################################################################################################################################################

orig_stdout = sys.stdout

#################	Get arguments from terminal input	###############

parser = argparse.ArgumentParser(
      prog='python2 CRISPRleaderget.py',
      formatter_class=argparse.RawDescriptionHelpFormatter,
      epilog=textwrap.dedent('''\

      	Author: Murat Buyukyoruk
      	Associated lab: Wiedenheft lab

        CRISPRleaderget Version 1 help:

CRISPRleaderget is a python program developed and tested in Mint Linux operating system. CRISPRleaderget.py should run under any unix based operating system that has a working 'python v2.7' (Not compatible with python3). 

This script is developed to isolate CRISPR leaders from CRISPRDetect output file in order to generate CRISPR leader phylogeny and identification of motifs within the leader sequences. If all the 3rd party dependencies are installed and available in user/system $PATH this automated script should run without any issues.

Please refer to CRISPRDetect for detection of CRISPR arrays, which the output file is required to use as an input for this script. Link to CRISPRDetect: https://github.com/ambarishbiswas/CRISPRDetect_2.2

	Syntax:

	python2 CRISPRleaderget.py -j demo -i demo.fasta -l demo_CRISPRDetect -s I-E

		OR

	python CRISPRleaderget.py -j demo -i demo.fasta -l demo_CRISPRDetect -s All -f 200 -cdc 0.95 -r Yes -cm avg_clade -ct 1.5 -t 50 

CRISPRleaderget dependencies:
	Cd-hit
	MAFFT
	SamTools
	FastTree
	TreeCluster.py
	WebLogo
	Color_tree
	Maxalign

What is this script for?

CRISPRleaderget isolates leader sequences of CRISPR arrays identified by CRISPRDetect v2.4. Fetched leader sequences used to generate NON-redundant list of leader sequences (cdhit). These leader sequences then aligned (MAFFT), gappy/poor alignments were removed (MaxAlign) and realigned (MAFFT). Realignment used to generate phylogeny (FastTree) in order to identified clusters of similar leader sequences (TreeCluster and color_tree[optional]). Each Cluster of leaders processed to identify motifs with Weblogo or additional scripts (not included in this script).

Input Paramaters (REQUIRED):
----------------------------
	-j/--job_name		TEXT			Specify a name to assign a job name. All files generated in this script will include this job name (i.e. bac_arc_db).

	-i/--input		FASTA			Specify a fasta file that has been used to generate CRISPRDetect output file. FASTA file requires headers starting with accession number. (i.e. >NZ_CP006019 [fullname])

	-l/--log		CRISPRDetect_outfile	Specify a CRISPRDetect output file that contains info about the found arrays in searched fasta file. (Please DO NOT use the .fp, .gff outputs. As default CRISPRDetect generates an output file name defined at -o option.)

	-s/--subtype		TEXT			Specify a subtype to be processed. Please use same annotation defined in the CRISPRDetect output file which is usually I-E, I-F, III-A etc. Type 'All' if you want to process all CRISPR leaders.

Parameters [optional]:
----------------------
	-f/--flank		200			This is the default length of leader sequence to be processed. Leader sequence is fetcehd from the left flank of CRISPR array which is considered according to the direction of the CRISPR loci (i.e. right flank if the array is in reverse).

	-cdc/--cdhit_c		0.95			This is the default for the sequence identity threshold used in cdhit to generate non-redundant list of CRISPR leader sequences.

	-r/--repeat		No			Specifies to include the first repeat to the automated process or not. Provides better alignment and phylogeny when the first repeat is included.

	-cm/--cluster_method	avg_clade		This script uses average clade method for identifying the clusters within generated leader phylogeny. However, this option can be changed. Please refer to TreeCluster options for more.

	-ct/--cluster_threshold	0.9		This script uses 0.9 as clustering threshold for identifying the clusters within generated leader phylogeny. However, this option can be changed. Please refer to TreeCluster options for more.

Basic Options:
--------------
	-t/--thread		Threads		MAFFT is able to utilize multiple threads. Specify number of parallel processes CRISPRDleaderget should use. Dedault is 8.

	-h/--help		HELP			Shows this help text and exits the run.
         '''))
parser.add_argument('-j', '--job_name', required=True, type=str, dest='j_n', help='Specify database name (eg. bac_arc_db).\n')
parser.add_argument('-i', '--input', required=True, dest='fasta_name', help='Specify a fasta input.\n')
parser.add_argument('-l', '--log', required=True, dest='filename', help='Specify a CRISPRDetect_outfile input.')
parser.add_argument('-s', '--subtype', required=True, type=str, dest='subtype_name', help='Specify subtype. (All, I-A, I-B, I-F, etc)\n')
parser.add_argument('-f', '--flank', required=False, type=int,dest='flank', default=200, help='Specify length (bp) of flanking regions (Default is 200).\n')
parser.add_argument('-cdc', '--cdhit_c', required=False, type=float, dest='cd_hit', default=0.95, help='Set a threshold value for cd-hit (Default is 0.95).\n')
parser.add_argument('-r', '--repeat', required=False, type=str, dest='get_1st_repeat', default='No', help='Include the first repeat, Yes or No? (Default is No)\n')
parser.add_argument('-cm', '--cluster_method', required=False, type=str, dest='cluster_method', default='avg_clade', help='Specify the clustering method (Default is avg_clade).\n')
parser.add_argument('-ct', '--cluster_threshold', required=False, type=float, dest='cluster_treshold', default=0.9, help='Set a threshold value for clustering (Default is 0.9).\n')
parser.add_argument('-t', '--thread', required=False, type=float, dest='no_of_threads', default=8, help='Specify thread (Default is 8).\n')

results = parser.parse_args()
j_n = results.j_n
fasta_name = results.fasta_name
filename = results.filename
flank = results.flank
cd_hit = results.cd_hit
get_1st_repeat = results.get_1st_repeat
cluster_method = results.cluster_method
cluster_treshold = results.cluster_treshold
no_of_threads = int(results.no_of_threads)
subtype_name = results.subtype_name

################## Check Dependencies ##################

if (distutils.spawn.find_executable("cdhit")) == None:
	print("\nERROR: cdhit is not installed. Please make sure cdhit is defined in path.\n")
	sys.exit()

if (distutils.spawn.find_executable("mafft")) == None:
	print("\nERROR: MAFFT is not installed. Please make sure MAFFT is defined in path.\n")
	sys.exit()
	
if (distutils.spawn.find_executable("samtools")) == None:
	print("\nERROR: samtools is not installed. Please make sure samtools is defined in path.\n")
	sys.exit()

if (distutils.spawn.find_executable("fasttree")) == None:
	print("\nERROR: fasttree is not installed. Please make sure fasttree is defined in path.\n")
	sys.exit()

if (distutils.spawn.find_executable("TreeCluster.py")) == None:
	print("\nERROR: TreeCluster.py is not installed. Please make sure TreeCluster.py is defined in path.\n")
	sys.exit()

if (distutils.spawn.find_executable("weblogo")) == None:
	print("\nERROR: weblogo is not installed. Please make sure weblogo is defined in path.\n")
	sys.exit()

if (distutils.spawn.find_executable("color_tree")) == None:
	print("\nERROR: color_tree is not installed. Please make sure color_tree is defined in path.\n")
	sys.exit()

proc = subprocess.Popen("locate -br '^maxalign.pl$'", shell=True, stdout=subprocess.PIPE, )
path_2_maxalign = (proc.communicate()[0].split('\n')[0])

if path_2_maxalign != "":
	pass

else:
	print("\nERROR: maxalign.pl is not installed. Please make sure maxalign.pl is installed and accessible.\n")
	sys.exit()

print("\nAll dependencies are found.\n")

##################	Identify Job name based on the defined parameters	###############

if (get_1st_repeat == "No"):
	job_name = j_n + '_' + subtype_name + '_'+ str(flank) + 'bp'
elif (get_1st_repeat == "Yes"):
	job_name = j_n + '_'+ subtype_name + '_'+str(flank) + 'bp_repeat'

##################	Identify variables	###############

access_num = []
full_name = []
lower_bound = []
higher_bound = []
new_high_bound =[]
sub_type = []
arr_orientation = [] # REVERSE = 0, FORWARD=1, UNCONFIRMED=2
line_nums = []
total_rep = []
acc = []
subtype = []
new_acc =[]
new_acc_rc =[]
new_acc_fw = []
condition = []
sequence = []
line_nums_sub = []
seq_line = []
repeat = []
comp_name = []
name = []
clade = []
ori_count = 0
new_count = 0
new_new_count = 0
new_new_new_count = 0
count = 0
some_count = 0
more_count = 0

isim = []
clade_num = []
final_countdown = 0


################	Fetch sequences from CRISPRDetect outfile	###############

os.system('echo "\nJob name: "' + job_name)
os.system('echo "\n>>>Getting ' + subtype_name + ' CRISPR leader (1st Repeat: ' + get_1st_repeat + ') sequences from CRISPRDetect outfile<<<\n"')

class Params1:

	def __init__(self, accession_number, lower_bound, higher_bound, arr_orientation, sub_type, new_high_bound, full_name):
		self.accession_number = accession_number
		self.lower_bound = lower_bound
		self.higher_bound = higher_bound
		self.arr_orientation = arr_orientation
		self.sub_type = sub_type
		self.new_high_bound = new_high_bound
		self.full_name = full_name

	def getParam(self, filename):
		global count
		with open(filename,'rU') as file:
			for line in file:
				count += 1
				arr = line.split()
				if (len(arr) != 0):
					if (arr[0] == "Array"):
						pos_complete = arr[2]
						pos_arr = pos_complete.split('-')
						pos_low_bound = int(pos_arr[0])
						pos_high_bound = int(pos_arr[1])
						self.lower_bound.append(pos_low_bound)
						self.higher_bound.append(pos_high_bound)

					if (line[0] == ">" ):
						arr = line.split("\t")
						full_name = arr[0]
						self.full_name.append(full_name)

					if (line[0] == ">" ):
						line = line[:1] + " " + line[1:]
						arr = line.split()
						acc_num = arr[1]
						self.accession_number.append(acc_num)
						if (arr[-1] == "Forward"):
							self.arr_orientation.append(1)
						elif (arr[-1] == "Reverse"): 
							self.arr_orientation.append(0)
						elif (arr[-1] == "Unconfirmed"): 
							self.arr_orientation.append(2)	
						else:
							print ("ERROR in arr_orientation: Expected 0 - 1, got: " + self.arr_orientation)

					if (line[0] == "="):
						line_nums.append(count)
						if (len(line_nums) % 2 == 0 and len(line_nums) != 0):
							total_rep.append((line_nums[-1] - line_nums[-2]) - 1)

					if (line[0:7] == "# Array"):
						sub = arr [4]
						self.sub_type.append(sub)

					if(line[0] == " "): #FORWARD
						arrr = line.split()
						non_decimal = re.compile(r'[^\d.]+')
						arr[0] = non_decimal.sub('',arrr[0])
						if (len(arr[0]) != 0):
							plus = int(line[0:10])
							if(pos_low_bound == plus):
								new_high_bound = int(arr[1])
								self.new_high_bound.append(new_high_bound)
							elif(pos_low_bound == plus+1):
								new_high_bound = int(arr[1])
								self.new_high_bound.append(new_high_bound)

	def retParam(self):
		self.accession_number, self.lower_bound, self.higher_bound, self.arr_orientation, self.sub_type, self.new_high_bound, self.full_name

def runCommand(acc_num, lower_bound, higher_bound, arr_orientation, tot_rep, sub_type, new_high_bound):
	k = 0
	os.system('mkdir temp')
	for (a_n, l_b, h_b, a_o, t_r, s_t, n_h) in zip(acc_num, lower_bound, higher_bound, arr_orientation, tot_rep, sub_type, new_high_bound):
		if (get_1st_repeat == "No"):
			if (s_t == subtype_name):
				if (a_o == 1):
					if (l_b > flank):
						cmd1 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b-flank) + '-' + str(l_b)
						os.system( cmd1 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
					if (l_b < flank):
						cmd2 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b)
						os.system( cmd2 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
				if (a_o == 0):
					if (h_b > flank):
						cmd3 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(l_b) + '-' + str(l_b+flank)
						os.system( cmd3 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
					if (h_b < flank):
						cmd4 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(l_b) + '-' + str(l_b+flank)
						os.system( cmd4 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
				if (a_o == 2):
					if (l_b > flank):
						cmd5 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b-flank) + '-' + str(l_b)
						os.system( cmd5 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
						cmd6 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(h_b) + '-' + str(h_b+flank)
						os.system( cmd6 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
					if (l_b < flank):
						cmd7 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b)
						os.system( cmd7 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
						cmd8 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(h_b) + '-' + str(h_b+flank)
						os.system( cmd8 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
		if (get_1st_repeat == "Yes"):
			if (s_t == subtype_name):
				if (a_o == 1):
					if (l_b > flank):
						cmd9 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b-flank) + '-' + str(l_b+n_h)
						os.system( cmd9 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
					if (l_b < flank):
						cmd10 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b+n_h)
						os.system( cmd10 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
				if (a_o == 0):
					if (h_b > flank):
						cmd11 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(l_b-n_h) + '-' + str(l_b+flank)
						os.system( cmd11 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
					if (h_b < flank):
						cmd12 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(l_b-n_h) + '-' + str(l_b+flank)
						os.system( cmd12 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
				if (a_o == 2):
					if (l_b > flank):
						cmd13 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b-flank) + '-' + str(l_b+n_h)
						os.system( cmd13 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
						cmd14 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(h_b-n_h) + '-' + str(h_b+flank)
						os.system( cmd14 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
					if (l_b < flank):
						cmd15 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b+n_h)
						os.system( cmd15 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
						cmd16 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(h_b-n_h) + '-' + str(h_b+flank)
						os.system( cmd16 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
		if (subtype_name == 'All'):
			if (get_1st_repeat == "No"):
				if (a_o == 1):
					if (l_b > flank):
						cmd17 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b-flank) + '-' + str(l_b)
						os.system( cmd17 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
					if (l_b < flank):
						cmd18 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b)
						os.system( cmd18 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
				if (a_o == 0):
					if (h_b > flank):
						cmd19 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(l_b) + '-' + str(l_b+flank)
						os.system( cmd19 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
					if (h_b < flank):
						cmd20 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(l_b) + '-' + str(l_b+flank)
						os.system( cmd20 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
				if (a_o == 2):
					if (l_b > flank):
						cmd21 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b-flank) + '-' + str(l_b)
						os.system( cmd21 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
						cmd22 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(h_b) + '-' + str(h_b+flank)
						os.system( cmd22 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
					if (l_b < flank):
						cmd23 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b)
						os.system( cmd23 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
						cmd24 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(h_b) + '-' + str(h_b+flank)
						os.system( cmd24 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
			if (get_1st_repeat == "Yes"):
				if (a_o == 1):
					if (l_b > flank):
						cmd25 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b-flank) + '-' + str(l_b+n_h)
						os.system( cmd25 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
					if (l_b < flank):
						cmd26 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b+n_h)
						os.system( cmd26 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
				if (a_o == 0):
					if (h_b > flank):
						cmd27 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(l_b-n_h) + '-' + str(l_b+flank)
						os.system( cmd27 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
					if (h_b < flank):
						cmd28 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(l_b-n_h) + '-' + str(l_b+flank)
						os.system( cmd28 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
				if (a_o == 2):
					if (l_b > flank):
						cmd29 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b-flank) + '-' + str(l_b+n_h)
						os.system( cmd29 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
						cmd30 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(h_b-n_h) + '-' + str(h_b+flank)
						os.system( cmd30 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
					if (l_b < flank):
						cmd31 = 'samtools faidx ' + fasta_name + ' ' + a_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b+n_h)
						os.system( cmd31 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
						cmd32 = 'samtools faidx -i ' + fasta_name + ' ' + a_n + ':' + str(h_b-n_h) + '-' + str(h_b+flank)
						os.system( cmd32 + '> temp/temp_fasta_file' + str(k) + '_.seq')
						k += 1
	os.system('cat temp/*.seq >' + job_name + '.fasta')
	os.system('rm temp/*.seq')

def out_log(acc_num, lower_bound, higher_bound, arr_orientation, tot_rep, sub_type, new_high_bound, full_name):
	i = 0
	for (a_n, l_b, h_b, a_o, t_r, s_t, n_h, f_n) in zip(acc_num, lower_bound, higher_bound, arr_orientation, tot_rep, sub_type, new_high_bound, full_name):
		if (get_1st_repeat == "No"):
			if (s_t == subtype_name):
				if (a_o == 1):
					if (l_b > flank):
						cmd1 =  f_n + ':' + str(l_b-flank) + '-' + str(l_b) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd1 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
					if (l_b < flank):
						cmd2 =  f_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd2 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
				if (a_o == 0):
					if (h_b > flank):
						cmd3 =  f_n + ':' + str(l_b) + '-' + str(l_b+flank) + '/rc' + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd3 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
					if (h_b < flank):
						cmd4 =  f_n + ':' + str(l_b) + '-' + str(l_b+flank) + '/rc' + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd4 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
				if (a_o == 2):
					if (l_b > flank):
						cmd5 =  f_n + ':' + str(l_b-flank) + '-' + str(l_b) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd5 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
						cmd6 = f_n + ':' + str(h_b) + '-' + str(h_b+flank) + '/rc' + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd6 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
					if (l_b < flank):
						cmd7 = f_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd7 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
						cmd8 = f_n + ':' + str(h_b) + '-' + str(h_b+flank) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd8 + '" > temp/temp_' + str(i) + '_.file')
						i += 1

		if (get_1st_repeat == "Yes"):
			if (s_t == subtype_name):
				if (a_o == 1):
					if (l_b > flank):
						cmd9 =  f_n + ':' + str(l_b-flank) + '-' + str(l_b+n_h) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd9 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
					if (l_b < flank):
						cmd10 =  f_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b+n_h) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd10 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
				if (a_o == 0):
					if (h_b > flank):
						cmd11 =  f_n + ':' + str(l_b-n_h) + '-' + str(l_b+flank) + '/rc' + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd11 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
					if (h_b < flank):
						cmd12 =  f_n + ':' + str(l_b-n_h) + '-' + str(l_b+flank) + '/rc' + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd12 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
				if (a_o == 2):
					if (l_b > flank):
						cmd13 =  f_n + ':' + str(l_b-flank) + '-' + str(l_b+n_h) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd13 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
						cmd14 = f_n + ':' + str(h_b-n_h) + '-' + str(h_b+flank) + '/rc' + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd14 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
					if (l_b < flank):
						cmd15 = f_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b+n_h) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd15 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
						cmd16 = f_n + ':' + str(h_b-n_h) + '-' + str(h_b+flank) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd16 + '" > temp/temp_' + str(i) + '_.file')
						i += 1

		if (subtype_name == 'All'):
			if (get_1st_repeat == "No"):
				if (a_o == 1):
					if (l_b > flank):
						cmd17 =  f_n + ':' + str(l_b-flank) + '-' + str(l_b) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd17 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
					if (l_b < flank):
						cmd18 =  f_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd18 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
				if (a_o == 0):
					if (h_b > flank):
						cmd19 =  f_n + ':' + str(l_b) + '-' + str(l_b+flank) + '/rc' + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd19 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
					if (h_b < flank):
						cmd20 =  f_n + ':' + str(l_b) + '-' + str(l_b+flank) + '/rc' + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd20 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
				if (a_o == 2):
					if (l_b > flank):
						cmd21 =  f_n + ':' + str(l_b-flank) + '-' + str(l_b) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd21 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
						cmd22 = f_n + ':' + str(h_b) + '-' + str(h_b+flank) + '/rc' + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd22 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
					if (l_b < flank):
						cmd23 = f_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd23 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
						cmd24 = f_n + ':' + str(h_b) + '-' + str(h_b+flank) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd24 + '" > temp/temp_' + str(i) + '_.file')
						i += 1

			if (get_1st_repeat == "Yes"):
				if (a_o == 1):
					if (l_b > flank):
						cmd25 =  f_n + ':' + str(l_b-flank) + '-' + str(l_b+n_h) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd25 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
					if (l_b < flank):
						cmd26 =  f_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b+n_h) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd26 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
				if (a_o == 0):
					if (h_b > flank):
						cmd27 =  f_n + ':' + str(l_b-n_h) + '-' + str(l_b+flank) + '/rc' + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd27 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
					if (h_b < flank):
						cmd28 =  f_n + ':' + str(l_b-n_h) + '-' + str(l_b+flank) + '/rc' + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd28 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
				if (a_o == 2):
					if (l_b > flank):
						cmd29 =  f_n + ':' + str(l_b-flank) + '-' + str(l_b+n_h) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd29 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
						cmd30 = f_n + ':' + str(h_b-n_h) + '-' + str(h_b+flank) + '/rc' + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd30 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
					if (l_b < flank):
						cmd31 = f_n + ':' + str(l_b - l_b + 1) + '-' + str(l_b+n_h) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd31 + '" > temp/temp_' + str(i) + '_.file')
						i += 1
						cmd32 = f_n + ':' + str(h_b-n_h) + '-' + str(h_b+flank) + '\t| Repeat:' + str(t_r) + '\t| Subtype:' + (s_t)
						os.system( 'echo "' + cmd32 + '" > temp/temp_' + str(i) + '_.file')
						i += 1

	os.system('cat temp/*.file >' + job_name + '_name.txt')
	os.system('rm temp/*.file')
	os.system('rm *fai')
	os.system('rm -r temp')
	os.system('echo ">>" >> ' + job_name+ '.fasta')

##################	Grep info from CRISPRDetect and concatanate it with sequence names	###############

log = job_name + '_name.txt'
raw = job_name + '.fasta'

class Params2:

	def __init__(self, acc, subtype, new_acc_rc, new_acc_fw, new_acc , condition, sequence, seq_line, repeat, comp_name):
		self.acc = acc
		self.subtype = subtype
		self.new_acc_fw = new_acc_fw
		self.new_acc_rc = new_acc_rc
		self.new_acc = new_acc
		self.condition = condition
		self.sequence = sequence
		self.seq_line = seq_line
		self.repeat = repeat
		self.comp_name = comp_name


	def getParam(self, log):
		global ori_count
		with open(log,'rU') as file:
			for line in file:
				ori_count += 1
				if (line[0] == ">" ):
					line = line[:1] + " " + line[1:]
					arr = line.split("|")
					n = arr[0].split('> ')
					m = n[1].split("\t")
					comp_name = m[0]
					self.comp_name.append(comp_name)
					brr = arr[0].split()
					acc = brr[1]
					self.acc.append(acc)
					info = arr[2].split(':')
					a = info[1].split('\n')
					subtype = a[0]
					self.subtype.append(subtype)
					rep = arr[1].split(':')
					b = rep[1].split('\t')
					repeat = b[0]
					self.repeat.append(repeat)
					
	def newParam(self, fasta_name):
		global new_count
		with open(fasta_name,'rU') as file:
			for line in file:
				new_count += 1
				if (line[0] == ">" ): 
					line = line[:1] + " " + line[1:]
					brr = line.split()
					new_acc = brr[1]
					self.new_acc.append(new_acc)

					if '/' in new_acc:
						crr = new_acc.split('/')
						new_acc_rc = crr[0]
						self.new_acc_rc.append(new_acc_rc)
					else:
						new_acc_fw = brr[1]
						self.new_acc_fw.append(new_acc_fw)
				
				if (line[0] == ">"):
					line_nums_sub.append(new_count)
					if (len(line_nums_sub) != 0 and len(line_nums_sub) != 1):
						seq_line = line_nums_sub[-1] - line_nums_sub[-2] - 1
						self.seq_line.append(seq_line)

				if (line[0] != ">"):
					sequence = line
					self.sequence.append(sequence)

	def retParam(self):
		self.acc, self.subtype, self.new_acc_fw, self.new_acc_rc, self.condition, self.sequence, self.new_acc, self.seq_line, self.repeat, self.comp_name

def write(acc, subtype, new_acc_rc, new_acc_fw, condition, sequence, seq_line, repeat, comp_name):
	f = open(job_name + "_info.fasta", 'w')
	sys.stdout = f
	j = 0

	for i in range(len(acc)):
		if (acc[i] == new_acc[i]):
			print (">",comp_name[i]," Subtype ",subtype[i]," Repeat ",repeat[i],sep="")
			print(*sequence[j:(j+seq_line[i])],sep="")
			j = j+seq_line[i]

		else:
			print (">",comp_name[i]," Subtype ",subtype[i]," Repeat ",repeat[i],sep="")
			print(*sequence[j:(j+seq_line[i])],sep="")
			j = j+seq_line[i]	
	sys.stdout = orig_stdout

def process1():

##################	Preparing fasta files for downstream processes: Removing special characters from names	###############

	os.system("""sed -e 's/:/ /g' -e 's/,/ /g' -e "s/'//g" -e 's/(/ /g' -e 's/)/ /g' """ + job_name + '_info.fasta > ' + job_name + '_bcdhit.fasta')
	time.sleep(1)
	os.system("echo 'No of " + subtype_name +" CRISPR leader sequences BEFORE cd-hit:' & grep -o '>' " + job_name + "_bcdhit.fasta | wc -l" )

##################	cd-hit of sequences to remove identical sequences	###############

	os.system('cd-hit -i ' + job_name + '_bcdhit.fasta -o ' + job_name + '_cdhit_' + str(cd_hit) + '.fasta -c ' + str(cd_hit) + '&>/dev/null')
	time.sleep(1)
	os.system("echo 'No of CRISPR leader sequences AFTER cd-hit " + str(cd_hit) +":' & grep -o '>' " + job_name + '_cdhit_' + str(cd_hit) + '.fasta | wc -l')

##################	Multiple sequence alignment with MAFFT	###############

	os.system('echo "\n>>>ALIGNING SEQUENCES<<<\n"')

	os.system('echo "Settings: --maxiterate 1000 --genafpair "')

	os.system('mafft --maxiterate 1000 --threadtb '+ str(no_of_threads) +'  --threadit '+ str(no_of_threads) +' --genafpair --thread ' + str(no_of_threads) + ' ' + job_name + '_cdhit_' + str(cd_hit) + '.fasta > ' + job_name + '_aligned.fasta')
	time.sleep(1)
##################	Max-alignment of aligned sequences	###############

 	os.system('echo "\n>>>MAX_ALIGNING SEQUENCES<<<\n"')

 	os.system('perl ' + path_2_maxalign + ' -d -v=1 -o=5 ' + job_name + '_aligned.fasta')

 	os.system('mv heuristic.fsa ' + job_name + '_max_aligned.fasta')
 	time.sleep(1)
 	os.system("echo 'No of CRISPR leader sequences AFTER maxalign:' & grep -o '>' " + job_name + '_max_aligned.fasta | wc -l')

##################	Re-alignment of max-aligned sequences	###############

	os.system('echo "\n>>>Re_ALIGNING SEQUENCES<<<\n"')

	os.system('mafft --maxiterate 1000 --threadtb '+ str(no_of_threads) +'  --threadit '+ str(no_of_threads) +' --genafpair --thread ' + str(no_of_threads) + ' ' + job_name + '_max_aligned.fasta > ' + job_name + '_realigned.fasta')

#################	Phylogenetic Tree construction with Fasttree	###############

	os.system('echo "\nCONSTRUCTING TREE...\n"') 

	os.system('fasttree -quote -gamma -spr 4 -mlacc 2 -slownni -nt ' + job_name + '_realigned.fasta > tree.' + job_name)
	time.sleep(1)
	os.system('echo "\n>>>TREE GENERATED<<<\n"')

#################	Defining clades in tree with TreeCluster	###############

	os.system('echo "Defining Clades and coloring phylogenetic tree based on identified clades..." "\t Settings: -t "' +  str(cluster_treshold) + ' -m ' + str(cluster_method))

	os.system('TreeCluster.py -t ' + str(cluster_treshold) + ' -m ' + cluster_method + ' -i tree.' + job_name + ' -o ' + job_name + '_cluster_info_' + str(cluster_treshold) + '.txt')
	time.sleep(1)

#################	Creating colored phylogenetic tree based on the identified clades	###############

	sys.stdout = orig_stdout


	os.system("sed 's/$/i/' " + job_name + '_cluster_info_' + str(cluster_treshold) + '.txt > temp_file1.txt')
	time.sleep(0.5)
	os.system("tail -n +2 temp_file1.txt > temp_file2.txt")
	time.sleep(0.5)
	os.system('rm temp_file1.txt')
	time.sleep(0.5)
	os.system('mv temp_file2.txt color_tree_pattern.txt')

file_name = 'color_tree_pattern.txt'

class Params4:


		

	def __init__(self, isim, clade_num ):
		self.isim = isim
		self.clade_num = clade_num

	def getParam(self, file_name):
		global final_countdown
		with open(file_name,'rU') as file:
			for line in file:
				final_countdown += 1
				if (len(line) != 0):
					if (line[0] == "'"):
						arr = line.split("'")
						isim = arr[1]
						self.isim.append(isim)
						brr = arr[2].split("\t")
						crr = brr[1].split("i\n")
						clade_num = int(crr[0])
						self.clade_num.append(clade_num)

	def retParam(self):
		self.isim, self.clade_num

def color(isim, clade_num):
	
	color = ['some', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080', '#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000', '#808000', '#008000', '#800080', '#008080', '#000080']
	low = int(min(clade_num))
	high = int(max(clade_num))
	m = 0
	for i in range (low, high+1):
		for k in range(len(clade_num)):
			if (i != 0):
				if (i == clade_num[k]):
					os.system("sed -e 's/\t" + str(i) + "i/\t"+str(color[i])+"/g' " +file_name+" > temp1.txt")
					time.sleep(0.1)
					os.system('rm ' + file_name)
					time.sleep(0.1)
					os.system('mv temp1.txt ' + file_name)
					time.sleep(0.1)
	os.system("""sed -e "s/'//g" """ + file_name + ' > temp2.txt')
	os.system('rm ' + file_name)
	time.sleep(0.1)
	os.system('mv temp2.txt ' + file_name)
	os.system('color_tree -bt -p ' + file_name + ' tree.' + job_name + ' > ' + job_name + '.nexus')
	time.sleep(0.1)

def tidy():

	os.system("echo 'No of " + subtype_name +" CRISPR leader sequences in the generated tree:' & grep -o '>' " + job_name + "_realigned.fasta | wc -l" )

################	Organizing generated files before fetching clade sequences from TreeCluster outfile	###############

	os.system('mkdir 2_Leader_seq_' + subtype_name)
	time.sleep(0.1)
	os.system('mv ' + job_name + '.fasta 2_Leader_seq_' + subtype_name)

	os.system('mv *bcdhit.fasta 2_Leader_seq_' + subtype_name)

	os.system('mv ' + job_name + '_info.fasta 2_Leader_seq_' + subtype_name)

	os.system('mv ' + job_name + '_name.txt 2_Leader_seq_' + subtype_name)

	os.system('mkdir 3_cdhit')
	time.sleep(0.1)
	os.system('mv *.clstr 3_cdhit')

	os.system('mkdir 4_alignment')
	time.sleep(0.1)
	os.system('mv *aligned.fasta 4_alignment')

	os.system('mkdir 5_max_align')
	time.sleep(0.1)
	os.system('mv heuristic_exclude_headers.txt 5_max_align')

	os.system('mv heuristic_include_headers.txt 5_max_align')

	os.system('mv optimal* 5_max_align')

	os.system('mv report.txt 5_max_align')

	os.system('mkdir 6_fasttree')
	time.sleep(0.1)
	os.system('mv tree.' + job_name + " 6_fasttree")

	os.system('mv ' + job_name + '.nexus 6_fasttree')

	os.system('mv color_tree_pattern.txt 6_fasttree')

	os.system("sed -e 's/ /*/g' -e 's/(/*/g' -e 's/)/*/g' -e 's/;/*/g' -e 's/=/*/g' " + job_name + '_cdhit_' + str(cd_hit) + '.fasta >  temp_cdhit.fasta')

	os.system("sed -e 's/ /*/g' -e 's/(/*/g' -e 's/)/*/g' -e 's/;/*/g' -e 's/=/*/g' " + job_name + '_cluster_info_' + str(cluster_treshold) + '.txt' + ' >  temp_info.txt')

##################	Fetching sequences from each identified clade	###############

clades = 'temp_info.txt'
ori = 'temp_cdhit.fasta'

class Params3:

	def __init__(self, name, clade ):
		self.name = name
		self.clade = clade

	def getParam(self, clades):
		global new_new_new_count
		with open(clades,'rU') as file:
			for line in file:
				new_new_new_count += 1
				if (len(line) != 0):
					if (line[0] == "'"):
						arr = line.split("'")
						name = arr[1]
						self.name.append(name)
						brr = arr[2].split("\t")
						crr = brr[1].split("\n")
						clade = int(crr[0])
						self.clade.append(clade)
	
	def retParam(self):
		self.name, self.clade

def call_clusters(name, clade):
	low = int(min(clade))
	high = int(max(clade))
	m = 0
	for i in range (low, high+1):
		for k in range(len(clade)):
			if (i == clade[k]):
				cmd = 'samtools faidx ' + ori + ' ' + name[k]
				os.system(cmd + ' > ' + job_name + '_full_clade_n_'+str(m)+'_'+str(i)+'.fasta')
				os.system("sed -e 's/*/ /g' " + job_name + '_full_clade_n_'+str(m)+'_'+str(i)+'.fasta > ' + job_name + '_full_clade_'+str(m)+'_'+str(i)+'.fasta')
				time.sleep(0.01)
				os.system("rm " + job_name + '_full_clade_n_'+str(m)+'_'+str(i)+'.fasta')
				m += 1
				time.sleep(0.01)
	for l in range (low, high+1):
		if (l != 0):
			cmd_cat = 'cat *_' + str(l) + '.fasta > ' + job_name + '_full_'+str(l)+'_clade.fasta'
			os.system(cmd_cat)
			time.sleep(2)
			os.system('rm *_'+ str(l) + '.fasta')

##################	Aligning clades and generating weblogos for each clade	###############

			os.system('mafft --maxiterate 1000 --threadtb '+ str(no_of_threads) +'  --threadit '+ str(no_of_threads) +' --genafpair --thread ' + str(no_of_threads) + ' ' + job_name + '_full_'+str(l)+'_clade.fasta > ' + job_name + '_full_'+str(l)+'_clade_aligned.fasta')
			time.sleep(0.1)
			os.system('weblogo -f ' + job_name + '_full_'+str(l)+'_clade_aligned.fasta -D fasta -o ' + job_name + '_full_'+str(l)+'_clade_aligned.eps -A dna -s large --fontsize 15 -C "#D7191C" G "G" -C "#FDAE61" A "A" -C "#2C7BB6" C "C" -C "#ABD9E9" T "T"')

	time.sleep(2)

	os.system('\necho "Identified "' + str(l) + ' clade/s\n')
	os.system('\necho "Generated individual .fasta files, aligned sequences and created weblogos for each clade\n"')

def process2():

##################	Organizing files generated by clade detection and Weblogo	###############

	os.system('rm *fai')

	os.system('rm temp_cdhit.fasta')

	os.system('rm temp_info.txt')

	os.system('cat *clade.fasta > isolated_clusters.fasta')

	os.system("echo 'Total no of " + subtype_name +" CRISPR leader sequences clustered by TreeCluster:' & grep -o '>' isolated_clusters.fasta | wc -l" )

	os.system('rm isolated_clusters.fasta')

	os.system('mkdir 7_treecluster')
	time.sleep(0.1)
	os.system('mv *clade.fasta 7_treecluster')

	os.system('mv ' + job_name + '_cluster_info_' + str(cluster_treshold) + '.txt' + ' 7_treecluster')

	os.system('mv *cdhit_' + str(cd_hit) + '.fasta 3_cdhit')

	os.system('mkdir 8_aligned_clades')
	time.sleep(0.1)
	os.system('mv *clade_aligned.fasta 8_aligned_clades/')

	os.system('mkdir 9_weblogo_clades')

	os.system('mv *.eps 9_weblogo_clades/')

##################	Move all generated folders into Jobname folder	###############	

	os.system('mkdir ' + job_name)
	time.sleep(0.1)
	os.system('mv 2_Leader_seq_' + subtype_name + ' 3_cdhit 4_alignment 5_max_align 6_fasttree 7_treecluster 8_aligned_clades 9_weblogo_clades ' + job_name)

if __name__ == '__main__':
	Parameters1 = Params1(access_num, lower_bound, higher_bound, arr_orientation, sub_type, new_high_bound, full_name)
	Parameters1.getParam(filename)
	Parameters1.retParam()
	runCommand(access_num, lower_bound, higher_bound, arr_orientation, total_rep, sub_type, new_high_bound)
	out_log(access_num, lower_bound, higher_bound, arr_orientation, total_rep, sub_type, new_high_bound, full_name)
	Parameters2 = Params2(acc, subtype, new_acc_rc, new_acc_fw, new_acc, condition, sequence, seq_line, repeat, comp_name)
	Parameters2.getParam(log)
	Parameters2.newParam(raw)
	Parameters2.retParam()
	write(acc, subtype, new_acc_rc, new_acc_fw, condition, sequence, seq_line, repeat, comp_name)
	process1()
	Parameters4 = Params4(isim, clade_num)
	Parameters4.getParam(file_name)
	Parameters4.retParam()
	color(isim, clade_num)
	tidy()
	Parameters3 = Params3(name, clade)
	Parameters3.getParam(clades)
	Parameters3.retParam()
	call_clusters(name, clade)
	process2()
	
	






