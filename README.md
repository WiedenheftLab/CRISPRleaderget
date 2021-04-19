# CRISPRleaderget

CRISPRleaderget Version 1 help:

CRISPRleaderget is a python program developed and tested in Mint Linux operating system. CRISPRleaderget.py should run under any unix based operating system that has a working 'python v2.7' (Not compatible with python3). 

This script is developed to isolate CRISPR leaders from CRISPRDetect output file in order to generate CRISPR leader phylogeny and identification of motifs within the leader sequences. If all the 3rd party dependencies are installed and available in user/system $PATH this automated script should run without any issues.

Please refer to CRISPRDetect for detection of CRISPR arrays, which the output file is required to use as an input for this script. Link to CRISPRDetect: https://github.com/ambarishbiswas/CRISPRDetect_2.2

	Syntax:

	python2 CRISPRleaderget.py -j demo -i demo.fasta -l demo_CRISPRDetect -s I-E

		OR

	python CRISPRleaderget.py -j demo -i demo.fasta -l demo_CRISPRDetect -s All -f 200 -cdc 0.95 -r Yes -cm avg_clade -ct 1.5 -t 50 

Installation

CRISPRleaderget does not require any installation. However, it uses 3rd party tools to perform automated pipeline. Please make sure that the following 3rd party tools are installed in your system.

CRISPRleaderget dependencies:

	Cd-hit		https://github.com/weizhongli/cdhit
	MAFFT		https://mafft.cbrc.jp/alignment/software/source.html
	SamTools	http://www.htslib.org/download/
	FastTree	http://www.metagenomics.wiki/tools/phylogenetic-tree/construction/fasttree
	TreeCluster.py	https://github.com/niemasd/TreeCluster
	WebLogo		https://pypi.org/project/weblogo/
	Color_tree	https://github.com/mooreryan/color_tree
	Maxalign	http://www.cbs.dtu.dk/services/MaxAlign/

Once all the dependencies are installed, please check that they are successfully installed and available in the user/system PATH by typing the following:

	cdhit -h
	mafft -h
	samtools
	fasttree -h
	TreeCluster.py -h
	weblogo -h
	color_tree -h
	
For MaxAlign, direct to the folder where the sowftware is installed check if it is installed successfull with the following command (No need to include it to user/system PATH):
	
	perl maxalign.pl

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
