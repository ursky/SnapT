#!/usr/bin/env bash
VERSION="0.1"

##############################################################################################################################################################
# SnapT: Small ncRNA annotation pipeline for transcriptomic data
# This pipeline aligns transcriptomic (or metatranscriptomic) reads to the genome (or metagenome), and finds transcripts that cannot be explained with coding
# regions. These transcripts are further process and curated until the final intergenic and anti-sense small ncRNAs are reported.
#
##############################################################################################################################################################


help_message () {
	echo ""
	echo "Usage: SnapT [options] -1 reads_1.fastq -2 reads_2.fastq -g genome.fa -o output_dir"
	echo ""
	echo "SnapT options:"
	echo "	-1 STR		forward transcriptome (or metatranscriptome) reads"
	echo "	-2 SRT		reverse transcriptome (or metatranscriptome) reads (optional)"
	echo "	-g STR		genome (or metagenome) fasta file"
	echo "	-a STR		genome (or metagenome) annotation gtf file (optional)"
	echo "	-o STR          output directory"
	echo ""
	echo "Aligment options:"
	echo "	-t INT          number of threads (default=1)"
	echo "	-r STR		rna-strandness [single-end: R or F; paired-end: RF or FR]. Only for strand-specific RNAseq"
	echo "	-I INT		min insert size (def: $MIN_INSERT_VALUE)"
	echo "	-X INT		max insert size (def: $MAX_INSERT_VALUE)"
	echo "	-m INT		gap distance to close transcripts (def: $GAP_VALUE)"
	echo ""
	echo "	--version | -v	show current SnapT version"
	echo "";}



comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }

# loading supplementary scripts
SNAPT_PATH=$(which snapt.sh)
BIN_PATH=${SNAPT_PATH%/*}
SOFT=${BIN_PATH}/snapt-scripts



########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################


# option defaults
STRAND_MAP_VALUE=auto
MIN_INSERT_VALUE=0
MAX_INSERT_VALUE=500
GAP_VALUE=50
STRAND_COUNT_VALUE=reverse
COUNT_ID_ATTRIBUTE=gene_id
THREADS=1

OUT=none
READS1=none
READS2=none
GENOME=none
ANNOTATION=none


# load in params
OPTS=`getopt -o hvt:1:2:g:a:o:r:I:X:m --long help,version -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
                -o) OUT=$2; shift 2;;
		-1) READS1=$2; shift 2;;
		-2) READS2=$2; shift 2;;
		-g) GENOME=$2; shift 2;;
		-a) ANNOTATION=$2; shift 2;;
		-t) THREADS=$2; shift 2;;
		-r) STRAND_MAP_VALUE=$2; shift 2;;
		-I) MIN_INSERT_VALUE=$2; shift 2;;
		-X) MAX_INSERT_VALUE=$2; shift 2;;
		-m) GAP_VALUE=$2; shift 2;;
                -h | --help) help_message; exit 0; shift 1;;
		-v | --version) echo SnapT v=${VERSION}; exit 0; shift 1;;
		--skip-refinement) refine=false; shift 1;;
                --) help_message; exit 1; shift; break ;;
                *) break;;
        esac
done

########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# check if all parameters are entered
if [[ $OUT == none ]] || [[ $READS1 == none ]] || [[ $GENOME == none ]]; then 
	comm "Some non-optional parameters (-1 -g -o) were not entered"
	help_message; exit 1
fi


# Checks for correctly configured snapt-scripts folder
if [ ! -s $SOFT/print_comment.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure that the snapt-scripts folder is in the same directory as snapt.sh"
fi

########################################################################################################
########################             BEGIN sRNA DISCOVERY PIPELINE!             ########################
########################################################################################################
announcement "BEGIN PIPELINE!"
comm "setting up output folder and copying relevant information..."
if [ ! -d $OUT ]; then
        mkdir $OUT
	if [ ! -d $OUT ]; then error "cannot make $OUT"; fi
else
        warning "Warning: $OUT already exists. It is recommended that you clear this directory to prevent any conflicts"
	#rm -r ${OUT}/*
fi


if [[ $ANNOTATION == none ]]; then
	comm "Annotation not given. Making OFT annotation with PRODIGAL"
	

 
 


########################################################################################################
########################     sRNA DISCOVERY PIPELINE SUCCESSFULLY FINISHED!!!   ########################
########################################################################################################
announcement "sRNA DISCOVERY PIPELINE FINISHED SUCCESSFULLY!"

