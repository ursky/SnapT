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
	echo "	-r STR		rna-strandness: R or F for single-end, RF or FR for paired-end (default=FR). Only for strand-specific RNAseq"
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
STRAND_MAP_VALUE=FR
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


########################################################################################################
########################         ALIGN RNA READS TO GENOME WITH HISAT2          ########################
########################################################################################################
announcement "ALIGN RNA READS TO GENOME WITH HISAT2"

comm "Building Hisat2 index from reference genome"
mkdir ${OUT}/hisat2
cp $GENOME ${OUT}/hisat2/genome.fa
hisat2-build ${OUT}/hisat2/genome.fa ${OUT}/hisat2/hisat2_index
if [[ $? -ne 0 ]]; then error "Hisat2 index could not be build. Exiting..."; fi


comm "Aligning $READS1 and $READS2 to $GENOME with hisat2"
if [[ $READS2 == none ]];then
	hisat2 -p 20 --verbose --no-spliced-alignment\
	 --rna-strandness "$STRAND_MAP_VALUE"\
	 -I "$MIN_INSERT_VALUE"\
	 -X "$MAX_INSERT_VALUE"\
	 -x ${OUT}/hisat2/hisat2_index\
	 -U $READS1 \
	 -S ${OUT}/hisat2/alignment.sam
else
	hisat2 -p 20 --verbose --no-spliced-alignment\
	 --rna-strandness "$STRAND_MAP_VALUE"\
	 -I "$MIN_INSERT_VALUE"\
	 -X "$MAX_INSERT_VALUE"\
	 -x ${OUT}/hisat2/hisat2_index\
	 -1 $READS1 -2 $READS2 \
	 -S ${OUT}/hisat2/alignment.sam
fi
if [[ $? -ne 0 ]]; then error "Hisat2 alignment failed. Exiting..."; fi


comm "Converting aligment.sam to .bam format"
samtools view -bS ${OUT}/hisat2/alignment.sam > ${OUT}/hisat2/alignment.bam
if [[ $? -ne 0 ]]; then error "Samtools sam to bam conversion failed. Exiting..."; fi


comm "Sorting hisat2 alignment.sam by coordinate"
samtools sort ${OUT}/hisat2/alignment.bam -o ${OUT}/hisat2/alignment.sorted.bam
if [[ $? -ne 0 ]]; then error "Samtools sorting failed. Exiting..."; fi


comm "Building IGV index from hisat2 alignment_sorted.bam"
samtools index ${OUT}/hisat2/alignment.sorted.bam
if [[ $? -ne 0 ]]; then error "Samtools indexing failed. Exiting..."; fi



########################################################################################################
########################              ASSEMBLE DE-NOVO TRANSCRIPTS              ########################
########################################################################################################
announcement "ASSEMBLE DE-NOVO TRANSCRIPTS"

comm "Building de-novo transcripts and transcriptome expression file"
mkdir ${OUT}/transcript_assembly
stringtie ${OUT}/hisat2/alignment.sorted.bam \
	-G $ANNOTATION -m $GAP_VALUE \
	-o ${OUT}/transcript_assembly/stringtie_transcripts.gtf \
	-A ${OUT}/transcript_assembly/gene_abundance.tab \
	-C ${OUT}/transcript_assembly/expressed_transcripts.gtf -e -B
if [[ $? -ne 0 ]]; then error "Stringtie transcript assembly failed. Exiting..."; fi


mkdir ${OUT}/gene_analysis
mv ${OUT}/transcript_assembly/stringtie_transcripts.gtf ${OUT}/gene_analysis
mv ${OUT}/transcript_assembly/*ctab ${OUT}/gene_analysis
mv ${OUT}/transcript_assembly/*tab ${OUT}/gene_analysis


comm "Building Reference-based transcriptome expression file"
stringtie ${OUT}/hisat2/alignment.sorted.bam \
	-o ${OUT}/transcript_assembly/reference_transcripts.gtf \
	-G $ANNOTATION -m $GAP_VALUE
if [[ $? -ne 0 ]]; then error "Stringtie reference transcript quantitation failed. Exiting..."; fi


########################################################################################################
########################                        PROGRESS END                    ########################
########################################################################################################



#Note: look into stringmerge for unifying transcripts
echo ========Comparing reference-based transcriptome expression file against reference annotation file========
if [ "$#" -eq 4 ]; then #counts number of input files given and adjusts variables offset -1 (single-end mode)
       ~/scratch/Diego/gffcompare-0.9.9b.Linux_x86_64/gffcompare -r "$2" "$4".ref-based_transcriptome_assembly_expression.gtf -o "$4".ref-based_gffcompare_transcriptome_assembly
elif [ "$#" -eq 5 ]; then #counts number of input files given and adjusts variables offset -1 (single-end mode)
       ~/scratch/Diego/gffcompare-0.9.9b.Linux_x86_64/gffcompare -r "$3" "$5".ref-based_transcriptome_assembly_expression.gtf -o "$5".ref-based_gffcompare_transcriptome_assembly
fi

mkdir "ref-based_transcriptome"
mv *ref-based_transcriptome_assembly* ref-based_transcriptome
mv *gffcompare* ref-based_transcriptome
echo ========Reference-based Transcriptome comparison completed========

module load python/2.7.10

###Pull out novel transcripts based on class code###
echo ========Finding sRNAs in Reference-based transcriptome========
if [ "$#" -eq 4 ]; then #counts number of input files given and adjusts variables offset -1 (single-end mode)
        python ../antisense_sRNA_gff_maker.py ref-based_transcriptome/"$4".ref-based_gffcompare_transcriptome_assembly*.gtf > "$4".ref-based_antisense_sRNA.gtf | python ../intergenic_sRNA_gff_maker.py ref-based_transcriptome/"$4".ref-based_gffcompare_transcriptome_assembly*.gtf > "$4".ref-based_intergenic_sRNA.gtf
elif [ "$#" -eq 5 ]; then #counts number of input files given and adjusts variables offset +1 (paired-end mode)
        python ../antisense_sRNA_gff_maker.py ref-based_transcriptome/"$5".ref-based_gffcompare_transcriptome_assembly*.gtf > "$5".ref-based_antisense_sRNA.gtf | python ../intergenic_sRNA_gff_maker.py ref-based_transcriptome/"$5".ref-based_gffcompare_transcriptome_assembly*.gtf > "$5".ref-based_intergenic_sRNA.gtf
fi
mkdir "ref-based_sRNAs"
mv *antisense_sRNA* ref-based_sRNAs
mv *intergenic_sRNA* ref-based_sRNAs
echo ========reference-based sRNAs found========

###de novo analysis
echo ========Building de novo transcriptome expression file========
if [ "$#" -eq 4 ]; then #counts number of input files given and adjusts variables offset -1 (single-end mode)
       ~/scratch/Diego/stringtie-1.3.2b.Linux_x86_64/stringtie alignments/"$4".alignment.cordsorted.bam -o "$4".denovo_transcriptome_assembly_expression.gtf -m "$GAP_VALUE"
elif [ "$#" -eq 5 ]; then #counts number of input files given and adjusts variables offset +1 (paired-end mode)
       ~/scratch/Diego/stringtie-1.3.2b.Linux_x86_64/stringtie alignments/"$5".alignment.cordsorted.bam -o "$5".denovo_transcriptome_assembly_expression.gtf -m "$GAP_VALUE"
fi
echo ========Building de novo transcriptome.gtf file finished========

echo ========Comparing de novo transcriptome expression file against reference annotation file========
if [ "$#" -eq 4 ]; then #counts number of input files given and adjusts variables offset -1 (single-end mode)
       ~/scratch/Diego/gffcompare-0.9.9b.Linux_x86_64/gffcompare -r "$2" "$4".denovo_transcriptome_assembly_expression.gtf -o "$4".denovo_gffcompare_transcriptome_assembly
elif [ "$#" -eq 5 ]; then #counts number of input files given and adjusts variables offset -1 (paired-end mode)
       ~/scratch/Diego/gffcompare-0.9.9b.Linux_x86_64/gffcompare -r "$3" "$5".denovo_transcriptome_assembly_expression.gtf -o "$5".denovo_gffcompare_transcriptome_assembly
fi
mkdir "denovo_transcriptome"
mv *denovo* denovo_transcriptome
echo ========de novo Transcriptome comparison completed========

###Pull out novel transcripts based on class code###
echo ========Finding sRNAs in de novo transcriptome========

if [ "$#" -eq 4 ]; then #counts number of input files given and adjusts variables offset -1 (single-end mode)
        python ../antisense_sRNA_gff_maker.py denovo_transcriptome/"$4".denovo_gffcompare_transcriptome_assembly*.gtf > "$4".denovo_antisense_sRNA.gtf | python ../intergenic_sRNA_gff_maker.py denovo_transcriptome/"$4".denovo_gffcompare_transcriptome_assembly*.gtf > "$4".denovo_intergenic_sRNA.gtf
elif [ "$#" -eq 5 ]; then #counts number of input files given and adjusts variables offset +1 (paired-end mode)
        python ../antisense_sRNA_gff_maker.py denovo_transcriptome/"$5".denovo_gffcompare_transcriptome_assembly*.gtf > "$5".denovo_antisense_sRNA.gtf | python ../intergenic_sRNA_gff_maker.py denovo_transcriptome/"$5".denovo_gffcompare_transcriptome_assembly*.gtf > "$5".denovo_intergenic_sRNA.gtf
fi
mkdir "denovo_sRNAs"
mv *antisense_sRNA* denovo_sRNAs
mv *intergenic_sRNA* denovo_sRNAs	

 
 


########################################################################################################
########################     sRNA DISCOVERY PIPELINE SUCCESSFULLY FINISHED!!!   ########################
########################################################################################################
announcement "sRNA DISCOVERY PIPELINE FINISHED SUCCESSFULLY!"

