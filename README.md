**To update to current version, run `conda install -c ursky snapt=0.4`

# SnapT - Small NcRNA Annotation Pipeline for (meta)Transcriptomic data

 SnapT is a **S**mall **n**on-coding RNA **a**nnotation **p**ipeline for **T**ranscriptomic or metatranscriptomic data. SnapT leverages transcriptomic or metatranscriptimic data to find, annotate, and quantify intergenic and anti-sense sRNA transcripts. To do this, SnapT aligns reads from a stranded RNAseq experiment to the reference (meta)genome, and then assembles the reads into transcripts. The transcripts are then intersected with the genome annotation as well as open reading frames to select for only transcripts that fall on non-coding regions, and further filtered to produce a final set of predicted small ncRNAs: 
 1. Intergenic transcripts must be at least 30nt away from any gene or ORF on both strands 
 2. Antisense transcripts must be 30nt away from any gene on their strand, but overlap with a gene on the opposite strand by at least 10nt. 
 3. Small peptides (<100nt) are not counted as a genes if they are encoded in a transcript that is more than 3 times their length. 
 4. The non-coding transcripts cannot contain any ORF greater than 1/3 of their length.
 4. Predicted non-coding transcripts near contig edges are discarded due to mis-annotation potential. 
 5. Small ncRNAs must be between 50nt and 500nt in length
 6. The transcripts must not have signifficant homology with any protein in the NCBI_NR database (query cover>30%, Bitscore>50, evalue<0.0001, and identity>30%).
 7. The transcripts must not have homology with any model (except for small RNAs) in the Rfam non-coding RNA database.

## SNAPT PIPELINE WORKFLOW
 ![SnapT small ncRNA annotation pipeline](https://i.imgur.com/pLdpRls.png)


## INSTALLATION

#### Conda installation:
 To start, download [miniconda2](https://conda.io/miniconda.html) and install it:
 ``` bash
 wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh #FOR LIXUX
 bash Miniconda2-latest-Linux-x86_64.sh
 ```
 
 Then simply install SnapT from the `ursky` conda channel (supports Linux64 and OsX):
 ``` bash
 conda install -c ursky snapt
 ```
 
 Download and index the NCBI NR protain database (you will need to input the nr.dmnd index into snapt with the -D or --nr option). The downloaded file should be 43GB (108GB unzipped), and indexing should take under an hour.
 ```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz
mv nr nr.faa
diamond makedb --in nr.faa -d nr
```
 Finally, download and index the Rfam non-coding RNA database (you will need to input the Rfam.cm database into snapt with the -R or --rfam option). The downloaded file should be 35MG (245MB unzipped), and indexing should take a few seconds.
 ```
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.1/Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm
```

#### Manual installation:
 You may want to manually install SnapT if you want better control over your environment, if you are installing on non-conventional system, or you just really dislike conda. In any case, you will need to manually install the [relevant prerequisite programs](https://github.com/ursky/SnapT/blob/master/conda_pkg/meta.yaml). When you are ready, download or clone this ripository and add the `SnapT/bin/` directory to to the `$PATH` or copy the `SnapT/bin/` contents into a directory that is under `PATH`. Thats it! 
 
 
## USAGE
 Example run of Snapt:
 ```
 snapt -1 READS/ALL_1.fastq -2 READS/ALL_2.fastq -g metagenomic_assembly.fasta -a metagenomic_assembly.gff -l 3000 -o SNAPT_OUT -t 48 --nr ../DATABASES/NCBI_nr/nr.dmnd --rfam ../DATABASES/rfam/Rfam.cm
 ```
 
 Help message
```bash
Usage: SnapT [options] -1 reads_1.fastq -2 reads_2.fastq -g genome.fa -o output_dir

SnapT options:
	-1 STR		forward transcriptome (or metatranscriptome) reads
	-2 SRT		reverse transcriptome (or metatranscriptome) reads (optional)
	-g STR		genome (or metagenome) fasta file
	-a STR		genome (or metagenome) annotation gtf/gff file (optional, but recommended)
	-l INT		minimum contig length (default=1000) for ncRNA annotation
	-o STR          output directory
	-t INT		number of threads (default=1)
	-d STR		NCBI_nr protein database DIAMOND index (see installation instructions for details)
	-rfam 		run analysis of ncRNAs detected against rfam database (rfam.cm) to annotate known ncRNAs/sRNAs (see installation instructions for details)

Aligment options:
	-r STR		rna-strandness: R or F for single-end, RF or FR for paired-end (default=FR)
	-I INT		min insert size (default=0)
	-X INT		max insert size (default=500)
	-m INT		gap distance to close transcripts (default=50)

	--version | -v	show current SnapT version
```

### Citing SnapT
SnapT is currently under early development. Stay tuned for future publication. 

### Acknowledgements
Authors of pipeline: [Gherman Uritskiy](https://github.com/ursky) and [Diego Gelsinger](https://github.com/dgelsin)

Principal Investigator: [Jocelyne DiRuggiero](http://bio.jhu.edu/directory/jocelyne-diruggiero/)

Institution: Johns Hopkins, [Department of Cell, Molecular, Developmental Biology, and Biophysics](http://cmdb.jhu.edu/) 

All feedback is welcome! For errors and bugs, please open a new Issue thread on this github page, and we will try to get things patched as quickly as possible. Please include the version of SnapT you are using (run `snapt -v`). For general questions about the conda impementation of this software, contact Gherman Uritskiy at guritsk1@jhu.edu. For general questions or suggestions about the pipeline itself, contact Diego Gelsinger at dgelsin1@jhu.edu. 

