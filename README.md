**Current SnapT version = 0.1**. To update, run `conda install -c ursky snapt=0.1`

# SnapT - **S**mall **n**cRNA **a**nnotation **p**ipeline for **t**ranscriptomic data

 SnapT is a small non-coding RNA discovery pipeline. SnapT leverages transcriptomic or metatranscriptimic data to find, annotate, and quantify intergenic and anti-sense sRNA transcripts. 

## SnapT pipeline workflow
 To make in draw.io


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
 
#### Manual installation:
 You may want to manually install SnapT if you want better control over your environment, if you are installing on non-conventional system, or you just really dislike conda. In any case, you will need to manually install the [relevant prerequisite programs](https://github.com/ursky/SnapT/blob/master/conda_pkg/meta.yaml). When you are ready, download or clone this ripository and add the `SnapT/bin/` directory to to the `$PATH` or copy the `SnapT/bin/` contents into a directory that is under `PATH`. Thats it! 
 
 
## USAGE
 Usage instructions and example run with SnapT. To be added.
```bash
Usage: SnapT [options] -1 reads_1.fastq -2 reads_2.fastq -g genome.fa -o output_dir

SnapT options:
	-1 STR		forward transcriptome (or metatranscriptome) reads
	-2 SRT		reverse transcriptome (or metatranscriptome) reads (optional)
	-g STR		genome (or metagenome) fasta file
	-a STR		genome (or metagenome) annotation gtf file (optional)
	-o STR          output directory

Aligment options:
	-t INT          number of threads (default=1)
	-r STR		rna-strandness: R or F for single-end, RF or FR for paired-end (default=FR). Only for strand-specific RNAseq
	-I INT		min insert size (def: 0)
	-X INT		max insert size (def: 500)
	-m INT		gap distance to close transcripts (def: 50)

	--version | -v	show current SnapT version
```

### Citing SnapT
SnapT is currently under early development. Stay tuned for future publication. 

### Acknowledgements
Authors of pipeline: [Gherman Uritskiy](https://github.com/ursky) and [Diego Gelsinger](https://github.com/dgelsin)

Principal Investigators: [Jocelyne DiRuggiero](http://bio.jhu.edu/directory/jocelyne-diruggiero/) and [James Taylor](http://bio.jhu.edu/directory/james-taylor/)

Institution: Johns Hopkins, [Department of Cell, Molecular, Developmental Biology, and Biophysics](http://cmdb.jhu.edu/) 

All feedback is welcome! For errors and bugs, please open a new Issue thread on this github page, and we will try to get things patched as quickly as possible. Please include the version of SnapT you are using (run `snapt -v`). For general questions about the conda impementation of this software, contact Gherman Uritskiy at guritsk1@jhu.edu. For general questions or suggestions about the pipeline itself, contact Diego Gelsinger at dgelsin1@jhu.edu. 

