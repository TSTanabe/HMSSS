# HMS-S-S
HMS-S-S: a tool for the identification of sulfur metabolism-related genes and analysis of operon structures in genome and metagenome assemblies

## What does HMS-S-S do?

## Installation
HMS-S-S is written in perl and is dependent on HMMER3, MySQL and Bioperl. Installation instructions are given below.
You can download HMS-S-S directly from GitHub.

1. Download the latest release from github: https://github.com/TSTanabe/HMSSS
2. Install all dependencies if not already installed (see below)
3. Unpack the files in a directory and open a terminal in this directory
4. Start the script with option -h for further instructions on handling or read this manual: ```$ perl HMSSS.pl -h```
5. Alternative to running the command line version a graphical user interface version can be started: ```$ perl HMSSS_gui.pl```

Dependencies may already be installed by default, but all packages are listed for completeness.
### Perl

The Perl interpreter is usually included in the default installation of Ubuntu, but can otherwise be downloaded and installed via
```$ sudo apt-get install perl```

For the **command line version** and the **gui version** following perl packages are required. Packages that are required but possibly not installed by the default installation can be downloaded from CPAN and are marked:

*Cwd
***DBI;**
***DBD::mysql**
*Error
*File::Basename
*File::Copy
*File::Find
*Getopt::Long
*Getopt::Std
*IO::Uncompress::Gunzip
*IO::Compress::Gzip
***Parallel::ForkManager**

For the **gui version** these additional dependencies are required:

*Tk
*Tk::Optionmenu
*Tk::ROText
*Tk::BrowseEntry
*Tk::Pane
*Tk::PNG
All TK libraries can be installed via 
```sudo apt-get install -y perl-tk```



### Bioperl
*Bio::Seq;
*Bio::SeqIO;
*Bio::SearchIO;
*Bio::SearchIO::FastHitEventBuilder;
*Bio::Tools::Run::StandAloneBlastPlus;
*Bio::DB::Taxonomy;
### Mysql/MariaDB

### HMMER3

### Prodigal

### NCBI Taxonomy Database

### NCBI Blast+

## Running HMS-S-S

## Running HMS-S-S with graphical user interface
