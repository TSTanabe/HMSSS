# HMS-S-S
HMS-S-S: a tool for the identification of sulfur metabolism-related genes and analysis of operon structures in genome and metagenome assemblies

## What does HMS-S-S do?

## Installation
HMS-S-S is written in perl and is dependent on HMMER3, MySQL and Bioperl. Installation instructions are given below.
You can download HMS-S-S directly from GitHub.

1. Download the latest release from github: https://github.com/TSTanabe/HMSSS
2. Install all dependencies if not already installed (see below)
3. Unpack the files in a directory and open a terminal in this directory
4. Start the script with option -h for further instructions on handling or read this manual: ```perl HMSSS.pl -h```
5. Alternative to running the command line version a graphical user interface version can be started: ```perl HMSSS_gui.pl```

##Prerequisites
Dependencies may already be installed by default, but all packages are listed for completeness.
### Perl

The Perl interpreter is usually included in the default installation of Ubuntu, but can otherwise be downloaded and installed via

```sudo apt-get install perl```

For the **command line version** and the **gui version** following perl packages are required. Packages that are required but possibly not installed by the default installation can be downloaded from CPAN and are marked:

* Cwd
* **DBI;**
* **DBD::mysql**
* Error
* File::Basename
* File::Copy
* File::Find
* Getopt::Long
* Getopt::Std
* IO::Uncompress::Gunzip
* IO::Compress::Gzip
* **Parallel::ForkManager**

For the **gui version** these additional dependencies are required:

* Tk
* Tk::Optionmenu
* Tk::ROText
* Tk::BrowseEntry
* Tk::Pane
* Tk::PNG

All TK libraries can be installed via 
```sudo apt-get install -y perl-tk```



### Bioperl
Bioperl is a toolbox of perl packages written specifically for computational molecular biology applications. These include but are not limited to reading sequences from fasta formats and processing BLAST or HMMER reports. To perform a hmmsearch and analysis HMSSS requires some dependencies in the path. Brief instructions on installation are given, but can also be obtained from the installation notes of the packages itself. Further information can be found at https://bioperl.org. The installation can be done with the following command:

```sudo apt-get install bioperl```

Alternatively, the packages can also be installed via CPAN. The following packages from bioperl are required:
* Bio::Seq;
* Bio::SeqIO;
* Bio::SearchIO;
* Bio::SearchIO::FastHitEventBuilder;
* Bio::Tools::Run::StandAloneBlastPlus;
* Bio::DB::Taxonomy;


### Mysql/MariaDB

MySQL and MariaDB are completely interchangeable relational database management systems that even share the same syntax and commands. However, the source code of MariaDB is freely available. HMS-S-S uses a database based on MariaDB, but could just as well use a MySQL database. For installation and configuration under Ubuntu are some steps necessary, starting with:
```
  $ sudo apt-get install mariadb-server
  $ mysql --version        #Check for correct installation
  mysql  Ver 15.1 Distrib 10.3.32-MariaDB, for debian-linux-gnu
  (i686) using readline 5.2
  $ sudo mysql -u root -p
  Enter password:	#Enter your ubuntu password
  
  #Choose a Username a password for your database account
  > GRANT ALL ON *.* TO Username IDENTIFIED BY "password";
  > GRANT FILE ON *.* TO Username;
  > QUIT
```
To log in with root, the password of your computer is required. After that, you can create a local user with your own username and an independent password. With GRANT ALL the user is created. with GRANT FILE the write and read rights are assigned.

Now it is possible to log in with the local newly created user with the following command. Afterwards a database must be created, whereby the name can be freely selected again. This database will later be used to store the search results. With the last command the database schema is configured. this step must be executed from the same directory where schema.sql is stored, or specify the complete path to this file.

```
  $ mysql -u Username -p
  Enter password:
  
  > CREATE DATABASE database_name;
  > QUIT
  
  #Enter your custom username, password and database_name
  $ mysql -u Username -p database_name < schema.sql

```




### HMMER3

HMMER is a software package that uses profile hidden markov models (HMM) to detect homologous sequences. Compatible HMMs can be obtained from databases such as Pfam or Interpro, or created from a custom sequence alignment via the HMMER package itself.
The website http://hmmer.org/ provides further reading on the HMMER algorithm and functions. Installation can either be done from the website or by:

```sudo apt-get install hmmer```

### Prodigal
Prodigal (Prokaryotic Dynamic Programming Genefinding Algorithm) is a microbial (bacterial and archaeal) gene finding program, that learns the ORF prediction directly by the provided data. It can be installed via:

```sudo apt-get install prodigal```

### NCBI Taxonomy Database

For assemblies derived from NCBI the NCBI Taxonomy database can be used to add taxonomic information. The database taxdmp.zip can be downloaded from NCBI FTP server: https://ftp.ncbi.nih.gov/pub/taxonomy/
The unpacked taxdump folder containing the .dmp files has to be moved to HMS-S-S source folder.

### NCBI Blast+
NCBI BLAST+ enables local usage of the BLAST algorithm. It is not necessarily required by HMS-S-S but can be used as an addition to HMMER.

```sudo apt-get install ncbi-blast+```

## Running HMS-S-S
HMSSS takes many options directly from the command line without the graphical interface. Except for the login data of the local database which are essential, all other options are not mandatory and can be executed with default values. It is recommended to give a path to the directory to be analyzed. Analysis can be started with the command:

```perl HMSSS.pl -user Username -password Password -database Database_name -fasta_directory ~\directory ```

**Input:** Directory with nucleotide fasta files or protein fasta and corresponding gff3 formatted files.
           Suffix should be .fa, .fna, .fa.gz, .faa, .faa.gz, .gff or .gff.gz
**Output:** 

Possible options and the default values can be printed with the -h option.

```perl HMSSS.pl -h```

Which will output the following text.
 
```
Options:
   -h or -help		Print the help dialog
   
Database Connection
   -user		Mysql database user name
   -password		Mysql database password
   -database		Name of the Mysql database
   
Protein sequence search
   -algorithm		Name of the search algorithm: hmmer or blastp
   			(Default: hmmer)
   -cutoff_type 	Options: standard_cutoff or evalue_cutoff or
   		 	score_cutoff (Default: standard_cutoff)
   			standard_cutoff requires a threshold file as input
   			evalue_cutoff takes a constant evalue threshold from -cutoff
   			score_cutoff takes a constant score threshold from -cutoff
   			
   -cutoff		Score or E-value for search cutoff (Default: 10)
   -threshold_file	Path to tab separated file with threshold values for each HMM
   			or query protein.
   			Format: protein_type \t threshold (Default: /Thresholds)
   			
   -query_file		HMM library or fasta file for blastp (Default: HMMlib)
   -fasta_directory	Directory with fasta files to be analysed (Default: /Genomes)
   -gff3_directory	Directory with gff3 files if protein fasta is provided 
   			(Default: /Genomes)
   -cores		Number of cores used by hmmsearch (Default: 2)
   -cleanreports	Remove hmmsearch report files after finishing search
   -redosearch		Search all files in the directory even if it is already 
   			in the database
   -redoreports         Skip hmmsearch only use results from existing report files
   -prodigal		Translation nucleotide fasta with prodigal
   
Gene cluster prediction and keyword assignment
   -distance		Distance in nucleotides between two genes to considered
   			them as a gene cluster	(Default: 3500)
   -pattern_file	Keywords and patterns of synthenic genes
   			(Default: /Patterns.txt)
   -output_directory	Output directory for results from the searched fastas 
   			(Default: /Output)
            
Process control options
   -start_clustering	Start at the find synteny step
   -start_keywords	Start at the keyword assignment step
   -start_output	Start at the Output. Only outputs results from files in
   			the fasta_directory path
   -output_DB_all	Output of all entrys in the database
```

**HMSSS.ini** To simplify the input of frequently used options, they can be stored in the HMSSS.ini file. Parameters specified in this file will overwrite the default values of the program. If additional options are used in the command line that differ from the values specified in HMSSS.ini, the options from the command line are used by the program. The parameters should be written between the quotation marks in the .ini file.

**Workflow** The processing steps of HMS-S-S are as followed:


1. If prodigal is installed and HMSSS was started with -prodigal option, the ORF prediction and translation for nucleotide fasta files will be done and genomic features will be saved in a .gff3 file.
2. All protein fasta files in the directory will be searched by hmmsearch consecutively with the specified HMM library
3. Results will be stored in the local database, with the genomic features and the corresponding protein sequence 

4. Gene clusters are predicted. Two genes are considered to be syntenic if they are within a maximum distance of nucleotides from each other.
5. Keyword assignment. Given assemblies are searched each pattern of gene occurence defined in the Patterns.txt file. A genecluster matching one or several patterns is given one keyword for each match with a pattern.


## Running HMS-S-S with graphical user interface

HMS-S-S has a graphical user interface that supports the program flow. The interface requires no installation but depends on the TK package and can be started directly from the directory with the command:

```perl HMSSS_gui.pl```

All input FASTA should be located in a single folder. Nucleotide FASTA formats must have a .fa, .fna, .fa.gz or .fna.gz file suffix. For protein FASTA format a .faa.gz suffix is required, as well as a corresponding gff3 file with the same filename but a gff3.gz file suffix. HMS-S-S will iteratively analyse all files in the directory.

###
