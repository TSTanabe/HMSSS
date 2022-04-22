#!/usr/bin/perl

use strict;
use warnings;
use utf8;

use Getopt::Long;
use File::Basename;
use File::Copy;
use File::Find;
use Cwd qw(getcwd);

use DBI;
use DBD::mysql;


my $cwd = getcwd;
my $date = localtime();

my $db_user = "";
my $db_name = "";	
my $db_pass = "";

my $algorithm		="hmmer";
my $cutoff_type	="standard_cutoff";
my $cutoff		="10";
my $threshold_file	="$cwd/Thresholds";
my $query_file		="$cwd/HMMlib";
my $fasta_directory	="$cwd/Genomes";
my $gff3_directory	="$cwd/Genomes";
my $cores		="2";
my $cleanreports	="0";
my $redosearch		="0";
my $redoreports	="0";
my $prodigal		="0";

my $distance		="3500";
my $pattern_file	="$cwd/Patterns.txt";

my $output_directory 	="$cwd/Output/";
my $DB_all 		="0";

my $start_search 	="0";
my $start_clustering	="0";
my $start_keywords	="0";
my $start_output	="0";

my $help		="0";
#Initilize all parameters which usually would be defined through the GUI via the ini file
my $ini_hashref;

if(-f "$cwd/HMSSS.ini"){
$ini_hashref = &getIniList("$cwd/HMSSS.ini"); #Hashtable for standard parameters
}



###########Check if search parameter are initilized###############
if(defined $ini_hashref->{user}){$db_user = $ini_hashref->{user};}
if(defined $ini_hashref->{password}){$db_pass = $ini_hashref->{password};}
if(defined $ini_hashref->{database}){$db_name = $ini_hashref->{database};}

if(defined $ini_hashref->{algorithm}){
	if( $ini_hashref->{algorithm} eq "hmmer" || $ini_hashref->{algorithm} eq "blastp"){
	$algorithm = $ini_hashref->{algorithm};}
	else{&warning("Unrecognized search algorithm in ini.txt");}
}

if(defined $ini_hashref->{cutoff_type}){
	if( $ini_hashref->{cutoff_type} eq "standard_cutoff" || $ini_hashref->{cutoff_type} eq "evalue_cutoff" || $ini_hashref->{cutoff_type} eq "score_cutoff"){
	$cutoff_type = $ini_hashref->{cutoff_type}}
	else{&warning("Unrecognized cutoff type in ini.txt");}
}

if(defined $ini_hashref->{cutoff}){
	if($ini_hashref->{cutoff} =~ /^\d+$/){$cutoff = $ini_hashref->{cutoff}}
	else{&warning("Cutoff is not a number in ini.txt");}
}

if(defined $ini_hashref->{threshold_file}){
	if(-f $ini_hashref->{threshold_file}){$threshold_file = $ini_hashref->{threshold_file};}
	else{&warning("Threshold file from ini.txt was not found ".$ini_hashref->{threshold_file});}
}

if(defined $ini_hashref->{query_file}){
	if(-f $ini_hashref->{query_file}){$query_file = $ini_hashref->{query_file};}
	else{&warning("Query file from ini.txt was not found ".$ini_hashref->{query_file});}
}

if(defined $ini_hashref->{fasta_directory}){
	if(-d $ini_hashref->{fasta_directory}){$fasta_directory = $ini_hashref->{fasta_directory};}
	else{&warning("Fasta file directory from ini.txt was not found ".$ini_hashref->{fasta_directory});}
}
if(defined $ini_hashref->{gff3_directory}){
	if(-d $ini_hashref->{gff3_directory}){$gff3_directory = $ini_hashref->{gff3_directory};}
	else{&warning("GFF3 file directory from ini.txt was not found ".$ini_hashref->{gff3_directory});}
}

if(defined $ini_hashref->{cores}){
	if( $ini_hashref->{cores} =~ /^\d+$/){$cores = $ini_hashref->{cores}}
	else{&warning("Cores is not a number in ini.txt");}
}

if(defined $ini_hashref->{cleanreports}){
	if( $ini_hashref->{cleanreports} == 0 || $ini_hashref->{cleanreports} == 1){$cleanreports = $ini_hashref->{cleanreports}}
	else{&warning("Cleanreports option was not boolean in ini.txt");}
}

if(defined $ini_hashref->{redosearch}){
	if( $ini_hashref->{redosearch} == 0 || $ini_hashref->{redosearch} == 1){$redosearch = $ini_hashref->{redosearch}}
	else{&warning("Redosearch option was not boolean in ini.txt");}
}


if(defined $ini_hashref->{redoreports}){
	if($ini_hashref->{redoreports} == 0 || $ini_hashref->{redoreports} == 1){$redoreports = $ini_hashref->{redoreports}}
	else{&warning("Redoreports option was not boolean in ini.txt");}
}

if(defined $ini_hashref->{prodigal}){
	if($ini_hashref->{prodigal} == 0 || $ini_hashref->{prodigal} == 1){$prodigal = $ini_hashref->{prodigal}}
	else{&warning("Prodigal was not boolean in ini.txt");}
}

if(defined $ini_hashref->{distance}){
	if( $ini_hashref->{distance} =~ /^\d+$/){$cores = $ini_hashref->{distance}}
	else{&warning("Distance is not a number in ini.txt");}
}

if(defined $ini_hashref->{pattern_file}){
	if(-f $ini_hashref->{pattern_file}){$pattern_file = $ini_hashref->{pattern_file};}
	else{&warning("Pattern file from ini.txt was not found ".$ini_hashref->{pattern_file});}
}

if(defined $ini_hashref->{output_directory}){
	if(-d $ini_hashref->{output_directory}){$output_directory = $ini_hashref->{output_directory};}
	else{&warning("Output directory from ini.txt was not found. Creating new directory ".$ini_hashref->{output_directory});
	mkdir($ini_hashref->{output_directory});
	}
}
GetOptions ("user=s" => \$db_user,
            "password=s" => \$db_pass,
            "database=s" => \$db_name,            
            "algorithm=s" => \$algorithm,
            "cutoff_type=s" => \$cutoff_type,
            "cutoff=i"	=> \$cutoff,
            "threshold_file=s" => \$threshold_file,
            "query_file=s" => \$query_file,
            "fasta_directory=s"=> \$fasta_directory,
            "gff3_directory=s" => \$gff3_directory,
            "cores=i"	=> \$cores,
            "cleanreports" => \$cleanreports,
            "redosearch" => \$redosearch,
            "redoreports" => \$redoreports,
            "prodigal" => \$prodigal,
            "distance=i" => \$distance,
            "pattern_file=s" => \$pattern_file,
            "output_directory" => \$output_directory,
            
            #TODO Control of STAGE for starting point
            "start_search" => \$start_search,
            "start_clustering" =>\$start_clustering,
            "start_keywords" =>\$start_keywords,
            "start_output" => \$start_output,	
            #TODO Output option for whole output from database
            "output_DB_all" => \$DB_all,            
            "help" => \$help,
            "h" => \$help,
            )or die("Error in command line arguments\n");
if($help){&help();exit 0;}
#TODO Test auf bugs und stabilität bei er eingabe von parametern

if($help){&help();exit 0;}

###########Check if database connection is initilized#############
#Nach get Options DB connection prüfen, weil gar kein default wert im hardcode steht.
my $dbh;

if($db_user && $db_name && $db_pass){
$dbh = DBI->connect("DBI:mysql:database=$db_name", $db_user, $db_pass) or &warning("Database connection failed");
}else{print "Database login data is missing. \n Options: -user [username] -password [password] -database [name] \n have to be initilized \n";exit;}
#TODO plan for failed database connection to establish

unless(-d $fasta_directory  ){
warn "ERROR:\t No such directory found. Skip protein search step\n";
$start_clustering = 1;
$DB_all = 1;}
unless(-d $output_directory ){mkdir($output_directory);}



my $COMMAND;

############################# STAGE Search  ##################################################

unless($start_clustering || $start_keywords || $start_output){
$COMMAND = "perl bin/Search.pl -u $db_user -p $db_pass -b $db_name -a $algorithm -c $cutoff_type -v $cutoff -t $threshold_file -q $query_file -f $fasta_directory -g $gff3_directory -n $cores -r $redosearch -l $cleanreports -o $prodigal -k $redoreports ";
system($COMMAND);
}
############################# STAGE Cluster  ##################################################
unless($start_keywords || $start_output){
$COMMAND = "perl bin/Cluster.pl -d $distance -f $pattern_file -u $db_user -b $db_name -p $db_pass -m synteny";

system($COMMAND);
}
############################# STAGE keywords ##################################################
unless($start_output){
$COMMAND = "perl bin/Cluster.pl -d $distance -f $pattern_file -u $db_user -b $db_name -p $db_pass -m keyword";

system($COMMAND);
}
############################# STAGE Output   ##################################################
print "\n";
print "#"x 40;
print "\n# Starting output\n";
print "#"x 40;
print "\n";

&output($DB_all);
###########################################################################################
################ SUBROUTINES ##############################################################
################ FOR OUTPUT  ##############################################################
###########################################################################################
###
###
###
###
###

sub output{
my $DB_all = shift;
my ($size,$faaFiles);

if(-d $fasta_directory){
	$faaFiles = &getallDirs($fasta_directory,".faa\.gz");
	$size = scalar @{ $faaFiles };
}else{
	$size = 0;
	$DB_all = 1;
}


my $count = 1;
my @result_table;
my %hash_keywords;
my @headings = qw(genome_ID protein_ID cluster_ID contig start end strand HMM score cluster_type);

#Iterate all cluster_IDs and concat the keywords
my $query = "SELECT distinct cluster_ID,keyword from Functions;";
my $mysql = $dbh->prepare($query);
$mysql->execute();

while(@result_table = $mysql->fetchrow_array()){ #TODO check if this is correctly working
	if(defined $hash_keywords{$result_table[0]}){
	my $tmp = $hash_keywords{$result_table[0]};
	$tmp .= "-".$result_table[1];
	$hash_keywords{$result_table[0]} = $tmp;
	}else{
	$hash_keywords{$result_table[0]} = $result_table[1];
	}
}	
	
$date =~ s/ /_/g;	
print "Writing results to ".$output_directory."$date\_Results\.tsv\n";
open my $fh,">",$output_directory."$date\_Results\.tsv" or warn "Could not open output file $!";
print $fh join("\t",@headings); #headings for columns
print $fh "\n";
#retrieve results from database
#if db all is boolean == 1 dann alle ergebnisse unabhängig von der genome_ID;
#if fasta directory empty dann ebenfalls alle ergebnisse;
if($DB_all || $size == 0){
	$query = "SELECT Genomes.genome_ID, Proteins.protein_ID,Proteins.cluster_ID,Proteins.contig,Proteins.start,Proteins.end,Proteins.strand,Proteins.HMM,Proteins.score  FROM Proteins LEFT JOIN Genomes ON Proteins.genome_ID = Genomes.genome_ID ORDER BY Proteins.genome_ID,contig,start;";
	$mysql = $dbh->prepare($query);
	
		
		$mysql->execute();
		while(@result_table = $mysql->fetchrow_array()){
		my $string;
		if (defined $result_table[2] && defined $hash_keywords{$result_table[2]}){
		$string = join("\t",@result_table)."\t".$hash_keywords{$result_table[2]}."\n";
		}else{
		$result_table[2] = "-";
		$string = join("\t",@result_table)."\t"."\n";
		}
		print $fh $string;
		}
	
}else{
	$query = "SELECT Genomes.genome_ID, Proteins.protein_ID,Proteins.cluster_ID,Proteins.contig,Proteins.start,Proteins.end,Proteins.strand,Proteins.HMM,Proteins.score  FROM Proteins LEFT JOIN Genomes ON Proteins.genome_ID = Genomes.genome_ID WHERE Proteins.genome_ID = ? ORDER BY Proteins.genome_ID,contig,start ";
	$mysql = $dbh->prepare($query);
	foreach my $File (@$faaFiles){
		my $genomeID = &getGenomeID($File);
		
		$mysql->execute($genomeID);
		while(@result_table = $mysql->fetchrow_array()){
		my $string;
		if (defined $result_table[2] && defined $hash_keywords{$result_table[2]}){
		$string = join("\t",@result_table)."\t".$hash_keywords{$result_table[2]}."\n";
		}else{
		$result_table[2] = "-";
		$string = join("\t",@result_table)."\t"."\n";
		}
		print $fh $string;
		}
	}
}

	
close $fh;	
	
	
}	

###########################################################################################
################ SUBROUTINES ##############################################################
################ FOR PROCESSING  ##########################################################
###########################################################################################
###
###
###
###
###
sub getIniList{
my $File = shift;
my %Hashtable;
my $fh;


open ($fh,'<',$File) or warn "Couldn't initilize. File is missing\n";
while(<$fh>){
	unless($_ =~ m/^#/){
	chomp($_);
	my $string = $_;
	my $type = (split(/\s+/,$string))[0];
	my $value;
	my $substring = (split("=",$string))[-1];
	unless($substring =~ /\"\"/){$value = (split("\"",$substring))[-1];}
	$Hashtable{$type}=($value);
	}
}
close $fh;

return \%Hashtable;

}

sub checkDBIconnection{
if (!$dbh->ping){
$dbh = $dbh->clone() or die "cannot connect to db";
}
}

sub getGenomeID{
#gets a complete File pathway
my $file = shift;
my ($name,$dir,$ext) = fileparse($file, qr/\.[^.]*/);

#print "\nFound $name as File\n";
if($name =~ s/\.genes\.faa//){
return $name;} #JGI
elsif($name =~ /(GC(A|F){1}\_\d+).*/){
return $1;} #NCBI Genbank or RefSeq
elsif($name =~ s/\.\d\_\D+.+//){
return $name;} #NCBI unclear name
elsif(length $name > 32){
return substr $name,32;}
else{#print "\tConsidering $name as genomeID\n";
my ($name2,$dir,$ext) = fileparse($name, qr/\.[^.]*/);
return $name2};

}

sub getallFiles{
#Get Directory and Pattern for Files to search for recursively
my $directory = shift;
my @patterns = @_;
my @files;
my $pattern = join("|",@patterns);
	find ( sub {
		
		return unless -f;       #Must be a file
		return unless /$pattern/;  #Must match regex
		push @files, $File::Find::name;
	}, $directory );
	
return \@files;
}

sub getallDirs{
#Get Directory and Pattern for Files to search for recursively
#ATTENTION: Also returns Files if extension fits
my $directory = shift;
my $pattern = shift;
my @files;

	find ( sub {
	  # return unless -d;       #Must be a file
	    
	    return unless /$pattern/;  #Must match regex
	    push @files, $File::Find::name;
	}, $directory );
	
return \@files;
}

sub warning{
print "Warning: ",shift,"\n";
}

sub help{
my $message = <<'HELP_MESSAGE';
HMSSS version 1.0

Usage: perl HMSSS.pl 

The options can be defined in the Ini.txt file, but will be overwritten if a corresponding option is specified on the command line. Database connection options are required for the process and can be specified either in the Ini.txt file or on the command line.

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
   -redoreports	Skip hmmsearch only use results from existing report files
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

HELP_MESSAGE

print $message;
}
