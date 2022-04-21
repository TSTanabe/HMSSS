#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use File::Copy;
use File::Find;
use Getopt::Std;
use Cwd qw(getcwd);
use Error qw(:try);

use DBI;
use DBD::mysql;

use Bio::DB::Taxonomy;



my $date = localtime();
print("$date\n");
my %opts;
getopts('u:p:b:a:c:', \%opts);

my ($db_user, $db_name, $db_pass) = ($opts{u},$opts{b}, $opts{p});
my $dbh = DBI->connect("DBI:mysql:database=$db_name", $db_user, $db_pass) or die "\nERROR:\tDatabase connection failed! $db_user, $db_name, $db_pass\n $! ";

##################################################################
################ MAIN BODY #######################################
##################################################################
#
#
#

#Globale Variablen
my $directory = $opts{a};
my $cwd = $opts{c};


&getGenomesAssembly();
&getGenomesTaxonomy();



##################################################################
################ MAIN BODY #######################################
##################################################################
#
#
#

sub getGenomesAssembly{
#Input directory with genome assembly statistics
my $assFiles = &getallFiles($directory,"assembly_(report|stats).txt");

for my $file (@$assFiles){
print "Assembly stats from $file\n";
my ($genomeID,$Biosample,$Bioproject,$Refseq,$Genbank,$Taxid) = (0,0,0,0,0,0,0);
my $mysql;

$genomeID = &getGenomeID($file); #Consider genome ID from filename


open my $fh, "<",$file or die "$!";
	while(<$fh>){
		#print $_,"\n";
		if($_ =~ /Taxid:\s+(\d+)\s/){$Taxid=$1;}
		elsif($_ =~ /BioSample:\s+(SAMN\d+)\s/){$Biosample=$1;}
		elsif($_ =~ /BioProject:\s+(\w{5}\d+)\s/){$Bioproject=$1;}
		#elsif($_ =~ /Organism name:\s+(.*)\n/){$Organism=$1;}
		elsif($_ =~ /GenBank assembly accession:\s+(GCA\S+)\s/){$Genbank=$1;}
		elsif($_ =~ /RefSeq assembly accession:\s+(GCF\S+)\s/){$Refseq=$1;}
		elsif($Taxid and $Biosample and $Bioproject and $Refseq){last;}
	}
close $fh;

	$mysql = $dbh->prepare("SELECT genome_ID FROM Genomes WHERE genome_ID = ? OR genome_ID = ? OR genome_ID = ? ;");
	$mysql->execute($genomeID,$Refseq,$Genbank);
	
	$genomeID = $mysql->fetchrow_array();
	if($genomeID){
	$mysql = $dbh->prepare('UPDATE IGNORE Genomes SET NCBITaxon = ?, NCBIBioproject = ?, NCBIBiosample = ?, NCBIassembly = ? WHERE genome_ID = ?');
	$mysql->execute($Taxid,$Bioproject,$Biosample,$Genbank,$genomeID);
	}
}
	
}

#################################################################

sub getGenomesTaxonomy{
my ($num,@taxonid);

print "Starting taxonomy search\n";

my $taxonDB = Bio::DB::Taxonomy->new(-source => 'flatfile', -directory => "$cwd/source/taxdump", -nodesfile => "$cwd/source/taxdump/nodes.dmp" , -namesfile => "$cwd/source/taxdump/names.dmp");
#$num = $taxonDB->get_num_taxa();
#print "There are $num taxa stores in the TaxonDB $taxonDB\n";

#Get all taxid without Phylum in the ProjectRed DB
my $mysql = $dbh->prepare("SELECT NCBITaxon, genome_ID FROM Genomes WHERE Phylum is null;");
$mysql->execute();

#For all results 
	while(@taxonid = $mysql->fetchrow_array()){
	my %lineage;
	my $taxonid = $taxonid[0];
	my $genomeID = $taxonid[1];
	my $taxon;
	print "Retrieving lineage for taxon id $taxonid\n";
	#Retrieve lineage of Phylogeny with rank from flatfileDB
		do{
		
		$taxon = $taxonDB->get_taxon(-taxonid => $taxonid);
		if(defined $taxon){
		$lineage{$taxon->rank} = $taxon->scientific_name;
		$taxonid = $taxon->parent_id;
		}else{next;}
		}until($taxon->rank eq "phylum");

	#Write correct ranks into Database
	#	strain species genus family order class phylum
	
	
	my $species = $lineage{'species'};
	unless($lineage{"strain"}){
	$lineage{"strain"} = "NULL";
	}elsif($species){
	$lineage{"strain"} =~ s/$species //g;
	}

		
	my $sqlphylogeny = $dbh->prepare('UPDATE Genomes SET strain = ?, Species = ?, Genus = ?, Family = ?, Ordnung = ?, Class = ?, Phylum = ? WHERE genome_ID = ?');
	$sqlphylogeny->execute($lineage{'strain'},$lineage{'species'},$lineage{'genus'},$lineage{'family'},$lineage{'order'},$lineage{'class'},$lineage{'phylum'},$genomeID);
	
	}
print "Finished taxon assignment\n";

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
