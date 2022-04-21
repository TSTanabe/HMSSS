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

use Parallel::ForkManager;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::SearchIO::FastHitEventBuilder;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::DB::Taxonomy;


#my $date = localtime();
#&printline("$date");
my %opts;
getopts('d:f:q:u:p:b:r:m:', \%opts); #g for new genomes // m for PSI-BLAST the current dataset
#print "\n\n",$opts{u}," ",$opts{b}," ",$opts{p}," ",$opts{d}," ", $opts{f}," ",$opts{q},"\n";

my ($db_user, $db_name, $db_pass) = ($opts{u},$opts{b}, $opts{p});
my $dbh = DBI->connect("DBI:mysql:database=$db_name", $db_user, $db_pass) or die "\nERROR:\tDatabase connection failed!\n $! ";

##################################################################
################ MAIN BODY #######################################
##################################################################
#
#
#
my $Phylogroup = $opts{q};
my $Specgroup = $opts{r};

my $distance = $opts{d};
my $PatternFile = $opts{f};
my $query = "SELECT DISTINCT Proteins.genome_ID FROM Proteins LEFT JOIN Genomes ON Proteins.genome_ID = Genomes.genome_ID WHERE cluster_ID IS NULL ";
if($Specgroup){$query .= "AND ".$Phylogroup." = \"".$Specgroup."\";"}
else{$query.=";"}
#print $query,"\n";
my $modus = $opts{m};
if($modus eq "synteny"){
print "\n";
print "#"x 40;
print "\n# Starting gene cluster prediction\n";
print "#"x 40;
print "\n";
&Clusterprediction($distance,$query);
}elsif($modus eq "keyword"){
print "\n";
print "#"x 40;
print "\n# Starting keyword assignment\n";
print "#"x 40;
print "\n";
&Clusteranalyser($PatternFile,$query);
}



##################################################################
################ SUBROUTINES BODY ################################
##################################################################
#
#
#
sub Clusterprediction{
my $distance = shift;
my $query = shift;
my $mysql;
my $genomes;
#print $query;
$mysql = $dbh->prepare($query);
#$mysql = $dbh->prepare('SELECT DISTINCT genome_ID FROM Proteins WHERE cluster_ID IS NULL;');
my $limit = $mysql->execute();
my $count = 1;
my $fork = Parallel::ForkManager->new(40);

while($genomes = $mysql->fetchall_arrayref(undef,500)){
	foreach my $genome( @$genomes){
	print "Genome $count of $limit\t Searching for gene clusters from genomic locations in ",$genome->[0],"\n";$count++;	
	#Hier die Fork starten falls nötig
	$fork->start and next;
	my $genI;
	my $genes;
	my $sql;
	
	my $cflag = 1;
	my $clusterindex = 1;
	my $genomeID = $genome->[0];
	my %Clusters;
	
	my $dbhf = $dbh->clone();

	#print "Searching for gene clusters from genomic locations in ",$genomeID,"\n";
	$sql = $dbhf->prepare('DELETE FROM Clusters WHERE genome_ID= ? ;');
	$sql->execute($genomeID);
		
	$sql = $dbhf->prepare("SELECT protein_ID, contig, start, end FROM Proteins WHERE genome_ID = ? AND NOT contig IS NULL AND NOT start IS NULL AND NOT end IS NULL ORDER BY contig, start;");
	$genI = $sql->execute($genomeID);
	$genI = $genI-1;
	
		while($genes = $sql->fetchall_arrayref(undef,5000)){
		$Clusters{$genes->[0]->[0]} = $genomeID."_".$clusterindex;
			for(my $i=1; $i<=$genI;$i++){
				#if contig1 eq contig0
				if($genes->[$i]->[1] eq $genes->[$i-1]->[1]
				   && ($genes->[$i]->[2]-$distance)<$genes->[$i-1]->[2]){

				$Clusters{$genes->[$i]->[0]} = $genomeID."_".$clusterindex;
				$cflag = 0;
				}else{
					if($cflag){
					delete $Clusters{$genes->[$i-1]->[0]};
					}
				$clusterindex++;
				$cflag = 1;
				}
			$Clusters{$genes->[$i]->[0]} = $genomeID."_".$clusterindex;
			}
			
			
		}
	$dbhf = $dbh->clone();
	#Write Clusters keys into Cluster table 
	my $clustertable = $dbhf->prepare('INSERT IGNORE INTO Clusters(cluster_ID,genome_ID) VALUES(?,?)');
	keys %Clusters;
	while(my($protID, $clusterID) = each %Clusters) {
		$clustertable->execute($clusterID,$genomeID);
	}
	
	#Write Clusters into Proteins table
	my $proteintable = $dbhf->prepare('UPDATE Proteins SET cluster_ID = ? WHERE protein_ID = ?');
	keys %Clusters;
	while(my($protID, $clusterID) = each %Clusters) {
		$proteintable->execute($clusterID,$protID);
	}
	
	undef(%Clusters);
	$dbhf->disconnect();
	$fork->finish;	
	}
	$fork->wait_all_children;
}
&checkDBIconnection();
}


sub Clusteranalyser{
my $pattern = shift;
my $query = shift;
#print $query;
my $mysql;
#Load the patterns
#
my @pattern_size_array;
my @pattern_clust_array;
my @pattern_key_array;
open my $fh, "<", $pattern or die "\nERROR:\t$pattern does not exist\n",$!;
my $line;
	while ( $line = <$fh> ) {
	chomp $line; 
	my ($size,$name, $rest) = split /\t/, $line, 3;
	push(@pattern_size_array,$size);
	push(@pattern_key_array,$name);
	push(@pattern_clust_array,[split(/\t/, $rest)]) ;
	
	}
close $fh;

#cluster_IDs selektieren
#$mysql = $dbh->prepare($query);
$mysql = $dbh->prepare('SELECT DISTINCT cluster_ID FROM Proteins WHERE NOT cluster_ID IS NULL;');
my $sumCounter = $mysql->execute();
my $counter = 1;

while(my @clusterID = $mysql->fetchrow_array()){
print " Cluster $counter of $sumCounter \t", $clusterID[0], " "; $counter++;


	#Returns a ref for an array full of refs to the HMM
	my @cluster;
	my $sql = $dbh->prepare('SELECT DISTINCT HMM from Proteins where cluster_ID = ? ');
	$sql->execute($clusterID[0]);
	my  $HMMs = $sql->fetchall_arrayref();
	foreach my $i(0..scalar(@{$HMMs}-1)){
		
		my $gene = $HMMs->[$i]->[0] ;
		$gene = (split(/\_/, $gene))[-1];
		push(@cluster,$gene);
	}
	print join("\t",@cluster),"\n";

	#iterates the HMMs of a single cluster

	my %asskeywords;
	foreach my $j(0..$#pattern_size_array){
	my $threshold = $pattern_size_array[$j];
	my $pattern = $pattern_clust_array[$j];
	my $key = $pattern_key_array[$j];
	#print "K $key Thrs $threshold\t";
	#print "\nKey is $key\t Threshold is: $threshold\t Pattern is \t",join("\t",@{$pattern_clust_hash{$key}}),"\n";
	#print "\nCluster is \t",join("\t",@cluster),"\n";
		foreach my $p (@$pattern){
		#print "Searching for $p thrs $threshold\n";
		if ( grep( /$p/, @cluster) ) {
			$threshold--;
			
			if($threshold==0){
			$asskeywords{$key} = $clusterID[0];
			}
		}
		}
	}
	
	# assign keywords to a single cluster
	keys %asskeywords; # reset the internal iterator so a prior each() doesn't affect the loop
	print "\t\t\tkeywords";
	while(my($k, $v) = each %asskeywords){
	print " $k ";
	&assignkeyword($v,$k);
	}
	print "\n";
}
}
sub checkDBIconnection{
if (!$dbh->ping){
$dbh = $dbh->clone() or die "cannot connect to db";
}
}

sub assignkeyword{
my $cluster = shift;
my $key = shift;

my $mysql = $dbh->prepare("INSERT INTO Functions (cluster_ID, keyword) SELECT ?,? FROM DUAL WHERE NOT EXISTS (SELECT * FROM Functions WHERE cluster_ID= ? AND keyword = ? LIMIT 1);");
$mysql->execute($cluster,$key,$cluster,$key);

#my $mysql = $dbh->prepare("INSERT INTO Functions (cluster_ID, keyword) VALUES (?,?);");
#$mysql->execute($cluster,$key);

}
