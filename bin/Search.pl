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

use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError);

use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::SearchIO::FastHitEventBuilder;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::DB::Taxonomy;


my $date = localtime();
#print("$date\n");
my %opts;
getopts('u:p:b:a:c:t:v:q:f:g:n:r:l:o:k:w:y:', \%opts);
#print "\n\n",$opts{u}," ",$opts{b}," ",$opts{p}," ",$opts{a}," ", $opts{c}," ",$opts{t},"\n";

my ($db_user, $db_name, $db_pass) = ($opts{u},$opts{b}, $opts{p});
my $dbh = DBI->connect("DBI:mysql:database=$db_name", $db_user, $db_pass) or die "\nERROR:\tDatabase connection failed! $db_user, $db_name, $db_pass\n $! ";

##################################################################
################ MAIN BODY #######################################
##################################################################
#
#
#

#Globale Variablen
my $search_algorithm = $opts{a};
my $cutoff_type = $opts{c};
my $cutoff_value = $opts{v};
my $thrs_path = $opts{t};
my $query_path = $opts{q};
my $fasta_dir = $opts{f}; #Fasta File storage directory
my $gff_dir = $opts{g}; #gff3 File storage directory
my $cores = $opts{n}; #cores used by hmmer and other
my $redo_search = $opts{r}; #either 1 or undefined, sometimes 0 // Ignore previous results already in database
my $clean_reports = $opts{l}; #same // deletes the hmm reports after the hmmsearch
my $prodigal = $opts{o}; #same // Prediction of ORFs before hmmsearch
my $redoreport = $opts{k}; #// Read only existing reports no hmmsearch
my $spec_model = $opts{w};
my $spec_thrs = $opts{y};

unless(-d $gff_dir){
$gff_dir = $fasta_dir;
}

&start_text();
unless(-d $fasta_dir){
warn "ERROR:\to such directory found. Skip protein search\n";exit;
}

if($prodigal == 1){
#Overwrite GFF directory
$gff_dir = $fasta_dir;
#Prodigal for ORF search
&translation();
&search();

}elsif($redoreport == 1){
#Same as search without hmmsearch
&start_text();
&redoreports();
}elsif($spec_model && $spec_thrs){
#single model threshold modification
&reiterateHMM();
}else{
#HMMsearch the fasta files and get the gff
&search();
}
##################################################################
################ SUBROUTINES BODY ################################
##################################################################
#
#
#

sub translation{
#gather all fna files / fna.gz
my $Files = &getallFiles($fasta_dir,"\.fa\$","\.fna\$");
my $gzFiles = &getallFiles($fasta_dir,"\.fa\.gz","\.fna\.gz");

my $count = 1;
my $sum = scalar @$Files + scalar @$gzFiles;
#Uncompressed File Translation
	foreach my $File (@$Files){
	print "\n\nTranslate uncompressed assembly $count of $sum\n"; $count++;
	my ($name,$dir,$ext) = fileparse($File, qr/\.[^.]*/);

	my $PRODIGAL = "prodigal -a $dir$name.faa -f gff -i $File -o $dir$name.features";
	system($PRODIGAL);

	&FaaToGff("$dir$name.faa");	#Custom Gff3 File writer and packer
	&packgz("$dir$name.faa");	#Pack fasta File
	&packgz("$dir$name.features");
	unlink("$dir$name.faa");
	unlink("$dir$name.features");

	}

#Compressed File Translation
	foreach my $gzFile (@$gzFiles){
	&unpackgz($gzFile);
	print "\n\nTranslate compressed assembly $count of $sum\n"; $count++;
	my $File = $gzFile =~ s/\.gz$//r;
	my ($name,$dir,$ext) = fileparse($File, qr/\.[^.]*/);

	my $PRODIGAL = "prodigal -a $dir$name.faa -f gff -i $File -o $dir$name.features";
	system($PRODIGAL);

	&FaaToGff("$dir$name.faa");	#Custom Gff3 File writer and packer
	&packgz("$dir$name.faa");	#Pack fasta File
	&packgz("$dir$name.features");
	unlink($File);
	unlink("$dir$name.faa");
	unlink("$dir$name.features");

	}

}

sub search{

$date = localtime();
print("$date\n");

#Ready Thresholds
my $Threshold_hashref;
if($cutoff_type eq 'standard_cutoff'){
$Threshold_hashref = &makeThresholdHash($thrs_path);
}

#Unless reiteration is requested select existing genome_IDs for omitting them in the search
my $genomeIDhash_ref;
unless($redo_search == 1){
my $mysql = $dbh->prepare("SELECT distinct genome_ID,genome_ID from Proteins ");
$mysql->execute();
$genomeIDhash_ref = $mysql->fetchall_hashref('genome_ID');
}


#covert plain fasta to compressed
my $plainfaaFiles = &getallDirs($fasta_dir,".faa\$");
foreach my $File (@$plainfaaFiles){
	&packgz($File);
	my ($name,$dir,$ext) = fileparse($File, qr/\.[^.]*/);
	if(-e "$dir$name\.gff"){&packgz("$dir$name\.gff");}
}

#Get all faa.gz files in the fasta directory
my $faaFiles = &getallDirs($fasta_dir,".faa\.gz");
my $size = scalar @{ $faaFiles };
my $count = 1;
if($size == 0){print "\tWarning: There where no files searched. Consider using -prodigal if there are nucleotide fasta files, or -redosearch if files should be searched again\n"}
foreach my $File (@$faaFiles){
	
	
	my $genomeID = &getGenomeID($File); #get the genomeID from the File name
	print "$count of $size $genomeID \n";$count++;
	#print $genomeID," is present? ",$genomeIDhash_ref->{$genomeID}->{'genome_ID'},"\n";
	#Omitting this File if already in database
	if(defined $genomeIDhash_ref->{$genomeID}){next;} #Omitting files already in the local DB

	my $faaFile = &unpackgz($File);
	my ($name,$dir,$ext) = fileparse($faaFile, qr/\.[^.]*/);


	#############################
	#Search the Files with HMMER#
	#############################
	print "\tStarting the ";
	if($search_algorithm eq 'hmmer'){
	print "HMM search\n";
	&HMMsearch($faaFile,$genomeID,$name,$dir,$query_path,$cores,$Threshold_hashref,$cutoff_value);
	}elsif($search_algorithm eq 'blastp'){
	print "Blast search\n";
	&BlastSearch($faaFile,$genomeID,$name,$dir,$query_path,$cores,$Threshold_hashref,$cutoff_value);
	}
	unlink($faaFile);

	##########################
	#Get the genomic features#
	##########################
	
	print "\tUpdating genomic features for $genomeID\n";
	my $gffFiles = &getallFiles($gff_dir,"$genomeID.+gff.gz");
	if(defined $gffFiles->[0]){
	
	my $gffFile = &unpackgz($gffFiles->[0]);
	&GenomicGeneralFeatures($gffFile,$genomeID);
	unlink($gffFile);
	}

}




}



sub redoreports{
#redo analysis of all reports with the full threshold set

my $thrs_hashref;
if($cutoff_type eq 'standard_cutoff'){
$thrs_hashref = &makeThresholdHash($thrs_path);
}
my $faaFiles = &getallDirs($fasta_dir,".HmmReport");
my $size = scalar @{ $faaFiles };
my $count = 1;

	foreach my $File (@$faaFiles){

		my $genomeID = &getGenomeID($File);
		print "$count of $size $genomeID \n";$count++;

		my ($name,$dir,$ext) = fileparse($File, qr/\.[^.]*/);
		
		&HMMreports($File,$genomeID,$name,$dir,$cores,$thrs_hashref,$cutoff_value);
		
		my $gffFiles = &getallFiles($gff_dir,"$genomeID.+gff.gz");
		if(defined $gffFiles->[0]){
		my $gffFile = &unpackgz($gffFiles->[0]);
		&GenomicGeneralFeatures($gffFile,$genomeID);
		unlink($gffFile);

		}

	}
}



sub reiterateHMM{
#redo HMMreports with a single modified model and score
my $directory = $fasta_dir;
my $HMM = $spec_model;
my $score = $spec_thrs;

my $ReportFiles = &getallDirs($directory,".HmmReport\$");


#fork of the process here possible?
#iterate all reports
my $count = 1;
my $size = scalar @{$ReportFiles};
foreach my $report(@$ReportFiles){
print "$count of $size\n";$count++;
print "\tReport $report\n";
my ($hitname,$hitscore,$hitsignificance);
my %protein_score;
my %protein_seq;
my $in;

#STEP 1 Get significant results for specific HMM
$in = Bio::SearchIO->new(-format=>'hmmer', -version => 3, -file  => $report);
$in->attach_EventHandler(Bio::SearchIO::FastHitEventBuilder->new);

while( my $result = $in->next_result ){
   unless($result->query_name() eq $HMM){next;}

     while(my $hit = $result->next_hit){
     $hitname = $hit->name();
     $hitscore = $hit->score();
     $hitsignificance = $hit->significance();
       if($score < $hitscore){
       $protein_score{$hitname} = $hitscore;
       #$protein_evalue{$hitname} = $hitsignificance;   
       }
     }


   }

   #STEP 2 Get the Sequences
   if(keys %protein_score>0){
   print "\tFound ".%protein_score." new hits\n";
   my $genomeID = &getGenomeID($report);
   my ($name,$dir,$ext) = fileparse($report, qr/\.[^.]*/);
   my $File = unpackgz($dir.$name.".faa.gz");

   my $seq;
   $in = Bio::SeqIO->new(-file => $File, -format => "fasta" );
     while($seq = $in->next_seq){
       if($protein_score{($seq->id)}){
       $protein_seq{($seq->id)}=$seq->seq;
       }
     
     }
   unlink($File);

   #STEP 3 Write to database
   my $mysql;
   #$mysql = $dbh->prepare('INSERT IGNORE INTO Proteins(protein_ID,genome_ID,HMM,score,sequence) VALUES(?,?,?,?,?)');
   $mysql = $dbh->prepare('INSERT INTO Proteins(protein_ID, genome_ID,HMM,score,sequence) VALUES(?,?,?,?,?) ON DUPLICATE KEY UPDATE HMM = IF(score < VALUES(score),VALUES(HMM),HMM), score =  IF(score < VALUES(score),VALUES(score),score);');
     while(my($protID, $hitscore) = each %protein_score) {
     my $identifier = $protID."-".$genomeID;
     print "\t$identifier $hitscore\n";
     if(length$identifier>32){$identifier = substr $identifier, 0,32;}
     print "\tExecute: ";
     print $mysql->execute($identifier,$genomeID,$HMM,$hitscore,$protein_seq{$protID});
     print"\n";
     }
   
   
   #STEP 4 Get the gff3 features
   
   my $gffFiles = &getallFiles($gff_dir,"$genomeID.+gff.gz");
		if(defined $gffFiles->[0]){
		my $gffFile = &unpackgz($gffFiles->[0]);
		&GenomicGeneralFeatures($gffFile,$genomeID);
	unlink($gffFile);
   
   }
   }else{print"\t No new hits\n";}
   

}
#fork off to here possible?


}
##################################################################
################ TERTIÄRE SUBROUTINEN ############################
##################################################################
#
#
#

sub HMMsearch{
my $File = shift;
my $genomeID = shift;
my $name = shift;
my $dir = shift;
my $HMMlib = shift;
my $cores = shift;
my $thrs_hashref = shift;
my $thrs_score = shift;
my ($HMMER,$status);


	#my ($name,$dir,$ext) = fileparse($File, qr/\.[^.]*/);
	my $outputFile = $dir.$name."\.HmmReport";
	
	print "\tExecuting hmmsearch for $name \n";
	
	$HMMER = "hmmsearch -E 0.0001 --cpu $cores $HMMlib $File > $outputFile";
	$status = system($HMMER);
		if($status > 0){
		warn "\nERROR:\t$HMMER\nDied with exit code:\t$status\n";
		next;
		}
		
	
	my $predHMM;
	my ($hitname,$hitscore,$hitevalue);
	my %protein_hmm;
	my %protein_score;
	my %protein_evalue;
	my %protein_seq;
	my $in;
	
	$in = Bio::SearchIO->new(-format=>'hmmer', -version => 3, -file   => $outputFile);
$in->attach_EventHandler(Bio::SearchIO::FastHitEventBuilder->new);
	
	print "\tExecuting assigning motifs for $name \n";
	
	while( my $result = $in->next_result ){
	$predHMM = $result->query_name();

		while(my $hit = $result->next_hit){
		$hitname = $hit->name();
		$hitscore = $hit->score();
		$hitevalue = $hit->significance();
		if((defined $thrs_hashref->{$predHMM} && $thrs_hashref->{$predHMM} < $hitscore) || ((!(defined $thrs_hashref->{$predHMM})) && $thrs_score < $hitscore) ){
			unless(defined $protein_score{$hitname}){
				$protein_hmm{$hitname} = $predHMM;
				$protein_score{$hitname} = $hitscore;
				$protein_evalue{$hitname} = $hitevalue;
			
			}elsif($protein_score{$hitname} < $hitscore){
				$protein_hmm{$hitname} = $predHMM;
				$protein_score{$hitname} = $hitscore;
				$protein_evalue{$hitname} = $hitevalue;
			}			
		}
		}
	
	}
	
#STEP 2 Get the Sequences
#This solution cannot be forked but now no need for backtranslate to fasta file from report
	
	print "\tRetrieve sequences for $name \n";
	
	my $seq;
	$in = Bio::SeqIO->new(-file => $File, -format => "fasta" );

	while($seq = $in->next_seq){
		if($protein_hmm{($seq->id)}){
			$protein_seq{($seq->id)}=$seq->seq;
		}
		
	}

#STEP 3 Get Genome ID
	my $mysql;

	&checkDBIconnection();
	$mysql = $dbh->prepare('INSERT IGNORE INTO Genomes(genome_ID) VALUES (?)');
	$mysql->execute($genomeID);
		
#STEP 4 Write to database

	print "\tWriting results to DB for $name \n";
	
	$mysql = $dbh->prepare('INSERT IGNORE INTO Proteins(protein_ID, genome_ID,HMM,score,sequence) VALUES(?,?,?,?,?) ON DUPLICATE KEY UPDATE HMM = IF(score < VALUES(score),VALUES(HMM),HMM), score =  IF(score < VALUES(score),VALUES(score),score);');
	while(my($protID, $HMM) = each %protein_hmm) {
	#if($protein_evalue{$protID} != 0){
	#my $evalue = $protein_evalue{$protID};
	#}else{$evalue = 0}
	 my $identifier = $protID."-".$genomeID;
    	 if(length$identifier>32){$identifier = substr $identifier, 0,32;}
		$mysql->execute($identifier,$genomeID,$HMM,$protein_score{$protID},$protein_seq{$protID});
	}
	
if($clean_reports){unlink($outputFile)}

}

sub BlastSearch{
my $File = shift;
my $genomeID = shift;
my $name = shift;
my $dir = shift;
my $querylib = shift;
my $cores = shift;
my $thrs_hashref = shift;
my $thrs_score = shift;

my $blast;
my $queryseqio = Bio::SeqIO->new(-file => $querylib, -format => "fasta");


#my ($name,$dir,$ext) = fileparse($File, qr/\.[^.]*/);
my $pred;
my ($hitname,$bits);
my %protein_hmm;
my %protein_score;
my %protein_seq;
my $in;
$blast = Bio::Tools::Run::StandAloneBlastPlus->new(
	-db_data => "$File",
	-create => 1);
if($blast->check_db($blast->db) eq 0){

	print "\tExecuting blastp for $name \n";
	
	$blast->blastp(-query => $queryseqio ,
	#-outfile => 'Blast.report',
	#outformat => 5,
	-method_args => [evalue => 0.001]);
	my $blastreport = new Bio::SearchIO(-file => $blast->blast_out);
	
	while (my $result = $blastreport->next_result){
	$pred = $result->query_name;
	
		while(my $hit = $result->next_hit){
		$hitname = $hit->name;
		$bits =  $hit->bits;
		
		if((defined $thrs_hashref->{$pred} && $thrs_hashref->{$pred} < $bits) || ((!(defined $thrs_hashref->{$pred})) && $thrs_score < $bits) ){
			unless(defined $protein_score{$hitname}){
				$protein_hmm{$hitname} = $pred;
				$protein_score{$hitname} = $bits;
			
			}elsif($protein_score{$hitname} < $bits){
				$protein_hmm{$hitname} = $pred;
				$protein_score{$hitname} = $bits;
			}			
		}
		
		}

	}
	$blast->cleanup;
	
	my $seq;
	$in = Bio::SeqIO->new(-file => $File, -format => "fasta" );

	print "\tRetrieve sequences for $name \n";
	
	while($seq = $in->next_seq){
		if($protein_hmm{($seq->id)}){
			$protein_seq{($seq->id)}=$seq->seq;
		}
		
	}

#STEP 3 Get Genome ID
	my $mysql;
	&checkDBIconnection();
	$mysql = $dbh->prepare('INSERT IGNORE INTO Genomes(genome_ID) VALUES (?)');
	$mysql->execute($genomeID);
		
#STEP 4 Write to database
	
	print "\tWriting results to DB for $name \n";
	
	$mysql = $dbh->prepare('INSERT IGNORE INTO Proteins(protein_ID, genome_ID,HMM,sequence,score) VALUES(?,?,?,?,?) ON DUPLICATE KEY UPDATE HMM = IF(score < VALUES(score),VALUES(HMM),HMM), score =  IF(score < VALUES(score),VALUES(score),score);');
	while(my($protID, $HMM) = each %protein_hmm) {
		my $identifier = $protID."-".$genomeID;
		if(length$identifier>32){$identifier = substr $identifier, 0,32;}
		$mysql->execute($identifier,$genomeID,$HMM,$protein_seq{$protID},$protein_score{$protID});
	}
		
}


}





sub GenomicGeneralFeatures{
my $GffFile = shift;
my $genomeID = shift;

	unless(defined $genomeID && defined $GffFile){next;}
	
	my $proteinsql = $dbh->prepare("SELECT DISTINCT protein_ID from Proteins WHERE genome_ID = ? ");
	$proteinsql->execute($genomeID);
	
	my %gff3;
	while(my @proteinID = $proteinsql->fetchrow_array()){
	
	$gff3{(split("-",$proteinID[0]))[0]} = 1; 
	}
	
	
	my ($row,$fh);
	open $fh, "<",$GffFile or die "$!";

	while($row=<$fh>) {
	chomp($row);
		if($row =~ m/^#/){ 
			next;
		}
		
		if($row =~ /ID=(cds-){0,1}(\S+?);/){
			if($gff3{($2)}){
				$gff3{($2)}=$row;
			}
		}elsif($row =~ /\tcds\t/){
		my @ar = split ("\t",$row);
		$gff3{($ar[0])}=$row;
		}		
		
	}
	$fh->close();
	
	&checkDBIconnection();
	$proteinsql = $dbh->prepare("UPDATE Proteins SET locustag=?,contig=?,start=?,end=?,strand=? WHERE protein_ID rlike ? AND genome_ID = ?"); #TODO untested where condition
	while(my($protID, $GFF) = each %gff3) {
			
		my @gff = split(/\t/,$GFF);
		
		#my $id = $protID."-".$genomeID; #if proteinID too long after concat: truncate
		#if(length($id)>32){substr($id,0,32);}
		#print "proteinID = $protID genomeID = $genomeID \n";
		$proteinsql->execute(
		&getLocusTag($gff3{$protID}),
		$gff[0],
		$gff[3],
		$gff[4],
		$gff[6],
		$protID,
		$genomeID
		)
	}

}


sub HMMreports{
my $File = shift;
my $genomeID = shift;
my $name = shift;
my $dir = shift;
my $cores = shift;
my $thrs_hashref = shift;
my $thrs_score = shift;
my ($HMMER,$status);

	
	my $predHMM;
	my ($hitname,$hitscore,$hitevalue);
	my %protein_hmm;
	my %protein_score;
	my %protein_evalue;
	my %protein_seq;
	my $in;
	
	$in = Bio::SearchIO->new(-format=>'hmmer', -version => 3, -file   => $File);
$in->attach_EventHandler(Bio::SearchIO::FastHitEventBuilder->new);
	
	print "\tExecuting assigning motifs for $name \n";
	
	while( my $result = $in->next_result ){
	$predHMM = $result->query_name();

		while(my $hit = $result->next_hit){
		$hitname = $hit->name();
		$hitscore = $hit->score();
		$hitevalue = $hit->significance();
		if((defined $thrs_hashref->{$predHMM} && $thrs_hashref->{$predHMM} < $hitscore) || ((!(defined $thrs_hashref->{$predHMM})) && $thrs_score < $hitscore) ){
			unless(defined $protein_score{$hitname}){
				$protein_hmm{$hitname} = $predHMM;
				$protein_score{$hitname} = $hitscore;
				$protein_evalue{$hitname} = $hitevalue;
			
			}elsif($protein_score{$hitname} < $hitscore){
				$protein_hmm{$hitname} = $predHMM;
				$protein_score{$hitname} = $hitscore;
				$protein_evalue{$hitname} = $hitevalue;
			}			
		}
		}
	
	}
	
#STEP 2 Get the Sequences
#backtranslate to fasta file from report
	
	print "\tRetrieve sequences for $name \n";
	
	#my $genomeID = &getGenomeID($File);
	#my ($name,$dir,$ext) = fileparse($File, qr/\.[^.]*/);
	my $Fasta = unpackgz($dir.$name."\.faa.gz");

	my $seq;
	$in = Bio::SeqIO->new(-file => $Fasta, -format => "fasta" );
	while($seq = $in->next_seq){
		if($protein_score{($seq->id)}){
		$protein_seq{($seq->id)}=$seq->seq;
		}
     
	}

#STEP 3 Get Genome ID
	my $mysql;

	&checkDBIconnection();
	$mysql = $dbh->prepare('INSERT IGNORE INTO Genomes(genome_ID) VALUES (?)');
	$mysql->execute($genomeID);
		
#STEP 4 Write to database

	print "\tWriting results to DB for $name \n";
	
	$mysql = $dbh->prepare('INSERT IGNORE INTO Proteins(protein_ID, genome_ID,HMM,score,sequence) VALUES(?,?,?,?,?) ON DUPLICATE KEY UPDATE HMM = ?, score = ?');
	while(my($protID, $HMM) = each %protein_hmm) {
	#if($protein_evalue{$protID} != 0){
	#my $evalue = $protein_evalue{$protID};
	#}else{$evalue = 0}
	 my $identifier = $protID."-".$genomeID;
	 if(length$identifier>32){$identifier = substr $identifier, 0,32;}
		$mysql->execute($identifier,$genomeID,$HMM,$protein_score{$protID},$protein_seq{$protID},$HMM,$protein_score{$protID});
	}
	
}

##################################################################
################ TERTIÄRE SUBROUTINEN ############################
##################################################################
#
#
#

sub FaaToGff{
my $File = shift;
open my $in, "<", $File;
#&packgz($File);
$File =~ s/\.faa$/\.gff/;
open my $out, ">", $File;

my $line;

while ($line = <$in>){
	if($line =~ s/^>//){
	#print $line;
	my @ar = split ("#",$line);
	my $contig = (split (/\_{1}\d+\W+$/,$ar[0]))[0];
	#my $contig = $ar[0];
	#$contig =~ s/\_{1}\d+\W+$//;
	chop $ar[0];
	print $out $contig,"\tprodigal\tcds\t",$ar[1],"\t",$ar[2],"\t0.0\t+\t0\tID=cds-",$ar[0],";",$ar[3],"\n";
	}

}


close $out;
close $in;
&packgz($File);
unlink($File); #removes the uncompressed gff File
}






##################################################################
################ TERTIÄRE SUBROUTINEN ############################
##################################################################
#
#
#
sub checkDBIconnection{
if (!$dbh->ping){
$dbh = $dbh->clone() or die "cannot connect to db";
}
}

sub packgz{
	
my $archive = shift;

	#for my $archive(@$archives){

	print "\tPacking ",$archive,"\n";
	gzip $archive => "$archive.gz";
	#}
return $archive;
}

sub unpackgz{
	
my $archive = shift;
	$archive =~ s/\.gz$// ;
	print "\tUnpacking ",$archive,".gz \n";
	gunzip "$archive.gz" => $archive;

return $archive;
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

sub makeThresholdHash{
my $File = shift;
my %Hashtable;
my $fh;
open ($fh,'<',$File) or warn "Couldn't open theshold file\n";
while(<$fh>){
my @array = split("\t",$_);

$Hashtable{$array[0]}=($array[1]);
}
close $fh;

return \%Hashtable;
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

sub getLocusTag{
	
	if($_[0] =~ /old\_locus\_tag=(\S*?)[\n|;]/){
		return $1;	
	}elsif($_[0] =~ /locus\_tag=(\S*?)[\n|;]/){
		return $1;	
	}else{
	return "";
	}
}

sub start_text{

print "\n";
print "#"x 40;
print "\n# Starting hmmsearch\n";
print "#"x 40;
print "\n";
}
