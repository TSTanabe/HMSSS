#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Tk;
use Tk::Optionmenu;
use Tk::ROText;
use Tk::BrowseEntry;
use Tk::Pane;
use Tk::PNG;

use File::Basename;
use File::Copy;
use File::Find;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

use Parallel::ForkManager;

use DBI;
use DBD::mysql;

use Cwd qw(getcwd);

use Bio::SeqIO;
use Bio::SearchIO;
use Bio::SearchIO::FastHitEventBuilder;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::DB::Taxonomy;


#### FRAMES ####
#### DONE Frame größen Verhalten anpassen an konstanten rahmen
#### Labels in bessere Beschreibung ändern
#### Rahmen weglassen oder verkleinern wo unnötig
#### Hilfsverbindung zu Datenbanken löschen
#### TODO Sortierung und Limit bei Sequenzlänge mysql
#### TODO synteny und collinearity unterscheiden
#### Ausgabe als Tabelle ohne Sequenz

#### DONE Reihenfolge der HMMER/BLAST suche ändern, damit besonders große Datensets per genom verarbeitet werden
#### DONE Genome die bereits in der Datenbank sind bei suche ignorieren außer es wird extra angefordert
#### Das cutoff Problem Lösung entwickeln. Entweder Reports reiterierbar machen für bestimmte modelle oder generell mit noise und trusted cutoffs arbeiten => Gelöst durch reiterierbare hmmreports
#### TODO Modelle bei Library kompilierung auswählbar machen um eingeschränkte suchen zu ermöglichen


my $cwd = getcwd;
my $dbh;
my ($db_user,$db_name,$db_pass);

my $mw = MainWindow->new(title => "HMSSS", -bg => 'white');



$mw->optionAdd('*font' => 'sans 11');
$mw->geometry(950 . 'x' .550);
$mw->update;



my $menu_frame = $mw->Frame(-relief => 'ridge', -borderwidth => 3);
$menu_frame->pack(-side => 'left', -fill => 'y', -expand => 0, -anchor => 'nw');

my $work_frame = $mw->Scrolled('Frame', -scrollbars => 'se',-sticky => 'wens', -relief => 'flat', -width => 750, -borderwidth => 1);
$work_frame->packPropagate(0);
$work_frame->pack(-side => 'right', -fill => 'both', -expand => 1, -anchor => 'nw');
	

my $ini_hashref = &getIniList("$cwd/HMSSS.ini"); #Hashtable for standard parameters

###	####Menu Frames loading####
my $dbcon_frame = $work_frame->Frame(-relief => 'flat', -width => 1000, -borderwidth => 3);
my $search_frame = $work_frame->Frame(-relief => 'flat', -width => 1000, -borderwidth => 3);
my $cluster_frame = $work_frame->Frame(-relief => 'flat', -width => 1000, -borderwidth => 3);
my $mysql_frame = $work_frame->Frame( -relief => 'ridge', -width => 1000, -borderwidth => 3);
my $result_frame = $work_frame->Frame(-relief => 'ridge', -width => 1000, -borderwidth => 3);

&menuframe();
&dbconframe($ini_hashref);
&searchframe($ini_hashref);
&clusterframe($ini_hashref);
&outputqueryframe();







###########################################################################################
################ MAINROUTINE ##############################################################
################ FOR GUI     ##############################################################
###########################################################################################
###

$mw->MainLoop();

exit(0);

	
###########################################################################################
################ SUBROUTINES ##############################################################
################ FOR GUI WIDGETS  #########################################################
###########################################################################################
###
###
#Widgets Subroutinen laden ohne pack, auf button press pack forget and pack. Sonst werden widgets neu erstellt.
###
###
###	
###	####Sub Menu 0 Frame loading####
#	#### Menu Frame	     ####
#
#
#
sub menuframe{

	$menu_frame->Button(-text => 'Connect', -relief => 'raised', -width => 12, -height => 4, -command => sub {
		&clearWorkframe();
		$dbcon_frame->pack(-side => 'top', -fill => 'none', -expand => 1, -anchor => 'nw');
	})->pack(-fill=>'both',-expand => 0 );
	
	
	$menu_frame->Button(-text => 'Search', -relief => 'raised', -width => 12, -height => 4, -command => sub {
		&clearWorkframe();
		$search_frame->pack(-side => 'left', -fill => 'none', -expand => 1, -anchor => 'nw');

	})->pack(-fill=>'both',-expand => 0);
	
	
	$menu_frame->Button(-text => 'Synteny', -relief => 'raised', -width => 12, -height => 4, -command => sub {
		&clearWorkframe();
		$cluster_frame->pack(-side => 'left', -fill => 'none', -expand => 1, -anchor => 'nw');
	})->pack(-fill=>'both',-expand => 0);
	
	
	$menu_frame->Button(-text => 'Request results', -relief => 'raised', -width => 12, -height => 4, -command => sub {
		&clearWorkframe();
		$mysql_frame->pack(-side => 'left', -fill => 'both', -expand => 1, -anchor => 'nw');

	})->pack(-fill=>'both',-expand => 0);


	$menu_frame->Button(-text => 'Show results', -relief => 'raised', -width => 12, -height => 4, -command => sub {
		&clearWorkframe();
		$result_frame->pack(-side => 'left', -fill => 'both', -expand => 1, -anchor => 'nw');

	})->pack(-fill=>'both',-expand => 0);
	
}


###	####Sub Menu 1 Frame loading####
#	#### Connect Frame	     ####
#
#
#

sub dbconframe{

my $ini_hashref = shift;
$db_user = $ini_hashref->{Username};
$db_name = $ini_hashref->{Database};	
$db_pass = $ini_hashref->{Password};

if($db_user && $db_name && $db_pass){
$dbh = DBI->connect("DBI:mysql:database=$db_name", $db_user, $db_pass);
}


	$dbcon_frame->Label(-text => "Username:")->grid(-row => 0, -column => 0);
	$dbcon_frame->Entry(
	-width => 20,
	-textvariable => \$db_user
	)->grid(-row => 0, -column => 1, -padx => "1m", -pady => "1m");
	
	$dbcon_frame->Label(-text => "Database:")->grid(-row => 1, -column => 0);
	$dbcon_frame->Entry(
	-width => 20,
	-textvariable => \$db_name
	)->grid(-row => 1, -column => 1, -padx => "1m", -pady => "1m");

	
	$dbcon_frame->Label(-text => "Password:")->grid(-row => 2, -column => 0);
	$dbcon_frame->Entry(
	-show => '*',
	-width => 20,
	-textvariable => \$db_pass
	)->grid(-row => 2, -column => 1, -padx => "1m", -pady => "1m");	
#####CONNECT BUTTON
my $connectbuttontext = "Connect";
	$dbcon_frame->Button(
	-textvariable => \$connectbuttontext,
	-width => 64,
	-command => 
	sub{
	 $connectbuttontext = "Connected";
	 $dbh = DBI->connect("DBI:mysql:database=$db_name", $db_user, $db_pass) or $connectbuttontext = "Connection failed ";
	 }
	)->grid(-row => 3, -column => 0, -columnspan => 2, -padx => "2m", -pady => "2m");

}

########################################	
###	####Sub Menu 2 Frame loading####
#	#### Search Frame	     ####
#
#
#

sub searchframe{
#Initilize default directories & files
my $ini_hashref = shift;
my $thrs_path = "$cwd/Thresholds";
if($ini_hashref->{Threshold_File}){$thrs_path = $ini_hashref->{Threshold_File}}
my $query_path = "$cwd/HMMlib";
if($ini_hashref->{Query_File}){$query_path = $ini_hashref->{Query_File}}
my $fasta_dir = "$cwd/Genomes";
if($ini_hashref->{Fasta_directory}){$fasta_dir = $ini_hashref->{Fasta_directory}}
my $gff_dir = "$cwd/Genomes";
if($ini_hashref->{GFF_directory}){$gff_dir = $ini_hashref->{GFF_directory}}
my $cores = 4;
if($ini_hashref->{Cores}){$cores = $ini_hashref->{Cores}}
my $ass_dir = "$cwd/Genomes";
if($ini_hashref->{Assemblystats_directory}){$ass_dir = $ini_hashref->{Assemblystats_directory}}


#Set parameters for hmmsearch
my $search_hmmsearch_frame = $search_frame->Frame(-relief => 'ridge', -borderwidth => 2);
$search_hmmsearch_frame->Label(-text => "Set search parameters")->grid(-row => 0, -column => 0, -columnspan => 3, -padx => "1m", -pady => "1m", -sticky => 'w');

my ($search_algorithm,$tsearch_algorithm);
	$search_hmmsearch_frame->Label(-text => "Search algorithm")->grid(-row => 1, -column => 0, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	my $algorithm_option = $search_hmmsearch_frame->Optionmenu(
	 -options => [[HMMER=>'hmmer'], [BLASTp=>'blastp']],
	 -command => sub {}, #Fires the selection event
	 -variable => \$search_algorithm,				   #Value of the selection
	 -textvariable => \$tsearch_algorithm,		   #Selected Text
	 -width => 20
	)->grid(-row => 1, -column => 1, -padx => "1m", -pady => "1m");


	
my ($cutoff_value, $cutoff_type,$tcutoff_type);
my ($thrsButton,$FileEntry,$thrsEntry);
	$search_hmmsearch_frame->Label(-text => "Cutoff")->grid(-row => 2, -column => 0, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	$thrsEntry = $search_hmmsearch_frame->Entry(
	-width => 40,
	-textvariable => \$cutoff_value
	)->grid(-row => 2, -column => 2, -padx => "1m", -pady => "1m");


	$search_hmmsearch_frame->Label(-text => "Threshold File")->grid(-row => 3, -column => 0, -padx => "1m", -pady => "1m", -sticky => 'w');
	$thrsButton = $search_hmmsearch_frame->Button(
	-text => "Select file",
	-width => 21,
	-command => sub { $thrs_path = $mw->getOpenFile();}
	)->grid(-row => 3, -column => 1, -padx => "1m", -pady => "1m");
	
	$FileEntry = $search_hmmsearch_frame->Entry(
	-width => 40,
	-textvariable => \$thrs_path
	)->grid(-row => 3, -column => 2, -padx => "1m", -pady => "1m");

	$search_hmmsearch_frame->Optionmenu(
	 -options => [ [List=>'standard_cutoff'], 
	 		#[Evalue=>'evalue_cutoff'],
	 		[Score=>'score_cutoff']],
	 -command => sub {
		 if($cutoff_type eq 'standard_cutoff'){
		 $cutoff_value = "10";
		 $thrsEntry->configure(-state => 'disabled');
		 $thrsButton->configure(-state => 'normal');
		 $FileEntry->configure(-state => 'normal');
		 }elsif($cutoff_type eq 'evalue_cutoff'){
		 $cutoff_value = '1e-2';
		 $thrsEntry->configure(-state => 'normal');
		 $thrsButton->configure(-state => 'disabled');
		 $FileEntry->configure(-state => 'disabled');
		 }elsif($cutoff_type eq 'score_cutoff'){
		 $cutoff_value = '10';
		 $thrsEntry->configure(-state => 'normal');
		 $thrsButton->configure(-state => 'disabled');
		 $FileEntry->configure(-state => 'disabled');
		 }
	 },
	 -variable => \$cutoff_type,				   #Value of the selection
	 -textvariable => \$tcutoff_type,			   #Selected Text
	 -width => 20
	)->grid(-row => 2, -column => 1, -padx => "1m", -pady => "1m");
	
		
	$search_hmmsearch_frame->Label(-text => "Fasta/HMM Library")->grid(-row => 4, -column => 0, -padx => "1m", -pady => "1m", -sticky => 'w');
	$search_hmmsearch_frame->Button(
	-text => "Select query",
	-width => 21,
	-command => sub { $query_path = $mw->getOpenFile();}
	)->grid(-row => 4, -column => 1, -padx => "1m", -pady => "1m");
	
	$search_hmmsearch_frame->Entry(
	-width => 40,
	-textvariable => \$query_path
	)->grid(-row => 4, -column => 2, -padx => "1m", -pady => "1m");
	
	
	
	$search_hmmsearch_frame->Label(-text => "Genomes Fasta")->grid(-row => 5, -column => 0, -padx => "1m", -pady => "1m",-sticky => 'w');
	$search_hmmsearch_frame->Button(
	-text => "Select directory",
	-width => 21,
	-command => sub { $fasta_dir = $mw->chooseDirectory}
	)->grid(-row => 5, -column => 1, -padx => "1m", -pady => "1m");
	
	$search_hmmsearch_frame->Entry(
	-width => 40,
	-textvariable => \$fasta_dir
	)->grid(-row => 5, -column => 2, -padx => "1m", -pady => "1m");
	
	

	$search_hmmsearch_frame->Label(-text => "Genomes GFF3")->grid(-row => 6, -column => 0, -padx => "1m", -pady => "1m", -sticky => 'w');
	$search_hmmsearch_frame->Button(
	-text => "Select directory",
	-width => 21,
	-command => sub { $gff_dir = $mw->chooseDirectory}
	)->grid(-row => 6, -column => 1, -padx => "1m", -pady => "1m");
	
	$search_hmmsearch_frame->Entry(
	-width => 40,
	-textvariable => \$gff_dir
	)->grid(-row => 6, -column => 2, -padx => "1m", -pady => "1m");
	
	
	
	$search_hmmsearch_frame->Label(-text => "CPU Cores")->grid(-row => 7, -column => 0,-padx => "1m", -pady => "1m", -sticky => 'w');
	$search_hmmsearch_frame->Entry(
	-width => 20,
	-textvariable => \$cores
	)->grid(-row => 7, -column => 2, -padx => "2m", -pady => "1m",-sticky => 'w');
	
	
	$search_hmmsearch_frame->Label(-text => "Options")->grid(-row => 8, -column => 0,-padx => "1m", -pady => "1m", -sticky => 'w');
	my @hmm_search = (0,0,0,0); #[1] Also search genomes already in db [2] Discard HmmReport Files afterwards
	$search_hmmsearch_frame->Checkbutton(-text => 'Ignore previous results', -onvalue => '1', -command => sub{
	&selectColumnsquery(0,$Tk::widget->cget('-variable'),\@hmm_search);
	}
	)->grid(-row => 8, -column => 1, -padx => "1m", -pady => "1m", -sticky => 'w');
	$search_hmmsearch_frame->Checkbutton(-text => 'Discard HMM reports', -onvalue => '1', -command => sub{
	&selectColumnsquery(1,$Tk::widget->cget('-variable'),\@hmm_search);
	}
	)->grid(-row => 8, -column => 2, -padx => "1m", -pady => "1m", -sticky => 'w');
	$search_hmmsearch_frame->Checkbutton(-text => 'Prodigal ORF prediction', -onvalue => '1', -command => sub{
	&selectColumnsquery(2,$Tk::widget->cget('-variable'),\@hmm_search);
	}
	)->grid(-row => 9, -column => 1, -padx => "1m", -pady => "1m", -sticky => 'w');
	$search_hmmsearch_frame->Checkbutton(-text => 'Read existing reports only', -onvalue => '1', -command => sub{
	&selectColumnsquery(3,$Tk::widget->cget('-variable'),\@hmm_search);
	}
	)->grid(-row => 9, -column => 2, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	#SEARCH BUTTON

$search_hmmsearch_frame->Label(-text => "Search may take some minutes. Progress is displayed in the command window")->grid(-row => 10,  -column => 0,-columnspan => 3, -pady => "1m", -padx => "1m", -sticky => 'w');

	$search_hmmsearch_frame->Button(
	-text => "Start Search",
	-width => 64,
	-command => sub{
			my $hash_reference = {
			algorithm => "$search_algorithm",
			cutoff_type => "$cutoff_type",
			cutoff => "$cutoff_value",
			thrs_file => "$thrs_path",
			query_file => "$query_path",
			fasta_dir => "$fasta_dir",
			gff_dir => "$gff_dir",
			cores => "$cores",
			redosearch => $hmm_search[0], #Include hmmsearch for all files eventhough they are already in DB
			cleanreports => $hmm_search[1], #discard the hmmreports afterwards
			prodigal => $hmm_search[2], #translate fna files
			redoreport => $hmm_search[3], #only take the report information from previous searches
			spec_thrs => 0,
			spec_hmm => 0,
			};
			&search($hash_reference);
			
			}
	)->grid(-row => 11, -column => 0, -columnspan => 3, -padx => "1m", -pady => "2m", -sticky => 'we');

#$search_hmmsearch_frame->pack(-side => 'top', -fill => 'y', -expand => 1, -anchor => 'nw');
$search_hmmsearch_frame->grid(-row => 1, -column => 0, -sticky => 'we');




##### Search Taxonomy frame
#
#Works only with NCBI taxdump. Collect the NCBITaxon identifier and matches the phylogeny of NCBI Taxon database to it
	
my $search_taxonomy_frame = $search_frame->Frame(-relief => 'ridge', -borderwidth => 2);	
$search_taxonomy_frame->Label(-text => "Collect assembly statistics and taxonomy information")->grid(-row => 2, -column => 0, -columnspan => 3, -padx => "1m", -padx => "1m", -pady => "1m", -sticky => 'w');
	
	$search_taxonomy_frame->Label(-text => "Assembly statistics")->grid(-row => 3, -column => 0, -padx => "1m", -pady => "1m", -sticky => 'w');
	$search_taxonomy_frame->Button(
	-text => "Select directory",
	-width => 21,
	-command => sub { $ass_dir = $mw->chooseDirectory}
	)->grid(-row => 3, -column => 1, -padx => "1m", -pady => "1m");
	
	$search_taxonomy_frame->Entry(
	-width => 40,
	-textvariable => \$ass_dir
	)->grid(-row => 3, -column => 2, -padx => "1m", -pady => "1m");
	
$search_taxonomy_frame->Label(-text => "")->grid(-row => 4, -columnspan => 3, -sticky => 'w');
$search_taxonomy_frame->Label(-text => "Search may take some minutes. Progress is displayed in the command window")->grid(-row => 5, -column => 0, -columnspan => 3, -padx => "1m", -pady => "1m", -sticky => 'w');

	$search_taxonomy_frame->Button(
	-text => "Get taxa",
	-width => 64,
	-command => sub{
			&assemblystats($ass_dir)} #Command subroutine from button for assembly search
	)->grid(-row => 6, -column => 0, -columnspan => 3, -padx => "1m", -pady => "1m", -sticky => 'we');

#$search_taxonomy_frame->pack(-side => 'top', -fill => 'y', -expand => 1, -anchor => 'nw');
$search_taxonomy_frame->grid(-row => 2, -column => 0, -sticky => 'we');









##### Search reanalyse frame
##
##
## Go through all reports for a single HMM profile with adjusted threshold
##
##

my $search_reanalyse_frame = $search_frame->Frame(-relief => 'ridge', -borderwidth => 2);	


$search_reanalyse_frame->Label(-text => "Recollect hits from previous reports for specified model with modified score threshold")->grid(-row => 3, -column => 0, -columnspan => 3, -padx => "1m", -pady => "1m", -sticky => 'w');

$search_reanalyse_frame->Label(-text => "Reanalyse model")->grid(-row => 4, -column => 0, -padx => "1m", -pady => "1m", -sticky => 'w');

my $browse_widget;
my $browse_options;
my $browse_hmm = "";
$browse_widget = $search_reanalyse_frame->BrowseEntry(
 	-label => '',
	-variable => \$browse_hmm,
 	-listcmd => sub { $browse_options = &getHMMList($query_path);
 			$browse_widget->delete(0,'end');
 			foreach(@$browse_options){
 			$browse_widget->insert('end',$_)
 			}
 			
 		},
 #-choices => $browse_options,
 	)->grid(-row => 4, -column => 1, -padx => "0m", -pady => "1m", -sticky => 'we');
	
	$search_reanalyse_frame->Label(-text => "Threshold score")->grid(-row => 5, -column => 0, -padx => "1m", -pady => "1m", -sticky => 'w');

my $threshold_score = 10;	
	$search_reanalyse_frame->Entry(
	-width => 22,
	-textvariable => \$threshold_score
	)->grid(-row => 5, -column => 1, -padx => "1m", -pady => "1m", -sticky => 'w');

	$search_reanalyse_frame->Label(-text => "Select the directory where the HMM reports are saved ")->grid(-row => 6, -column => 0,-columnspan => 3, -padx => "1m", -pady => "1m", -sticky => 'w');
	
my $hmmreport_path = $fasta_dir;
	$search_reanalyse_frame->Label(-text => "Report files")->grid(-row => 7, -column => 0, -padx => "1m", -pady => "1m", -sticky => 'w');
	$thrsButton = $search_reanalyse_frame->Button(
	-text => "Select directory",
	-width => 21,
	-command => sub { $hmmreport_path = $mw->chooseDirectory;}
	)->grid(-row => 7, -column => 1, -padx => "1m", -pady => "1m");

	$FileEntry = $search_reanalyse_frame->Entry(
	-width => 40,
	-textvariable => \$hmmreport_path
	)->grid(-row => 7, -column => 2, -padx => "1m", -pady => "1m");

	$search_reanalyse_frame->Label(-text => "Gather all hits for specified hmm above threshold score")->grid(-row => 8, -column => 0, -columnspan => 3,-pady => "1m", -sticky => 'w');
	$search_reanalyse_frame->Button(
	-text => "Reanalyse",
	-width => 64,
	-command => sub{
	print "Variable ",$browse_hmm," Score ",$threshold_score,"\n";
	my $hash_reference = {
			algorithm => "$search_algorithm",
			cutoff_type => "$cutoff_type",
			cutoff => "$cutoff_value",
			thrs_file => "$thrs_path",
			query_file => "$query_path",
			fasta_dir => "$hmmreport_path",
			gff_dir => "$gff_dir",
			cores => "$cores",
			redosearch => $hmm_search[0], #Include hmmsearch for all files eventhough they are already in DB
			cleanreports => $hmm_search[1], #discard the hmmreports afterwards
			prodigal => $hmm_search[2], #translate fna files
			redoreport => $hmm_search[3], #only take the report information from previous searches
			spec_thrs => $threshold_score,
			spec_hmm => $browse_hmm,
			
			};
			
			&search($hash_reference);
			
			}
	)->grid(-row => 8, -column => 0, -columnspan => 3, -padx => "0m", -pady => "1m", -sticky => 'we');
	
#$search_reanalyse_frame->pack(-side => 'top', -fill => 'y', -expand => 1, -anchor => 'nw');
$search_reanalyse_frame->grid(-row => 3, -column => 0, -sticky => 'we');

####FIND SYNTHENIC GENES

#Input specific cluster or hmm
#number of genes in vicinity up and downstream

my $cluster_bot_frame = $search_frame->Frame(-relief => 'ridge', -borderwidth => 2);

$cluster_bot_frame->Label(-text => "Collect syntenic genes upstream/downstream of specific protein or genecluster")->grid(-row => 11, -columnspan => 3, -sticky => 'w');
#$cluster_bot_frame->Label(-text => "File containing the synteny patterns")->grid(-row => 12, -columnspan => 2, -sticky => 'w');
	my $var1;
	my $query_options = [[genecluster => "keyword"],[protein => "HMM"]];
	$cluster_bot_frame->Optionmenu(
	 -options => $query_options,
	 -variable => \$var1,				   
	 -width => 20,
	 -command => sub {
	 		  }
		)->grid(-row => 13, -column => 0, -padx => "1m", -pady => "1m");
		
	my $queryoption1;
	$cluster_bot_frame->Optionmenu(
	 -options => [[is => " = "],["like" => " rlike "]],
	 -variable => \$queryoption1,				   
	# -textvariable => \$tkeyword_group2,		   
	 -width => 10,
	 -command => sub {
	 		  }
		)->grid(-row => 13, -column => 1, -padx => "1m", -pady => "1m");
		
		
	my $query_widget1;
	my $query_variable1 = "-";
	$query_widget1 = $cluster_bot_frame->BrowseEntry(
	 -label => '',
	 
	 -variable => \$query_variable1,
	 -listcmd => sub {	my $browse_options = &getKeyList($var1);
				$query_widget1->delete(0,'end');
				$query_widget1->insert('end',"-");
	 			foreach(@$browse_options){
	 			$query_widget1->insert('end',$_);
	 			}
	 			
	 		},
	 #-choices => $browse_options,
	 	)->grid(-row => 13, -column => 2, -padx => "1m", -pady => "1m");


	
	my $genNum = 5;
	$cluster_bot_frame->Label(-text => "Number of genes upstream/downstream")->grid(-row => 16, -column => 0,-padx => "1m", -pady => "1m", -sticky => 'w');
	$cluster_bot_frame->Entry(
	-width => 10,
	-textvariable => \$genNum
	)->grid(-row => 16, -column => 1, -padx => "1m", -pady => "1m", -sticky => 'w');



$cluster_bot_frame->Button(
	-text => "Search syntenic genes",
	-width => 79,
	-command => 
	sub{
	 	
	 	&findsynteny($var1,$queryoption1,$query_widget1,$genNum);

	 }
	)->grid(-row => 17, -column => 0, -columnspan => 3, -padx => "0m", -pady => "1m", -sticky => 'we');
	
$cluster_bot_frame->grid(-row => 17, -column => 0, -sticky => 'we');



}


########################################	
###	####Sub Menu 3 Frame loading####
#	#### Cluster Frame	     ####
#
#
#

sub clusterframe{
my $ini_hashref =shift;
my $distance = 3500;
if($ini_hashref->{Distance}){$distance = $ini_hashref->{Distance}}
my $pattern_path = "$cwd/Patterns.txt";
if($ini_hashref->{Patterns_File}){$pattern_path = $ini_hashref->{Patterns_File}}
	
	
	
my $cluster_synteny_frame = $cluster_frame->Frame(-relief => 'ridge', -borderwidth => 2);	

	#Distance between genes
	$cluster_synteny_frame->Label(-text => "Allowed number of nucleotides between two genes to consider them as neighbouring")->grid(-row => 0, -columnspan => 2, -padx => "1m", -pady => "1m", -sticky => 'w');
	$cluster_synteny_frame->Label(-text => "Distance [nt]")->grid(-row => 1, -column => 0, -sticky => 'w');
	$cluster_synteny_frame->Entry(
	-width => 21,
	-textvariable => \$distance
	)->grid(-row => 1, -column => 1, -padx => "1m", -pady => "1m", -sticky => 'w');

$cluster_synteny_frame->Label(-text => "(Optional) Limit the organisms for clustering and synteny search")->grid(-row => 2, -columnspan => 3, -sticky => 'w');
$cluster_synteny_frame->Label(-text => "Classification:")->grid(-row => 3, -column => 0, -sticky => 'w');
$cluster_synteny_frame->Label(-text => "Group:")->grid(-row => 3, -column => 1, -sticky => 'w');


my ($searchcluster_group,$tsearchcluster_group);
my $searchclustergroup_options = [[All=> 0], [Phylum=>'Phylum'], [Class=>'Class'], [Order=>'Ordnung'], [Family=>'Family'],[Genus=>'Genus'],[Species=>'Species']];#Variables for the second Browse Entry menu

#Browse Entry for Group (right side) TODO disable if searchclustergroup is all => 0
my $browse_widget;
my $browse_options;
my $browse_variable = "0";
$browse_widget = $cluster_synteny_frame->BrowseEntry(
 -label => '',
 
 -variable => \$browse_variable,
 -listcmd => sub {	$browse_options = &getPhyloList($searchcluster_group);
 			$browse_widget->delete(0,'end');
 			foreach(@$browse_options){
 			$browse_widget->insert('end',$_)
 			}
 			
 		},
 #-choices => $browse_options,
 	)->grid(-row => 4, -column => 1, -padx => "0m", -pady => "1m", -sticky => 'we');
 	
 
 #Optionmenu with selection of the phylo group
$cluster_synteny_frame->Optionmenu(
 -options => $searchclustergroup_options,
 -variable => \$searchcluster_group,				   
 -textvariable => \$tsearchcluster_group,		   
 -width => 20,
 -command => sub {
 			$browse_widget->delete(0,"end");
 			$browse_variable = "0";
 
 		  }
	)->grid(-row => 4, -column => 0, -padx => "1m", -pady => "1m", -sticky => 'we');

	
	
	
	$cluster_synteny_frame->Label(-text => "File containing the synteny patterns")->grid(-row => 5, -columnspan => 2,  -padx => "1m", -pady => "1m", -sticky => 'w');
	
	$cluster_synteny_frame->Button(
	-text => "Select Pattern File",
	
	-command => sub { $pattern_path = $mw->getOpenFile();}
	)->grid(-row => 6, -column => 0, -padx => "1m", -pady => "0m", -sticky => 'we');
	
	$cluster_synteny_frame->Entry(
	-width => 40,
	-textvariable => \$pattern_path
	)->grid(-row => 6, -column => 1, -padx => "1m", -pady => "0m",-sticky => 'we');
	

	
#####SEARCH BUTTON	
$cluster_synteny_frame->Label(-text => "Operon prediction and keyword assignment may take some minutes.")->grid(-row => 7, -columnspan => 2,  -padx => "1m", -pady => "1m", -sticky => 'w');
$cluster_synteny_frame->Label(-text => "Progress is displayed in the command window")->grid(-row => 8, -columnspan => 2, -padx => "1m", -pady => "1m", -sticky => 'w');

	my $cluster_button_frame = $cluster_synteny_frame->Frame(-relief => 'flat', -borderwidth => 1);

		$cluster_button_frame->Button(
		-text => "Start Clustering",
		-width => 32,
		-command => sub{&cluster($distance,$searchcluster_group,$browse_variable,$pattern_path,"synteny")}
		#)->grid(-row => 1, -column => 0, -padx => "1m", -pady => "1m", -sticky => 'we');
		)->pack(-side => 'left', -fill => 'x', -expand => 1,-padx => "1m", -pady => "1m", -anchor => 'nw');
		$cluster_button_frame->Button(
		-text => "Start Keyword Assignment",
		-width => 32,
		-command => sub{&cluster($distance,$searchcluster_group,$browse_variable,$pattern_path,"keyword")}
		#)->grid(-row => 1, -column => 1, -padx => "1m", -pady => "1m", -sticky => 'we');
		)->pack(-side => 'right', -fill => 'x', -expand => 1, -padx => "1m", -pady => "1m", -anchor => 'nw' );
	$cluster_button_frame->grid(-row => 9, -column => 0, -columnspan => 2, -padx => "0m", -pady => "0m", -sticky => 'we');
	
	
$cluster_synteny_frame->grid(-row => 1, -column => 0, -columnspan => 2, -padx => "0m", -pady => "0m", -sticky => 'we');
	



}



########################################	
###	####Sub Menu 4 Frame loading####
#	#### Output Search Frame    ####
#
#
#
sub outputqueryframe{

##### Sub Frame 1
my $mysql_head_frame = $mysql_frame->Frame(-relief => 'flat', -width => 1000, -borderwidth => 3);
my @mysql_search = ("","","","","","","","",""); #proteinID, locustag,genomeID, HMM, keyword, start,end, strand, Species


	$mysql_head_frame->Label(-text => "Select columns to be displayed:")->grid(-row => 0, -columnspan => 3, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	
	$mysql_head_frame->Checkbutton(-text => 'Protein ID', -onvalue => 'Proteins.protein_ID', -command => sub{
	&selectColumnsquery(0,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 1, -column => 0, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	
	$mysql_head_frame->Checkbutton(-text => 'Locustag', -onvalue => 'Proteins.locustag', -command => sub{
	&selectColumnsquery(1,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 1, -column => 1, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	
	$mysql_head_frame->Checkbutton(-text => 'Genome ID', -onvalue => 'Proteins.genome_ID', -command => sub{
	&selectColumnsquery(2,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 1, -column => 2, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	
	$mysql_head_frame->Checkbutton(-text => 'Species', -onvalue => 'Genomes.Species', -command => sub{
	&selectColumnsquery(3,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 1, -column => 3, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	
	$mysql_head_frame->Checkbutton(-text => 'Protein type', -onvalue => 'Proteins.HMM', -command => sub{
	&selectColumnsquery(4,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 1, -column => 4, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	
	
	
	$mysql_head_frame->Checkbutton(-text => 'Cluster type', -onvalue => 'keyword', -command => sub{
	&selectColumnsquery(5,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 2, -column => 0, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	
		$mysql_head_frame->Checkbutton(-text => 'Contig', -onvalue => 'Proteins.Contig', -command => sub{
	&selectColumnsquery(6,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 2, -column => 1, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	
	$mysql_head_frame->Checkbutton(-text => 'Start', -onvalue => 'Proteins.start', -command => sub{
	&selectColumnsquery(7,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 2, -column => 2, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	
	$mysql_head_frame->Checkbutton(-text => 'End', -onvalue => 'Proteins.end', -command => sub{
	&selectColumnsquery(8,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 2, -column => 3, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	
	$mysql_head_frame->Checkbutton(-text => 'Strand', -onvalue => 'Proteins.strand', -command => sub{
	&selectColumnsquery(9,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 2, -column => 4, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	
	### Phylogy selections
	$mysql_head_frame->Checkbutton(-text => 'Phylum', -onvalue => 'Genomes.Phylum', -command => sub{
	&selectColumnsquery(10,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 3, -column => 0, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	$mysql_head_frame->Checkbutton(-text => 'Class', -onvalue => 'Genomes.Class', -command => sub{
	&selectColumnsquery(11,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 3, -column => 1, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	$mysql_head_frame->Checkbutton(-text => 'Order', -onvalue => 'Genomes.Ordnung', -command => sub{
	&selectColumnsquery(12,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 3, -column => 2, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	$mysql_head_frame->Checkbutton(-text => 'Family', -onvalue => 'Genomes.Family', -command => sub{
	&selectColumnsquery(13,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 3, -column => 3, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	$mysql_head_frame->Checkbutton(-text => 'Genus', -onvalue => 'Genomes.Genus', -command => sub{
	&selectColumnsquery(14,$Tk::widget->cget('-variable'),\@mysql_search);
	}
	)->grid(-row => 3, -column => 4, -padx => "1m", -pady => "1m", -sticky => 'w');
	
	
	
	
	
	
	
	
	
$mysql_head_frame->pack(-side => 'top', -fill => 'y', -expand => 1, -anchor => 'nw');	
	

##### Sub Frame 2
my $mysql_select_frame = $mysql_frame->Frame(-relief => 'flat', -width => 1000, -borderwidth => 3);
$mysql_select_frame->Label(-text => "")->grid(-row => 0, -column => 0, -sticky => 'w');
	
$mysql_select_frame->Label(-text => "(Optional) Select the genomes of the organisms to be searched:")->grid(-row => 1, -columnspan => 3, -column => 0, -sticky => 'w');

$mysql_select_frame->Label(-text => "Classification:")->grid(-row => 3, -column => 0, -padx => "1m", -sticky => 'w');
$mysql_select_frame->Label(-text => "Group:")->grid(-row => 3, -column => 2,-padx => "1m", -sticky => 'w');



my ($searchcluster_group,$tsearchcluster_group);
my $searchclustergroup_options = [[All=>'All'], [Phylum=>'Phylum'], [Class=>'Class'], [Order=>'Ordnung'], [Family=>'Family'], [Genus=>'Genus'],[Species=>'Species']];
#Variables for the second Browse Entry menu
my $browse_widget;
my $browse_options = "";
my $browse_variable = "All";
$browse_widget = $mysql_select_frame->BrowseEntry(
 -label => '',
 
 -variable => \$browse_variable,
 -listcmd => sub {	$browse_options = &getPhyloList($searchcluster_group);
 			$browse_widget->delete(0,'end');
 			foreach(@$browse_options){
 			$browse_widget->insert('end',$_)
 			}
 			
 		},
 #-choices => $browse_options,
 	)->grid(-row => 4, -column => 2, -padx => "1m", -pady => "1m");
 	
 	
$mysql_select_frame->Optionmenu(
 -options => $searchclustergroup_options,
 -variable => \$searchcluster_group,				   
 -textvariable => \$tsearchcluster_group,		   
 -width => 20,
 -command => sub {
 			$browse_widget->delete(0,"end");
 			$browse_variable = "All";
 
 		  }
	)->grid(-row => 4, -column => 0, -padx => "1m", -pady => "1m");
	

	
$mysql_select_frame->Label(-text => "(Optional) Limit the search to genomes with a specific pattern:")->grid(-row => 5, -columnspan => 3, -sticky => 'w');

my ($keyword_group2,$tkeyword_group2);
my $keywordgroup_options2 = [[genecluster => 'keyword']];
#Variables for the second Browse Entry menu
my $browse_widget2;
my $browse_options2;
my $browse_variable2 = "All";
$browse_widget2 = $mysql_select_frame->BrowseEntry(
 -label => '',
 -width => 20,
 -variable => \$browse_variable2,
 -listcmd => sub {	my $browser = &getKeyList($keyword_group2); 
			$browse_widget2->delete(0,'end');
 			foreach(@$browser){
 			$browse_widget2->insert('end',$_)
 			}
 			
 		},
 #-choices => $browse_options,
 	)->grid(-row => 6, -column => 2, -padx => "1m", -pady => "1m");


	my $keywordoptioncon;
	$mysql_select_frame->Optionmenu(
	 -options => [[is => " = "],["not" => " != "],["like" => " rlike "]],
	 -variable => \$keywordoptioncon,				   
	 -width => 10,
	 -command => sub {
	 		  }
		)->grid(-row => 6, -column => 1, -padx => "1m", -pady => "1m");
			
 	
$mysql_select_frame->Optionmenu(
 -options => $keywordgroup_options2,
 -variable => \$keyword_group2,				   
 -textvariable => \$tkeyword_group2,
 -width => 20,
 -command => sub {
 			$browse_widget2->delete(0,"end");
 			$browse_variable2 = "All";
 
 		  }
	)->grid(-row => 6, -column => 0, -padx => "1m", -pady => "1m");
	
$mysql_select_frame->Label(-text => "")->grid(-row => 7, -column => 0, -sticky => 'w');
	


$mysql_select_frame->pack(-side => 'top', -fill => 'y', -expand => 1, -anchor => 'nw');	




##### Sub Frame 3
my $mysql_query_frame = $mysql_frame->Frame(-relief => 'flat', -width => 1000, -borderwidth => 3);

$mysql_query_frame->Label(-text => "Specify the search the proteins where:")->grid(-row => 0, -columnspan => 2, -column => 0, -sticky => 'w');

	my $var1;
	my $query_options = [[genecluster => "keyword"],[protein => "HMM"]];
	$mysql_query_frame->Optionmenu(
	 -options => $query_options,
	 -variable => \$var1,				   
	 -width => 20,
	 -command => sub {
	 		  }
		)->grid(-row => 1, -column => 0, -padx => "1m", -pady => "1m");
		
	my $queryoption1;
	$mysql_query_frame->Optionmenu(
	 -options => [[is => " = "],["not" => " != "],["like" => " rlike "]],
	 -variable => \$queryoption1,				   
	# -textvariable => \$tkeyword_group2,		   
	 -width => 10,
	 -command => sub {
	 		  }
		)->grid(-row => 1, -column => 1, -padx => "1m", -pady => "1m");
		
		
	my $query_widget1;
	my $query_variable1 = "-";
	$query_widget1 = $mysql_query_frame->BrowseEntry(
	 -label => '',
	 
	 -variable => \$query_variable1,
	 -listcmd => sub {	my $browse_options = &getKeyList($var1);
				$query_widget1->delete(0,'end');
				$query_widget1->insert('end',"-");
	 			foreach(@$browse_options){
	 			$query_widget1->insert('end',$_);
	 			}
	 			
	 		},
	 #-choices => $browse_options,
	 	)->grid(-row => 1, -column => 2, -padx => "1m", -pady => "1m");
	my $connector1 = "-";	
	$mysql_query_frame->Optionmenu(
	 -options => [["-" => "-"],["and" => " AND "],["or" => " OR "]],
	 -variable => \$connector1,
	# -textvariable => \$tkeyword_group2,		   
	 -width => 10,
	 -command => sub {
	 		  }
		)->grid(-row => 1, -column => 3, -padx => "1m", -pady => "1m");
		
		
		
	my $var2;
	$mysql_query_frame->Optionmenu(
	 -options => $query_options,
	 -variable => \$var2,				   
	 -width => 20,
	 -command => sub {
	 		  }
		)->grid(-row => 2, -column => 0, -padx => "1m", -pady => "1m");
		
	my $queryoption2;
	$mysql_query_frame->Optionmenu(
	 -options => [[is => " = "],["not" => " != "],["like" => " rlike "]],
	 -variable => \$queryoption2,				   
	 -width => 10,
	 -command => sub {
	 		  }
		)->grid(-row => 2, -column => 1, -padx => "1m", -pady => "1m");
		
		
	my $query_widget2;
	my $query_variable2 = "-";
	$query_widget2 = $mysql_query_frame->BrowseEntry(
	 -label => '',
	 
	 -variable => \$query_variable2,
	 -listcmd => sub {	my $browse_options = &getKeyList($var2);
				$query_widget2->delete(0,"end");
				$query_widget2->insert('end',"-");
	 			foreach(@$browse_options){
	 			$query_widget2->insert('end',$_);
	 			}
	 			
	 		},
	 #-choices => $browse_options,
	 	)->grid(-row => 2, -column => 2, -padx => "1m", -pady => "1m");
	
	my $connector2 = "-";
	$mysql_query_frame->Optionmenu(
	 -options => [["-" => "-"],["and" => " AND "],["or" => " OR "]],
	 -variable => \$connector2,				   
	# -textvariable => \$tkeyword_group2,		   
	 -width => 10,
	 -command => sub {
	 		  }
		)->grid(-row => 2, -column => 3, -padx => "1m", -pady => "1m");
		
		
		
 	my $var3;
	$mysql_query_frame->Optionmenu(
	 -options => $query_options,
	 -variable => \$var3,				   
	 -width => 20,
	 -command => sub {
	 		  }
		)->grid(-row => 3, -column => 0, -padx => "1m", -pady => "1m");
		
	my $queryoption3;
	$mysql_query_frame->Optionmenu(
	 -options => [[is => " = "],["not" => " != "],["like" => " rlike "]],
	 -variable => \$queryoption3,
	 -width => 10,
	 -command => sub {
	 		  }
		)->grid(-row => 3, -column => 1, -padx => "1m", -pady => "1m");
		
		
	my $query_widget3;
	my $query_variable3 = "-";
	$query_widget3 = $mysql_query_frame->BrowseEntry(
	 -label => '',
	 
	 -variable => \$query_variable3,
	 -listcmd => sub {	my $browse_options = &getKeyList($var3);
				$query_widget3->delete(0,"end");
				$query_widget3->insert('end',"-");
	 			foreach(@$browse_options){
	 			$query_widget3->insert('end',$_);
	 			}
	 			
	 		},
	 #-choices => $browse_options,
	 	)->grid(-row => 3, -column => 2, -padx => "1m", -pady => "1m");
	
	my $connector3 = "-";	
	 $mysql_query_frame->Optionmenu(
	 -options => [["-" => "-"],["and" => " AND "],["or" => " OR "]],
	 -variable => \$connector3,
	# -textvariable => \$tkeyword_group2,		   
	 -width => 10,
	 -command => sub {
	 		  }
		)->grid(-row => 3, -column => 3, -padx => "1m", -pady => "1m");
		
			
	 my $var4;
	$mysql_query_frame->Optionmenu(
	 -options => $query_options,
	 -variable => \$var4,				   
	 -width => 20,
	 -command => sub {
	 		  }
		)->grid(-row => 4, -column => 0, -padx => "1m", -pady => "1m");
		
	my $queryoption4;
	$mysql_query_frame->Optionmenu(
	 -options => [[is => " = "],["not" => " != "],["like" => " rlike "]],
	 -variable => \$queryoption4,
	 -width => 10,
	 -command => sub {
	 		  }
		)->grid(-row => 4, -column => 1, -padx => "1m", -pady => "1m");
		
		
	my $query_widget4;
	my $query_variable4 = "-";
	$query_widget4 = $mysql_query_frame->BrowseEntry(
	 -label => '',
	 
	 -variable => \$query_variable4,
	 -listcmd => sub {	my $browse_options = &getKeyList($var4);
				$query_widget4->delete(0,"end");
				$query_widget4->insert('end',"-");
	 			foreach(@$browse_options){
	 			$query_widget4->insert('end',$_);
	 			}
	 			
	 		},
	 #-choices => $browse_options,
	 	)->grid(-row => 4, -column => 2, -padx => "1m", -pady => "1m");
	 	
	my $connector4 = "-";	
	 $mysql_query_frame->Optionmenu(
	 -options => [["-" => "-"],["and" => " AND "],["or" => " OR "]],
	 -variable => \$connector4,				   
	# -textvariable => \$tkeyword_group2,		   
	 -width => 10,
	 -command => sub {
	 		  }
		)->grid(-row => 4, -column => 3, -padx => "1m", -pady => "1m");
		
			
	 	
$mysql_query_frame->pack(-side => 'top', -fill => 'y', -expand => 1, -anchor => 'nw');	

##### Sub Frame
###Search Button for gui supported search
my $mysqlquery = "";

my $GUIResultButton = $mysql_frame->Button(
	-text => 'Retrieve filtered results',
	-command => 
	sub{
	 	my $mysql_search_ref = \@mysql_search;
	 	my $w1 = [$var1,$var2,$var3,$var4];
	 	my $w2 = [$queryoption1, $queryoption2, $queryoption3, $queryoption4];
	 	my $w3 = [$query_variable1, $query_variable2, $query_variable3, $query_variable4];
	 	my $con = [$connector1, $connector2, $connector3, $connector4];

	 	$mysqlquery = &generatequery($mysql_search_ref,$w1,$w2,$w3,$con,$browse_variable,$searchcluster_group,$keyword_group2,$keywordoptioncon,$browse_variable2);
	 	#$textquery->Contents("$mysqlquery\n");
		my $heading =	&getHeadings($mysqlquery);
				&resultframe($mysqlquery,$heading);

	 }
	)->pack(-side => 'top', -fill => 'x', -expand => 0, -anchor => 'nw');
	



###Panel for custom query search

my $mysql_textquery_frame = $mysql_frame->Frame(-relief => 'flat', -width => 1000, -borderwidth => 5);

$mysql_textquery_frame->Label(-text => "(Alternative)")->pack(-side => 'top',-fill => 'y', -expand => 0, -anchor => 'nw');
$mysql_textquery_frame->Label(-text => "Altering search command or direct entering mysql injection for complex database requests:")->pack(-side => 'top',-fill => 'y', -expand => 0, -anchor => 'nw');
my $textquery;
my $generatetexttext = "Generate query";
$mysql_textquery_frame->Button(
	-textvariable => \$generatetexttext,
	-width => 15,
	-command => 
	sub{
	 	my $mysql_search_ref = \@mysql_search;
	 	my $w1 = [$var1,$var2,$var3,$var4];
	 	my $w2 = [$queryoption1, $queryoption2, $queryoption3, $queryoption4];
	 	my $w3 = [$query_variable1, $query_variable2, $query_variable3, $query_variable4];
	 	my $con = [$connector1, $connector2, $connector3, $connector4];

	 	$mysqlquery = &generatequery($mysql_search_ref,$w1,$w2,$w3,$con,$browse_variable,$searchcluster_group,$keyword_group2,$keywordoptioncon,$browse_variable2);
	 	$textquery->Contents("$mysqlquery\n");
	 	
	 }
	)->pack(-side => 'top', -expand => 0, -anchor => 'nw');

$textquery = $mysql_textquery_frame->Text( -height => 8, -width => 80)->pack(-side => 'bottom', -fill => 'both', -expand => 0, -anchor => 'sw');

$mysql_textquery_frame->pack(-side => 'top', -fill => 'none', -expand => 1, -anchor => 'nw');	


my $QueryResultButton = $mysql_textquery_frame->Button(

	-text => 'Retrieve filtered results by query',
	-command => 
	sub{
		my $query = $textquery->get('1.0','end');
		#print "My query:$query";
		
		unless(length $query < 5){
		my $heading =	&getHeadings($query);
				&resultframe($query,$heading);
		}else{
		my $mysql_search_ref = \@mysql_search;
	 	my $w1 = [$var1,$var2,$var3,$var4];
	 	my $w2 = [$queryoption1, $queryoption2, $queryoption3, $queryoption4];
	 	my $w3 = [$query_variable1, $query_variable2, $query_variable3, $query_variable4];
	 	my $con = [$connector1, $connector2, $connector3, $connector4];

	 	$mysqlquery = &generatequery($mysql_search_ref,$w1,$w2,$w3,$con,$browse_variable,$searchcluster_group,$keyword_group2,$keywordoptioncon,$browse_variable2);
	 	$textquery->Contents("$mysqlquery\n");
		my $heading =	&getHeadings($mysqlquery);
				&resultframe($mysqlquery,$heading);
		}
	 }
	)->pack(-side => 'bottom', -fill => 'x', -expand => 1, -anchor => 'n');
	
my $ResetButton;
$ResetButton = $mysql_frame->Button(

	-text => "Reset",
	-command => 
	sub{
		&clearWorkframe();
		$mysql_head_frame->destroy;
		$mysql_select_frame->destroy;
		$mysql_query_frame->destroy;
		$mysql_textquery_frame->destroy;
		$GUIResultButton->destroy;
		#$QueryResultButton->destroy;
		$ResetButton->destroy;
		&outputqueryframe();
		$mysql_frame->pack(-side => 'left', -fill => 'both', -expand => 1, -anchor => 'nw');
	 }
	)->pack(-side => 'bottom', -fill => 'x', -expand => 1, -padx => "1m", -pady => "3m", -anchor => 'sw');
}


###########################################################################################

sub resultframe{

my $query = shift;
my $heading_ref = shift;
&checkDBIconnection();
my $sql = $dbh->prepare($query);
$sql->execute();

$result_frame->destroy();
$result_frame = $work_frame->Frame(-relief => 'ridge', -width => 1000, -borderwidth => 3);
###
###Sub Frame for the display results of the query
###



my $result_text_frame = $result_frame->Frame();

$result_frame->Label( -text => "Results:")->pack(-side => 'top', -fill => 'none', -expand => 0, -anchor => 'w', );

my $textfield = $result_text_frame->ROText(-wrap => 'none')->grid(-row => 0, -column => 0, -sticky => 'nswe');
	
my $scrollbar_vertical = $result_text_frame->Scrollbar(
	-orient => 'v',
	-command => [yview => $textfield]
)->grid(-row => 0, -column => 1, -sticky => 'ns');
my $scrollbar_horizontal = $result_text_frame->Scrollbar(
	-orient => 'h',
	-command => [xview => $textfield]
)->grid(-row => 1, -column => 0, -sticky => 'we');

$textfield->configure(
	-yscrollcommand => ['set', $scrollbar_vertical],
	-xscrollcommand => ['set', $scrollbar_horizontal],
);

$result_text_frame->gridRowconfigure(0, -weight  => 1);
$result_text_frame->gridColumnconfigure(0, -weight  => 1);
$result_text_frame->pack(-side => 'top', -fill => 'both', -expand=>1, -anchor => 'nw');



###
###Sub Frame for the saving results of the query
###

my $result_save_frame = $result_frame->Frame();
my $save_path = "$cwd/Save";
my $save_dir = "$cwd/Save";
	$result_save_frame->Label(-text => "Directory:")->grid(-row => 2, -column => 0, -padx => "1m", -pady => "3m");
	
	
	$result_save_frame->Button(
	-text => "Save as fasta",
	-width => 11,
	-command => sub { 
			$save_path = $mw->getSaveFile(
			 -initialdir => $save_dir,
			 -filetypes => [["Fasta", [qw/.fasta .faa .fas/]], ["All files", ['*']]]);

			if(defined $save_path && $save_path ne ""){
			&SaveFasta($save_path,$query,$heading_ref);
			}
			}
	
	)->grid(-row => 2, -column => 1, -padx => "1m", -pady => "1m");
	
	#$result_save_frame->Entry(
	#-width => 50,
	#-textvariable => \$save_fasta_path
	#)->grid(-row => 2, -column => 2, -padx => "1m", -pady => "1m");
	
	$result_save_frame->Button(
	-text => "Save as csv",
	-width => 11,
	-command => sub { 
	
			$save_path = $mw->getSaveFile(
			 -initialdir => $save_dir,
			 -filetypes => [["Csv", [qw/.csv/]], ["All files", ['*']]]);

			if(defined $save_path && $save_path ne ""){
			&SaveCsv($save_path,$query,$heading_ref);
			}
	
			}
	)->grid(-row => 2, -column => 2, -padx => "1m", -pady => "1m");
	
	
$result_save_frame->pack(-side => 'top', -fill => 'x', -expand=> 0, -anchor => 'nw');


##Fill the textfield with content
##TK table funktioniert hier nicht weil viel zu groß und ressourcen raubend

$textfield->insert('end',join(" \t| ",@$heading_ref)."\n");
while(my @line = $sql->fetchrow_array()){
@line = map{$_//"\t"} @line;
my $line = join(" \t| ",@line);
$line .= "\n";
$textfield->insert('end',$line);

}


}


###########################################################################################
################ SUBROUTINES ##############################################################
################ FOR GUI PROCESSING  ######################################################
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
	my @array = split("\t",$_);
	$Hashtable{$array[0]}=($array[1]);
	}
}
close $fh;

return \%Hashtable;

}

sub getPhyloList{
my $p = shift;
my @array;

unless($p eq "All" || $p eq 0 || $p eq "" || $p eq "-"){
my $query = "SELECT DISTINCT $p FROM Genomes ORDER BY $p;";
&checkDBIconnection();
$query = $dbh->prepare($query);
$query->execute();
	while(my @phylogeny = $query->fetchrow_array()){
	push(@array,$phylogeny[0]);
	}
}
return \@array;
}



sub getKeyList{
my $p = shift;
my @array;


unless($p eq "All" || $p eq 0 || $p eq "" || $p eq "-"){
my $query ;

if($p eq "keyword"){
$query = "SELECT DISTINCT $p FROM Functions ORDER BY $p;";
}elsif($p eq "PHY"){
$query = "SELECT DISTINCT $p FROM Proteins ORDER BY $p;";
}elsif($p eq "HMM"){ #TODO delete this
$query = "SELECT DISTINCT HMM FROM Proteins ORDER BY HMM;";
}
&checkDBIconnection();
$query = $dbh->prepare($query);
$query->execute();
while(my @phylogeny = $query->fetchrow_array()){
push(@array,$phylogeny[0]);

}
}
return \@array;
}




sub clearWorkframe{

$dbcon_frame->packForget();
$search_frame->packForget();
$cluster_frame->packForget();
$mysql_frame->packForget();
$result_frame->packForget();
}



###########################################################################################
################ SUBROUTINES ##############################################################
################ FOR DATA PROCESSING ######################################################
###########################################################################################
###
###
###
###
###
sub search{
my $hash_ref = shift;

my $COMMAND = "perl bin/Search.pl -u $db_user -p $db_pass -b $db_name -a ".$hash_ref->{algorithm}." -c ".$hash_ref->{cutoff_type}." -v ".$hash_ref->{cutoff}." -t ".$hash_ref->{thrs_file}." -q ".$hash_ref->{query_file}." -f ".$hash_ref->{fasta_dir}." -g ".$hash_ref->{gff_dir}." -n ".$hash_ref->{cores}." -r ".$hash_ref->{redosearch}." -l ".$hash_ref->{cleanreports}." -o ".$hash_ref->{prodigal}." -k ".$hash_ref->{redoreport}." -w ".$hash_ref->{spec_hmm}." -y ".$hash_ref->{spec_thrs};

print $COMMAND,"\n";
system($COMMAND);
}



sub cluster{
my $distance = shift;
my $Phylogroup = shift;
my $Specgroup = shift;
my $PatternFile = shift;
my $modus = shift;

my $command = "perl bin/Cluster.pl -d $distance -f $PatternFile -u $db_user -b $db_name -p $db_pass -q $Phylogroup -r $Specgroup -m $modus";

system($command);


}


sub assemblystats{
my $dir = shift;
my $command = "perl bin/Assemblystats.pl -u $db_user -b $db_name -p $db_pass -a $dir -c $cwd";
system($command);
}

sub selectColumnsquery{
my $index = shift;
my $status = shift;
my $array = shift;

if($$status){
$array->[$index] = $$status;
}else{
$array->[$index] = 0;
}

}


sub generatequery{
my $mysql_search = shift;
my $var = shift;
my $queryoption = shift;
my $query_variable = shift;
my $connector = shift;

my $browse_variable = shift;
my $searchcluster_group = shift;

my $browse_options2 = shift;
my $keyword_option_con = shift;
my $keyword_group2 = shift;

my $query = "SELECT DISTINCT ";
	foreach my $i (@$mysql_search){
	if(defined $i && $i ne ""){
	 $query .= $i.",";
	 }
	}
	chop($query);
	if($query eq "SELECT DISTINCT"){ $query .= " *";}

$query .= " FROM Proteins LEFT JOIN Genomes ON Proteins.genome_ID = Genomes.genome_ID LEFT JOIN Functions ON Proteins.cluster_ID = Functions.cluster_ID ";
my $where = "WHERE ";

	foreach my $i (0..3){

		if($query_variable->[$i] ne "-"){
		$query .= $where; $where = "";
		$query .= $var->[$i];
		$query .= $queryoption->[$i];
		$query .= "\"".$query_variable->[$i]."\"";
			if($connector->[$i] ne "-"){
			$query .= $connector->[$i];
			}else{
			
			last;
			}
		
		}
		
	}

	
	if($searchcluster_group ne "All"){
	$query .= " AND $searchcluster_group = \"$browse_variable\" ";
	}
	if($keyword_group2 ne "All"){
	
	$query .= "AND Proteins.genome_ID IN (SELECT genome_ID from Functions LEFT JOIN Clusters ON Functions.cluster_ID = Clusters.cluster_ID WHERE keyword $keyword_option_con \"$keyword_group2\") ";
	
	}
	
	$query .= "ORDER BY contig,start;";
	return $query;
}


sub SaveFasta{

my $save_fasta_path = shift;
my $query = shift;

my @result_table;
my $seqio_obj;

my ($name,$dir,$ext) = fileparse($save_fasta_path, qr/\.[^.]*/);
unless($ext =~ /\.fa/){$save_fasta_path .="\.faa"}

$query =~ s/SELECT DISTINCT/SELECT DISTINCT Proteins.sequence,/;
&checkDBIconnection();
my $sql = $dbh->prepare($query);
$sql->execute();

$seqio_obj = Bio::SeqIO->new(-file => ">$save_fasta_path", -format => "fasta" );
	while(@result_table = $sql->fetchrow_array()){

	my $seq_obj = Bio::Seq->new();
	$seq_obj->seq(shift @result_table);
	@result_table = map{$_//""} @result_table;
	$seq_obj->display_id(join ("-",@result_table));
	$seqio_obj->write_seq($seq_obj);

	}

}


sub SaveCsv{

my $save_csv_path = shift;
my $query = shift;
my $heading_ref = shift;

my @result_table;
my $seqio_obj;

my ($name,$dir,$ext) = fileparse($save_csv_path, qr/\.[^.]*/);
unless($ext =~ /\.csv/){$save_csv_path .="\.csv"}


&checkDBIconnection();
my $sql = $dbh->prepare($query);
$sql->execute();

open my $fh,">",$save_csv_path or warn "Could not open csv file $!";
print $fh join("\t",@$heading_ref); #headings for columns
print $fh "\n";
	while(@result_table = $sql->fetchrow_array()){
	@result_table = map{$_//""} @result_table;
	print $fh join("\t",@result_table);
	print $fh "\n";

	}
close $fh;
}


###########################################################################################
################ SECONDARY SUBROUTINES ####################################################
################ FOR DATA PROCESSING ######################################################
###########################################################################################
###
###
###
###
###


###########################################################################################
################ SUBROUTINES ##############################################################
################ FOR TAXONOMY PROGRAM #####################################################
###########################################################################################
###
###
###
###
###



###########################################################################################
################ SECONDARY SUBROUTINES ####################################################
################ FOR CLUSTERING        ######################################################
###########################################################################################
###
###
###
###
###




sub findsynteny{
my $type = shift; #'keyword' or 'HMM'
my $queryoption = shift;
my $motif = shift; #HMM motif or keyword name
my $genNum = shift;
my $sql;


if($type eq "keyword"){
my $query = "SELECT DISTINCT cluster_ID FROM Clusters LEFT JOIN Functions ON Clusters.cluster_ID WHERE Functions.keyword $queryoption $motif;";
$sql = $dbh->prepare($query);
$sql->execute();
}elsif($type eq "HMM"){
$sql = $dbh->prepare('SELECT DISTINCT HMM from Proteins where cluster_ID = ? ');
$sql->execute();
}

}
###########################################################################################
################ TERTIONARY SUBROUTINES ###################################################
################ FOR DATA PROCESSING ######################################################
###########################################################################################
###
###
###
###
###

sub checkDBIconnection{
if (!$dbh->ping){
$dbh = $dbh->clone() or die "cannot connect to db";
}
}

sub unpackgz{
	
my $archive = shift;

	#for my $archive(@$archives){

	$archive =~ s/\.gz$// ;
	print "\tUnpacking ",$archive,".gz \n";
	gunzip "$archive.gz" => $archive;
	#}
return $archive;
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

sub getHMMList{
my $HMMLibrary = shift;
my @array;
open my $fh,"<",$HMMLibrary or warn;
while(my $row = <$fh>){
   chomp $row;
   unless($row =~ /^NAME/){next;}
   else{
   $row =~ s/NAME  //;
   push @array,$row;
   }
}
@array = sort { lc($a) cmp lc($b) } @array;
return \@array;
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

sub assignkeyword{
my $cluster = shift;
my $key = shift;

my $mysql = $dbh->prepare("INSERT INTO Functions (cluster_ID, keyword) SELECT ?,? FROM DUAL WHERE NOT EXISTS (SELECT * FROM Functions WHERE cluster_ID= ? AND keyword = ? LIMIT 1);");
$mysql->execute($cluster,$key,$cluster,$key);
}

sub getHeadings{
my $query = shift;
my $selection;
if($query =~ /^SELECT\s(DISTINCT)?(.*?)FROM/i){$selection = $2;} 
$selection =~ s/ //g;
#print $selection,"\n";
my @fields = split ("\,", $selection);

my @headings;

for(@fields){
#print $_," here is the field thing\n";
push (@headings, (split("\\.",$_))[-1] );

}

return \@headings;
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


