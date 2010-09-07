#!/usr/bin/perl -w

use strict;
use OTU;
use File::Spec;
use Getopt::Long;

my ( $samp_path, $sim_read2source_tab );
my $reference_tree = "SSU_BAC_ref_unaln_SSU_BAC_FT_pseudo_pruned_fmt.tree";
my $is_sim = 0;
GetOptions(
    'i=s' => \$samp_path,
    'sim' => \$is_sim,
    'r=s' => \$reference_tree,
    );

if( $is_sim ){
    print "You are running handler for simulations!\n";    
}
print "Processing $samp_path\n";
########################
# USER RUNTIME OPTIONS #
########################

#my @domains   = qw(BAC ARC);
my @domains   = qw( BAC );
my $ecutoff   = 0.000001;
my $trim      = 1; #binary. When on, splices non-homologous sequence linked to 16S read
my $masterdir    = "/Users/sharpton/projects/OTU/db2/"; #upper level directory
my $scripts_path = "/Users/sharpton/projects/OTU/code/"; #where the code is stored
my $seqlencut    = 100;
my $clust_cutoff = 0.05;
my $clust_method = "average";
my $tree_method  = "fasttree";

#######################################
# SET UP FLAT FILE DATABASE STRUCTURE #
#######################################
#Note: sim 1 (prune sim) is testing how read tree compares to source/ref tree
#sim 2 (readsource sim) is testing how reads and sources co cluster on same tree
#sim 3 is prune sim with padded ends; ultrashort and 2 reads/source
my $simtype = 3;

my $project = OTU->new();
if( $is_sim ){
    $project->set_simulation();
}
$project->set_otu_workdir( $scripts_path );
$project->set_sample( $samp_path );
$project->set_domains( \@domains );
$project->set_tree_method( $tree_method );
$project->build_db( $masterdir );
$project->set_blastdb( "BAC", "stap_16S_BAC.fa" );
#$project->set_blastdb( "ARC", "stap_16S_ARC.fa" );
$project->run_SSU_blast( $ecutoff );
$project->load_blastreports;
$project->grab_SSU_reads( $ecutoff, $trim ); 
$project->run_cmalign();
$project->format_alignments();
$project->run_align_qc( $seqlencut );

#turn these on contingent upon simulation needs
#$project->build_read2source_tab();
#$project->prune_alignment( $simtype ); #simulation specific function for sim 1

$project->run_tree(); 
#$project->align_to_seqs();

if( $simtype == 1 ){
  $project->sim_prune_tips( $simtype, $reference_tree ); #for sim 2 and now also sim 1
  $project->sim_tree_to_matrix();
  $project->sim_format_matrix_to_phylip();
  $project->sim_run_mothur( $clust_cutoff, $clust_method );
}
else{
  $project->prune_tips();
  $project->tree_to_matrix();
  $project->format_matrix_to_phylip();
  $project->run_mothur( $clust_cutoff, $clust_method );
}
