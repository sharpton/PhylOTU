#!/usr/bin/perl -w

use strict;
use OTU;
use File::Spec;
use Getopt::Long;
use IPC::System::Simple qw(capture $EXITVAL);

my ( $samp_path, $sim_read2source_tab );
#This is an optional simulation specific variable; not necessary for analysis of metagenomic libraries
my $reference_tree = "/netapp/home/sharpton/projects/OTU/db/reference/trees/SSU_BAC_ref_unaln_SSU_BAC_FT_pseudo_pruned_fmt.tree";
my $is_sim     = 0;
my $numblast   = 1;
  # -1 = run only blast (for parallel jobs)
  # 1  = run one blast
  # >1 = number of parallel blast jobs to run
my $numalnQC   = 1;
  # -1 = run only alignqc (for parallel jobs)
  # 1  = run one alignqc
  # >1 = number of parallel alignQC jobs to run
my $numTreeMat = 1;
  # -1 = run only tree_to_matrix (for parallel jobs)
  # 1  = run one tree_to_matrix
  # >1 = number of parallel tree_to_matrix jobs to run
my $startMat   = 0;  # Only needs to be specified when parallel tree_to_matrix jobs are running
my $endMat     = 0;  # Only needs to be specified when parallel tree_to_matrix jobs are running
my $T2M_domain = ""; # Only needs to be specified when parallel tree_to_matrix jobs are running
my $skipto     = ""; # skip to (A)lignment, alignment(Q)c, (T)ree, tree2(M)atrix, (C)lustering 
my $clust_cutoff = 0.15;
GetOptions(
    'i=s' => \$samp_path,
    'sim' => \$is_sim,
    'r=s' => \$reference_tree,
    'b=i' => \$numblast,
    'a=i' => \$numalnQC,
    'm=i' => \$numTreeMat,
    's=i' => \$startMat,
    'e=i' => \$endMat,
    'd=s' => \$T2M_domain,
    'k=s' => \$skipto,
    'cc:f'=> \$clust_cutoff,
    );
if( $numTreeMat > 0 && ( $startMat!=0 || $endMat!=0 )){
  print "tree_to_matrix not running in parallel mode, start and end step to full matrix printing\n";
  $startMat = 0;
  $endMat   = 0;
}
if( $is_sim ){
    print "You are running handler for simulations!\n";    
}
my $do_pruning = $numTreeMat;
print "Processing $samp_path\n";
########################
# USER RUNTIME OPTIONS #
########################

my @domains   = qw(BAC ARC);
#my @domains   = qw( BAC );
my $ecutoff   = 0.000001;
my $trim      = 2; #Should the reads be trimmed to just the 16S sequence via blast analysis?
                   #0 = no trimming, use complete read
                   #1 = trim to the top HSP coordinates
                   #2 = tile HSPs from top hit, use the best tile, and trim to this tile's coordinates
my $masterdir    = "/netapp/home/rebeccamae/PhylOTU/db2/"; #upper level directory
#my $masterdir    = "/netapp/home/sharpton/projects/OTU/db/"; #upper level directory
my $scripts_path = "/netapp/home/rebeccamae/PhylOTU/code/"; #where the code is stored
#my $scripts_path = "/netapp/home/sharpton/projects/OTU/PhylOTU/"; #where the code is stored
my $seqlencut    = 100;
my $mincoverage  = 2;
#my $clust_cutoff = 0.15;#moved above!
my $clust_method = "average";
my $tree_method  = "fasttree";
my $dist_cutoff  = 2*$clust_cutoff;  # used to print shortened tree_to_matrix list for ESPRIT, bigger than the clustering threshold for safety
my $tree_to_matrix_code = "c";
  # c = c++ code, output formatted for mothur (will implement ESPRIT if needed)
  # R = R script, output formatted for mothur
my $cluster_code        = "E";
  # M = mothur, requires matrix format from tree_to_matrix (R or c++)
  # E = ESPRIT, requires distance list and frequency list (c++ code only)
if( ($tree_to_matrix_code eq 'R') && ($cluster_code ne 'M') ){
  print "ERROR: clustering method ESPRIT requires output format not implemented in R, quitting\n";
  exit;
}
#######################################
# SET UP FLAT FILE DATABASE STRUCTURE #
#######################################
#Note: sim 1 (prune sim) is testing how read tree compares to source/ref tree
#sim 2 (readsource sim) is testing how reads and sources co cluster on same tree
#sim 3 is prune sim with padded ends; ultrashort and 2 reads/source
my $simtype = 1;
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
$project->set_blastdb( "ARC", "stap_16S_ARC.fa" );

#######################################
# Jump to the appropriate step
#######################################
if( $numalnQC   == -1  ){ print "Skipping to alignmentQC\n";    goto ALIGNQC; }
if( $numTreeMat == -1  ){ print "Skipping to tree_to_matrix\n"; goto TREE2MATRIX; }
if( $skipto     eq "A" ){ print "Skipping to alignment\n";      goto ALIGN; }
if( $skipto     eq "Q" ){ print "Skipping to alignentQC\n";     goto ALIGNQC; }
if( $skipto     eq "L" ){ print "Skipping to alignentQC col\n"; goto ALIGNQCCOL; }
if( $skipto     eq "T" ){ print "Skipping to tree creation\n";  goto TREE; }
if( $skipto     eq "M" ){ print "Skipping to tree_to_matrix\n"; goto TREE2MATRIX; }
if( $skipto     eq "C" ){ print "Skipping to clustering\n";     goto CLUSTER; }

#######################################
# Now get to work
#######################################

if( $numblast > 1 ){
  parallel( "blast" );
} else {
  $project->run_SSU_blast( $ecutoff );
  $project->load_blastreports;
  $project->grab_SSU_reads( $ecutoff, $trim ); 
  if( $numblast == -1 ){ exit; }  # parallel mode, only run blast then quit
}

ALIGN:
$project->run_cmalign();
$project->format_alignments();

ALIGNQC:
if( $numalnQC > 1 ){
  parallel ( "alignQC" );
} else {
  $project->run_align_qc( $seqlencut, $numalnQC );
  if( $numalnQC == -1 ){ exit; } # parallel mode, only run alignQC then quit
}
ALIGNQCCOL:
$project->run_align_ColQC( $mincoverage );

#turn these on contingent upon simulation needs
#$project->build_read2source_tab();
#$project->prune_alignment( $simtype ); #simulation specific function for sim 1

TREE:
$project->run_tree(); 
#$project->align_to_seqs();

TREE2MATRIX:
if( $is_sim && $simtype == 1 ){
#  $project->sim_prune_tips( $simtype, $reference_tree ); #for sim 2 and now also sim 1
#  $project->sim_tree_to_matrix();
#  $project->sim_format_matrix_to_phylip();
  $project->sim_tree_to_matrix_cpp( $simtype, $reference_tree, $startMat, $endMat, $cluster_code, $clust_cutoff, $do_pruning );
  $project->sim_run_mothur( $clust_cutoff, $clust_method );
} else {

  if( $tree_to_matrix_code =~ 'R'){
    $project->prune_tips();
    # Pruning in cpp version is done in tree_to_matrix, via the flag do_pruning>0
  }

  if( $numTreeMat > 1 ){
    $project->tree_to_matrix_cpp($startMat, $endMat, $cluster_code, $dist_cutoff, $do_pruning, $T2M_domain);
    parallel("tree_to_matrix");
  } else {
    if(      $tree_to_matrix_code =~ 'R' ){
      $project->tree_to_matrix_R();
      $project->format_matrix_to_phylip();
    } elsif( $tree_to_matrix_code =~ 'c' ){
      $project->tree_to_matrix_cpp($startMat, $endMat, $cluster_code, $dist_cutoff, $do_pruning, $T2M_domain);
    } else {
      print "Unknown tree_to_matrix_code version: $tree_to_matrix_code, quitting\n";
      exit;
    }
    if( $numTreeMat == -1 ){ exit; } # parallel mode, only run tree_to_matrix then quit
  }

CLUSTER:
  if(      $cluster_code eq 'M' ){
    $project->run_mothur( $clust_cutoff, $clust_method );
  } elsif( $cluster_code eq 'E' ){
    $project->run_Ecluster();
  }
}

##############################################################################
# Send out scripts to run steps in parallel
# The queue commands may need to be modified on other systems
##############################################################################
sub parallel
{
  # Split up the sample into subsamples and blast/alignQC/tree_to_matrix all the subsamples in parallel

  my $type = shift;
  my @subsample_files;
  my $option;
  my @qID;

  if( $type =~ "blast"){
    $option = " -b -1";
    @subsample_files = $project->split_query($type, $numblast);
  } elsif($type =~ "alignQC") {
    @subsample_files = $project->split_query($type, $numalnQC);
    $option = " -a -1";
  } elsif($type =~ "tree_to_matrix") {
    # The subsample_files text here already include the options -s START -e END -d DOMAIN
    @subsample_files = $project->split_tree($numTreeMat);
    $option = " -m -1";
  }

  # loop over each subsample of queries and submit them
  foreach my $subsample ( @subsample_files ){
    push(@qID, submit($subsample.$option) );
  }
  my $jobs_running = scalar(@qID);
  
  # Check if some or all jobs are finished
  while ( $jobs_running > 0 ){
    # Something is still running, so wait another minute
    sleep(60);
    # loop through jobs to see if any are done
    for( my $i=0; $i<scalar(@qID); $i++ ){
      # I'm only interested in files that are done but I haven't flagged them as 0 yet
      if( $qID[$i]!=0  &&  !`qstat | grep $qID[$i]` ){
	# This job recently finished, check if everything ran ok
	my $outfile = "logs/".$qID[$i].".all";
	if( !(`tail -1 $outfile` =~ m/RUN FINISHED/ )){
	  # The node crashed, resubmit this job
	  print "Rerunning file ".$i." as it seems the node failed\n";
	  $qID[$i] = submit($i.$option);
	} elsif (`grep -i kill $outfile`){
	  # code crashed, quit (or I could change this to retry....)
	  print $type." crashed or something else went wrong with file ".$i.", I quit\n";
	  exit(0);
	} else {
	  if (`grep -i error $outfile`){
	    print $type." crashed or something else went wrong with file ".$i.", but I will continue\n";
	  }
	  # Everything ran fine, remove this job from the list
	  print "Job ".$qID[$i]." on subsample ".$i." finished sucessfully\n";
	  $qID[$i] = 0;
	  $jobs_running--;
	}
      }
    }
  }
  
  # Stich all the results back together and put them in the expected path/file
  if( $type =~ "blast"){
    $project->stitch_blast(@subsample_files);
  } elsif($type =~ "alignQC") {
    $project->stitch_alignQC(@subsample_files);
  } elsif($type =~ "tree_to_matrix") {
    $project->stitch_matrix(@subsample_files);
  }
}


##############################################################
# Submit a blast job to the queue
# The queue commands may need to be modified on other systems
##############################################################
sub submit
{
  my $subsample = shift;
  
  # the run_otu_handler_parallel.sh script must be configured properly for your batch system
  my $command = "./run_otu_handler_small.sh -i ".$subsample;
  my $result = capture( "qsub $command" );
  print( "qsub $command" );
  # Should return a string like 
  # Your job 4827092 ("run_otu_handler_parallel.sh") has been submitted
  # grab the job # only
  my @values = split(' ', $result);
  if( $values[2] =~ /\d/ ){
    print( "\t JOB_ID=$values[2]\n");
    return $values[2];
  } else {
    return 0;
  }
  
}
