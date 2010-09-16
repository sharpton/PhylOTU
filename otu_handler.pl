#!/usr/bin/perl -w

use strict;
use OTU;
use File::Spec;
use Getopt::Long;
use IPC::System::Simple qw(capture $EXITVAL);

my ( $samp_path, $sim_read2source_tab );
my $reference_tree = "SSU_BAC_ref_unaln_SSU_BAC_FT_pseudo_pruned_fmt.tree";
my $is_sim = 0;
my $numblast = 1; #number of parallel blast jobs to run
  # -1 = run only blast and then quit (for parallel jobs)
  # 0  = skip blast and begin the code with run_cmalign
  # 1  = run one blast (non-parallel)
  # >1 = number of parallel blast jobs to run
my $numalnQC = 1; #number of parallel alignment QC jobs to run 
  # -1 = run only alignqc, skip blast, clustering, etc. (for parallel jobs)
  # 1  = run one alignqc (non-parallel)
  # >1 = number of parallel alignQC jobs to run
GetOptions(
    'i=s' => \$samp_path,
    'sim' => \$is_sim,
    'r=s' => \$reference_tree,
    'b=i' => \$numblast,
    'a=i' => \$numalnQC,
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
my $trim      = 2; #Should the reads be trimmed to just the 16S sequence via blast analysis?
                   #0 = no trimming, use complete read
                   #1 = trim to the top HSP coordinates
                   #2 = tile HSPs from top hit, use the best tile, and trim to this tile's coordinates
#my $masterdir    = "/netapp/home/rebeccamae/PhylOTU/db2/"; #upper level directory
my $masterdir    = "/Users/sharpton/projects/OTU/db/"; #upper level directory
#my $scripts_path = "/netapp/home/rebeccamae/PhylOTU/code/"; #where the code is stored
my $scripts_path = "/Users/sharpton/projects/OTU/PhylOTU"; #where the code is stored
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

#######################################
# Now get to work
#######################################

if( $numalnQC != -1 && $numblast != 0){  # -1 means run only alignQC
  if( $numblast > 1 ){
    parallel( "blast" );
  } elsif( abs($numblast) == 1 ){  # 1 normal, -1 parallel mode
    $project->run_SSU_blast( $ecutoff );
    $project->load_blastreports;
    $project->grab_SSU_reads( $ecutoff, $trim ); 
  }
  if( $numblast == -1 ){           # parallel mode, only run blast then quit
    exit;
  }
  $project->run_cmalign();
  $project->format_alignments();
}

if( $numalnQC > 1 ){
  parallel ( "alignQC" );
} else {                           # 1 normal, -1 parallel mode
  $project->run_align_qc( $seqlencut, $numalnQC );
  if( $numalnQC == -1 ){           # parallel mode, only run alignQC then quit
    exit;
  }
}

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

##############################################################################
# Send out scripts to run the blast step in parallel
# The queue commands may need to be modified on other systems
##############################################################################
sub parallel
{
  # Split up the sample into $numblast subsamples and blast/alignQC all the subsamples in parallel

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
	} elsif (`grep -i error $outfile`){
	  # code crashed, quit (or I could change this to retry....)
	  print $type." crashed or something else went wrong with file ".$i.", I quit\n";
	  exit(0);
	} else {
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
  my $command = "./run_otu_handler.sh -i ".$subsample;
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
