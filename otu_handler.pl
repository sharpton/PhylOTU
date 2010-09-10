#!/usr/bin/perl -w

use strict;
use OTU;
use File::Spec;
use Getopt::Long;
use IPC::System::Simple qw(capture $EXITVAL);

my ( $samp_path, $sim_read2source_tab );
my $reference_tree = "SSU_BAC_ref_unaln_SSU_BAC_FT_pseudo_pruned_fmt.tree";
my $is_sim = 0;
my $numpara = 1; #number of parallel blast jobs to run
GetOptions(
    'i=s' => \$samp_path,
    'sim' => \$is_sim,
    'r=s' => \$reference_tree,
    'p=i' => \$numpara,
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
my $masterdir    = "/netapp/home/rebeccamae/PhylOTU/db2/"; #upper level directory
my $scripts_path = "/netapp/home/rebeccamae/PhylOTU/code/"; #where the code is stored
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
if( $numpara > 1 ){
    parallel_blast( $ecutoff );
} else {
    $project->run_SSU_blast( $ecutoff );
    $project->load_blastreports;
    $project->grab_SSU_reads( $ecutoff, $trim ); 
}
if( $numpara < 1 ){
    # Negative values are a flag indicating this is 
    # being run in parallel for only the blast portion
    # Now that that is finshed, quit
    # (This is done below in the submit() function)
    exit;
}
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


##############################################################################
# Send out scripts to run the blast step in parallel
# The queue commands may need to be modified on other systems
##############################################################################
sub parallel_blast
{

    # Split up the sample into $numpara subsamples and blast all the subsamples in parallel
    my @subsample_files = $project->split_query($numpara);

    # loop over each subsample of queries and submit them
    my @qID;
    foreach my $subsample ( @subsample_files ){
	push(@qID, submit($subsample) );
    }
    my $jobs_running = scalar(@qID);

    # Check if some or all jobs are finished
    while ( $jobs_running > 0 ){
	# Something is still running, so wait 5 more minutes
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
		    $qID[$i] = submit($i);
		} elsif (`grep -i error $outfile`){
		    # blast crashed, quit (or I could change this to retry....)
		    print "Blast crashed or something else went wrong with file ".$i.", I quit\n";
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
    $project->stitch_blast(@subsample_files);
    
}


##############################################################
# Submit a blast job to the queue
# The queue commands may need to be modified on other systems
##############################################################
sub submit
{
    my $subsample = shift;

    # the run_otu_handler_parallel.sh script must be configured properly for your batch system
    my $command = "./run_otu_handler_parallel.sh -i ".$subsample." -p -1";
    my $result = capture( "qsub $command" );
    # Should return a string like 
    # Your job 4827092 ("run_otu_handler_parallel.sh") has been submitted
    # grab the job # only
    my @values = split(' ', $result);
    if( $values[2] =~ /\d/ ){
	return $values[2];
    } else {
	return 0;
    }

}
