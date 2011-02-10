#!/usr/bin/perl -w

#OTU.pm - The PhylOTU workflow manager
#Copyright (C) 2011  Thomas J. Sharpton 
#author contact: thomas.sharpton@gladstone.ucsf.edu
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#    
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#    
#You should have received a copy of the GNU General Public License
#along with this program (see LICENSE.txt).  If not, see 
#<http://www.gnu.org/licenses/>.

package OTU;

use strict;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::TreeIO;
use Bio::SearchIO;
use Bio::Search::Tiling::MapTiling;
use Bio::Phylo::IO 'parse';
use File::Basename;
use File::Path qw( make_path );
use File::Copy;
use IPC::System::Simple qw(capture $EXITVAL);
use Data::Dumper;
use POSIX;

=head2 new

 Title   : new
 Usage   : $project = OTU->new()
 Function: initializes a new OTU analysis project
 Example : $project = OTU->new();
           $project->set_sample("Sample_Name", "Sample_Data_Path"); #Add a sample
           $project->build_db( "database_path" ); #create flat file database
           $project->process_sample(); #run the entire OTU pipeline
 Returns : OTU project object
 Args    : None

=cut

sub new{
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self  = {};
    $self->{"blastdb"}    = undef; #hashref:  blast seq db filenames
    $self->{"breports"}   = undef; #hashref:  blasthit report filenames
    $self->{"bhits"}      = undef; #hashref:  seqs with blasthit filenames
    $self->{"domains"}    = undef; #arrayref: sets or phylogenetic domain labels
    $self->{"sample"}     = undef; #hashref:  sample file, full path, and name
    $self->{"reads"}      = undef; #check if depricated, delete
    $self->{"db"}         = undef; #hashref: a hash of flat file database dir locations
    $self->{"workdir"}    = undef; #scalar: path to OTU scripts
    $self->{"treemethod"} = undef; #fasttree, xrate, raxml, etc.
    $self->{"simulation"} = 0;     #binary: turn on for simulation specific processing
    bless($self);
     return $self;
}

=head2 set_otu_workdir

 Title   : set_otu_workdir
 Usage   : $project->set_otu_workdir( "path_to_scripts" )
 Function: Stores file name and location of sample currently being processed
 Example : my $scripts_path = "/projects/OTU/scripts/";
           $project->set_otu_workdir( $scripts_path )
 Returns : The path to the work directory for this pipeline (scripts location)
 Args    : Scalar of full path to the directory containing the OTU pipeline scripts

=cut

sub set_otu_workdir{
  my $self = shift;
  my $path = shift;
  $self->{"workdir"} = $path;
  return $self->{"workdir"};
}

sub set_tree_method{
  my $self   = shift;
  my $method = shift;
  $self->{"treemethod"} = $method;
  return $self;
}

=head2 set_sample

 Title   : set_sample
 Usage   : $project->set_sample("samp_name", "samp_data_path")
 Function: Stores file name and location of sample currently being processed
 Example : my $samp_data = "/data/OTU/Test_Reads.fa";
           $project->set_sample( $samp_data )
 Returns : OTU project object with sample information stored
 Args    : Full file path of sample data

=cut

sub set_sample{
    my ( $self, $path )    = @_;
    my ($samplenm, $dir, $suffix) = fileparse( $path, qr/\.[^.]*/ );
    $self->{"sample"}->{"file"}   = $samplenm . $suffix;
    $self->{"sample"}->{"name"}   = $samplenm;
    $self->{"sample"}->{"path"}   = $path;  #full path to raw data, inc. file
    return $self;
}

=head2 set_domains

 Title   : set_domains
 Usage   : $project->set_domains( @domains )
 Function: Stores an array of unique labels which indicate the particular
           phylogenetic domains (or sets) to partition the analysis across.
           Note that these labels must appear in the reference data filenames.
 Example : my @domains   = qw(BAC ARC EUK);
           $project->set_domains( @domains );
 Returns : OTU project object with sample information stored
 Args    : An array of domain (set) identifiers

=cut

sub set_domains{
    my ($self, $refarray ) = @_;
    $self->{"domains"} = $refarray;
    return $self;
}

=head2 build_db

 Title   : build_db
 Usage   : $project->build_db( $masterdir )
 Function: Creates the directory structure for the flat file database.
           This database serves as the input/output point for all 
           OTU analysis procedures, thus this method should be one of the 
           first executed in any custom pipeline. This database structure 
           is rigid, thus manipulating the files and directories within
           it may break the analysis pipeline.
 Example : my $masterdir = "/projects/OTU/db/";
           $project->set_domains( $masterdir );
 Returns : OTU project object with a stored lookup table of database paths
 Args    : A scalar of the path to the directory within which the database
           will be stored

=cut

sub build_db{
    my ( $self, $masterdir ) = @_;    
    if( !( -e $masterdir ) ){
	mkdir( $masterdir );
    }
    my $sample    = $self->{"sample"}->{"name"};
    my $sampledir = $masterdir . "samples/" . $sample;
    $self->{"db"} = {
	#BLASTDBS
	blastdb   => $masterdir . "/blastdbs/",
	blastout  => $masterdir . "/blastout/",
	#REFERENCE DATA
	profile   => $masterdir . "/reference/profiles/",
	ref_align => $masterdir . "/reference/aligns/",
	#SAMPLE DATA        
	reads     => $sampledir . "/raw/",
	SSU_reads => $sampledir . "/SSU/",
	all_align => $sampledir . "/aligns/raw/",
	cm_scores => $sampledir . "/aligns/scores/",
	qc_align  => $sampledir . "/aligns/qc/",
	qc_seqs   => $sampledir . "/all_qc_seqs/",
	tree      => $sampledir . "/trees/",
	matrix    => $sampledir . "/matrix/",
	otudir    => $sampledir . "/otus/",
    };
    foreach my $path ( keys( %{ $self->{"db"} } ) ){
	make_path( $self->{"db"}->{$path} );
    }
    #push raw sample data into proper db spot - for now we copy!
    unless( $self->{"sample"}->{"path"} eq $self->{"db"}->{"reads"} . $self->{"sample"}->{"file"} ){
	copy ($self->{"sample"}->{"path"} , $self->{"db"}->{"reads"} . $self->{"sample"}->{"file"} );
    }
    return $self->{"db"};
}


=head2 set_blastdb

 Title   : set_blastdb
 Usage   : $project->set_blastdb( "domain_label", "SSU_blast_database_file" )
 Function: Maps domain (set) label to a blastdb which as *already* been placed
           into the flat file database and properly formatted with formatdb.
           Note that this script will not do these last two steps on your behalf.
           Also note that you do not need to indicate the full path of the database,
           only the name of the file as it is found in the flat file db blastdb dir.
           At the moment, there is no batch run - you need to run this method 
           iteratively to map blastdbs to multiple different sets.
 Example : my $domain_label = "BAC";
           my $SSU_blast_database_file = "stap_16S_BAC.fa";
           $project->set_blastdb( "BAC", "stap_16S_BAC.fa");
 Returns : OTU project object with this set to blastdb mapping stored
 Args    : A scalar of the set label as defined in @domains and
           a scalar of the name of the blastdb that maps to this label. 

=cut

sub set_blastdb{
    my ( $self, $set, $dbname ) = @_;
    $self->{"blastdb"}->{$set} = $dbname;
    return $self;
}

=head2 format_blast_db

 Title   : format_blast_db
 Usage   : $project->format_blast_db
 Function: Format nt FASTA file into blast database. Here we use it to 
           format the STAP databases. Autoskips if this is already done.
 Example : $project->format_blast_db
 Returns : None
 Args    : None

=cut

sub format_blast_db{
  my $self = shift;
  for my $set ( @{ $self->{"domains"} } ){
    my $sourcepath  = $self->{"workdir"} . "data/stap_16S_" . "$set" . ".fa";
    my ( $source, @suffixlist) = fileparse( $sourcepath );
    my $database    = $self->{"db"}->{"blastdb"}  . $self->{"blastdb"}->{$set};
    if( -e $database ){
      warn "Blast database exists at $database, so skipping formatdb\n";
      next;
    }
    else{
      copy ( $sourcepath, $database );
      print "Formatting blast database " . $self->{"blastdb"}->{$set} . "\n";
      my @args = ("-i $database", "-p F");
      my $results = capture( "formatdb @args" );
      if( $EXITVAL != 0 ){
	warn("Error running formatdb in format_blast_db!\n");
	exit(0);
      }
    }
  }
  return $self;
}

=head2 build_cm_model

 Title   : format_blast_db
 Usage   : $project->build_cm_model( $set, $ref_aln )
 Function: Builds a CMmodel given a reference stockholm alignment. Copies the alignment
           to the $masterdir/reference/aligns/ directory. Dumps model into 
           $masterdir/reference/profiles/ directory. 
 Example : $project->build_cm_model( "BAC", "SSU_BAC_ref.stk" );
 Returns : None
 Args    : The name of the domain for which the model is being built and the absolute path
           to the location of the reference stockholm alignment

=cut

sub build_cm_model{
  my $self    = shift;
  my $ref_aln = shift;
  my $set     = shift; 
  #Come back in and grab the ref alignment name to enable multiple models in single db
  #my ( $name, $path, $suffix ) = fileparse( $ref_aln, qr/\.[^.]*/ );
  my $db_ref_aln =  $self->{"db"}->{"ref_align"} . "SSU_" . $set . "_ref.stk";
  if( -e $db_ref_aln ){
    print "A reference alignment with the same name as the one you're using already exists in the PhylOTU datase at $db_ref_aln. Please remove this alignment before proceeding and then rerun. I don't want to overwrite a file you might need.\n";
    die;
  }
  copy ( $ref_aln , $db_ref_aln );
  my $model = $self->{"db"}->{"profile"} . "SSU_" . $set . ".mod";
  if( -e $model ){
    print "A CMmodel with the same name as the one you are trying to build already exists in the PhylOTU datase at $model. Please remove this model before proceeding and then rerun. I don't want to overwrite a file you might need.\n";
    die;
  }
  print "Building CMmodel $model...\n";
  my @args = ("--rf", "--ere 1.4", "$model", "$db_ref_aln");
  my $results = capture( "cmbuild @args" );
  if( $EXITVAL != 0 ){
    warn("Error running cmbuild in build_cm_model!\n");
    exit(0);
  }
  return $self;
}

=head2 load_blastreports

 Title   : load_blastreports
 Usage   : $project->load_blastreports
 Function: A convienence tool that loads a read V. SSUdb blastreport for 
           downstream processing. Prevents having to run blast again.
 Example : $project->load_blastreports
 Returns : OTU project object with the blast reports stored
 Args    : None

=cut

sub load_blastreports{
    my $self = shift;
    for my $set ( @{ $self->{"domains"} } ){
	my $breport = $self->{"db"}->{"blastout"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . ".blast";
	print "loading report for $set at $breport\n";
	$self->{"breports"}->{$set} = $breport;
    }
    return $self;
}

=head2 run_SSU_blast

 Title   : run_SSU_blast
 Usage   : $project->run_SSU_blast();
 Function: A blastn wrapper. Searches sample reads against each sets SSU blastdb.
 Example : $project->run_SSU_blast();
 Returns : None
 Args    : None

=cut

sub run_SSU_blast{
    my $self      = shift;
    my $blaste    = shift;
    my $blasttype = "blastn";
    my $input     = $self->{"db"}->{"reads"} . $self->{"sample"}->{"file"};
    my $set_ct    = 0;
    foreach my $set ( keys( %{ $self->{"blastdb"} } ) ){
	$set_ct++;
	my $database    = $self->{"db"}->{"blastdb"}  . $self->{"blastdb"}->{$set};
	my $blastreport = $self->{"db"}->{"blastout"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . ".blast";
	$self->{"breports"}->{$set} = $blastreport;
	my @args = ("-p blastn", "-i $input", "-d $database", "-o $blastreport", "-m 8", "-e $blaste");
	my $results = capture( "blastall @args" );
	if( $EXITVAL != 0 ){
	    warn("Error running $blasttype in run_SSU_blast!\n");
	    exit(0);
	}
    }
    return $set_ct;
}


=head2 grab_SSU_reads

 Title   : grab_SSU_reads
 Usage   : $project->grab_SSU_reads( $ecutoff, $trim )
 Function: Open blastoutputs for each set and determine which sample reads
           exhibit significant similarity to SSU blast database sequences.
           The hit score is used to determine which set a read belongs 
           to; sample reads can only belong to one set! If trim is 1 or 2 (not 0), 
           then reads with non-SSU-homologous sequences are trimmed down 
           to just the homologous subsequence. 
           Trim = 0 results in no triming.
           Trim = 1 uses the top HSP coordinates to trim the sequence. 
           Trim = 2 uses the top HSP tiling coordinates to trim.           
 Example : my $ecutoff =  0.00001;
           my $trim    = 1;
           $project->grab_SSU_reads( $ecutoff, $trim );
 Returns : 1 upon success
 Args    : A scalar of the evalue cutoff used to determine significance and
           a scalar of whether to employ read trimming and what trimming method

=cut

#NOTE: need to build filter to ensure that if no hits are found, don't keep processing
#the domain

sub grab_SSU_reads{
    my ($self, $ecutoff, $qc) = @_;
    foreach my $set ( keys( %{ $self->{"breports"} } ) ){
	my $blastreport = $self->{"breports"}->{$set};
	print "Processing BLAST results for $set found in $blastreport\n";
#	my $in = new Bio::SearchIO(-format => 'blast', -file => $blastreport);
	my $in = new Bio::SearchIO(-format => 'blasttable', -file => $blastreport);
        RESULT: while( my $result = $in->next_result ) {
	    while( my $hit = $result->next_hit ) {
	      if( $qc == 2 ){
		#evaluate the top hsp tile from the top blast hit
		my $tiling = Bio::Search::Tiling::MapTiling->new($hit);
		my @top_tile_hsps = $tiling->next_tiling( 'query' );
		#process the top hsp to initialize the tile
		my $top_hsp = pop( @top_tile_hsps );
		next unless( $top_hsp->evalue() < $ecutoff );
		my $qstart = $top_hsp->start('query');
		my $qstop  = $top_hsp->end('query');
		my $score  = $top_hsp->bits;
		my $qstrand = $top_hsp->strand('query');
		my $hstrand = $top_hsp->strand('hit');
		my $nhsps   = 1;
		#Uncomment the below for troubleshooting
		#print join( "\t", "." , ".", $result->query_name, $hit->name, $top_hsp->start('query'), $top_hsp->end('query'), $qstart, $qstop, "\n");		  
		#loop through the remaining hsps in the tile
		foreach my $hsp( @top_tile_hsps ){
		  next unless ( $hsp->evalue() < $ecutoff );
		  #update the query coordinates for trimming
		  if( $qstart > $hsp->start('query') ){
		    $qstart = $hsp->start('query');
		  }
		  if( $qstop < $hsp->end('query') ){
		    $qstop = $hsp->end('query');
		  }
		  $score = $score + $hsp->bits;
		  #make sure the tile mapping kept consistent tile context (strands are the same across all hsps w/in tile)
		  if( ( $hsp->strand('query') != $qstrand ) || ( $hsp->strand('hit') != $hstrand) ){
		    warn( "ERROR IN grab_SSU_reads: Strands differ between HSPs in the same tile, passing analysis on " . 
			  $result->query_name. " (top hit was " . $hit->name . ")\n");
		    next RESULT;
		  }
		  #Uncomment the below for troubleshooting
		  #print join( "\t", "." , ".", $result->query_name, $hit->name, $hsp->start('query'), $hsp->end('query'), $qstart, $qstop, "\n");	  
		}
		my $avg_score = $score / $nhsps;
		$self->{"bhits"}->{$result->query_name} = { 
							   'set'     => $set,
							   'qstart'  => $qstart,
							   'qstop'   => $qstop,
							   'score'   => $avg_score,
							   'qstrand' => $qstrand,
							   'hstrand' => $hstrand,
							  };
		next RESULT;
	      }
	      else{
		while( my $hsp = $hit->next_hsp ) {
		  #we're only going to evaluate the top hit
		  if( $hsp->evalue() < $ecutoff ){ 
		    my $qstart  = $hsp->start('query');
		    my $qstop   = $hsp->end('query');
		    #			my $score   = $hit->raw_score;
		    my $score   = $hsp->bits;
		    my $qstrand = $hit->strand("query");
		    my $hstrand = $hit->strand("hit");
		    #TURNED OFF FOR BLASTTABLE FORMATTING PURPOSES
		    #my $qlen   = $result->query_length;
		    #if($qstart > $qlen || $qstop > $qlen){
		    #die("HSP query coordinates out of bounds: query length is $qlen, while qhsp start is $qstart and qhsp stop is $qstop for ". $result->query_name . "\n"
		    #);
		    #}
		    #grab the passing items, store the best set's hsp results
		    if( defined( $self->{"bhits"}->{$result->query_name} ) ){
		      print Dumper ($self->{"bhits"}->{$result->query_name} );
		      print "set is $set\n";
		      print " comparing results for " . $result->query_name . "\n";
		      my $old_score = $self->{"bhits"}->{$result->query_name}->{"score"};
		      print "  score: $score\n";
		      print "  old:   $old_score\n";
		      if ($score > $old_score ){
			$self->{"bhits"}->{$result->query_name} = { 
								   'set'     => $set,
								   'qstart'  => $qstart,
								   'qstop'   => $qstop,
								   'score'   => $score,
								   'qstrand' => $qstrand,
								   'hstrand' => $hstrand,
								  };
		      }
		    }
		    else{
		      $self->{"bhits"}->{$result->query_name} = { 
								 'set'     => $set,
								 'qstart'  => $qstart,
								 'qstop'   => $qstop,
								 'score'   => $score,
								 'qstrand' => $qstrand,
								 'hstrand' => $hstrand,
								};
		    }
		    next RESULT;
		  }
		}
	      }
	    }
	  }
      }
    #Grab reads that pass thresholds
    print "Grabbing passing libary reads\n";
    if($qc == 1 ){
      print "Quality control is ON; trimming to top HSP coordinates\n";
    }
    if($qc == 2){
      print "Quality control is ON; tiling top hit HSPs, trimming linked coordinates\n";
    }
    my $seq_in = Bio::SeqIO->new(-format  => 'fasta', 
				 -file  => $self->{"db"}->{"reads"} . $self->{"sample"}->{"file"},
	);
    my %outseqs = ();
    foreach my $set( @{ $self->{"domains"} } ){
	my $output = $self->{"db"}->{"SSU_reads"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . ".fa";
	my $seq_out = Bio::SeqIO->new(-format => 'fasta', -file  => ">$output");
	$outseqs{$set} = $seq_out;
    }
    while( my $seq = $seq_in->next_seq() ){
	my $id = $seq->display_id();
	#NEED TO FORMAT FOR PATTERN MATCHING, but not for exact match...?
	#$id =~ s/\|/\\\|/g;
	#$id =~ s/\./\\\./g;
	foreach my $passing_id( keys( %{ $self->{"bhits"} } ) ) { 
	    my %hsp_data = %{ $self->{"bhits"}->{$passing_id} };
	    my $outset   = $hsp_data{"set"};
	    #Consider if this should be exact match or not...
	    if( $passing_id eq $id ){
		#print "match at $passing_id and $id\n";
		if($qc == 1 || $qc == 2){
		    #grab only homologous subsequences of $seq from top HSP
		    my $qstart = $hsp_data{'qstart'};
		    my $qstop  = $hsp_data{'qstop'};
		    #ADDED THE LENGTH CHECK HERE DUE TO BLASTTABLE, SEE ABOVE
		    my $qlen   = $seq->length();
		    if($qstart > $qlen || $qstop > $qlen){
			die("HSP query coordinates out of bounds: query length is $qlen, while qhsp start is $qstart and qhsp stop is $qstop for ". $seq->display_id() . "\n"
			    );
		    }
		    if($hsp_data{'qstrand'} == -1 ){
			$seq = $seq->revcom;
		    }
		    my $subseq = $seq->subseq($qstart, $qstop);
		    my $seqcp  = $seq;
		    $seqcp->seq($subseq);
		    if( $hsp_data{'hstrand'} == -1 ){
			$seqcp = $seqcp->revcom;
		    }
		    $outseqs{$outset}->write_seq($seqcp);
		}
		else{
		    if( $hsp_data{'hstrand'} == -1 ||  $hsp_data{'qstrand'} == -1 ){
			$seq = $seq->revcom;
		    }
		    $outseqs{$outset}->write_seq($seq);
		}
	    }
	}
    }    
}

sub reset_domains_by_SSU_size{
  my $self = shift;
  my @domains = @{ $self->{"domains"} };
  my @reset_domains = ();
  foreach my $set ( @domains) {
    my $ssufile = $self->{"db"}->{"SSU_reads"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . ".fa";
    my $filesize = -s $ssufile;
    if( $filesize == 0 ){
      print "There were no 16S rRNA sequences in your sample for the domain set $set. This domain will be ignored for all downstream analyses.\n";
      next;
    }
    else{
      print "Found 16S hits for the domain set $set.\n";
      push( @reset_domains, $set );
    }
  }
  $self->set_domains( \@reset_domains );
  return( \@reset_domains );
}


=head2 run_cmalign

 Title   : run_cmalign
 Usage   : $project->run_cmalign();
 Function: A cmalign (INFERNAL) wrapper. Aligns sample SSU reads for each set to the
           sets reference profile, merging the reference alignment to the output.
           The resulting alignment, alignment score table are stored
 Example : $project->run_cmalign();
 Returns : None
 Args    : None

=cut

sub run_cmalign{
    my ( $self )  = @_;
    my @domains   = @{ $self->{"domains"} };
    my $set_ct    = 0;
    foreach my $set ( @domains) {
        print "Aligning 16S reads to CMmodel for domain set $set...\n";
	my $sequences = $self->{"db"}->{"SSU_reads"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . ".fa";
	next unless( -e $sequences );
	my $outaln    = $self->{"db"}->{"all_align"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all.sto";
	my $cmscores  = $self->{"db"}->{"cm_scores"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_cmscores.tab";
	my $cmfile    = $self->{"db"}->{"profile"}   . "SSU_" . $set . ".mod";
	my $refaln    = $self->{"db"}->{"ref_align"} . "SSU_" . $set . "_ref.stk";
	my @args = ("--rf", "--dna","--hbanded","--sub", "-o $outaln", "--withali", "$refaln", "$cmfile", "$sequences");
#	open (SCORES, ">$cmscores" ) || die "Can't open $cmscores for read:$!\n";
	my $results = capture( "cmalign @args > $cmscores" );
	if( $EXITVAL != 0 ){
	    warn("Error running cmalign for $sequences against $cmfile using $refaln!\n");
	    exit(0);
	}
#	print SCORES $results;
#	close SCORES;
	$set_ct++;
    }
    return $set_ct;
}

=head2 format_alignments

 Title   : format_alignments
 Usage   : $project->format_alignments();
 Function: Converts cmalign stockholm alignments into fasta format, which is necessary
           for downstream steps. Also formats headers to remove weird cmalign
           characters ("\" and "-" are converted to "_")"
 Example : $project->format_alignments();
 Returns : None
 Args    : None

=cut

sub format_alignments{
    my ( $self ) = shift;
    my @domains   = @{ $self->{"domains"} };
    my $set_ct    = 0;
    foreach my $set ( @domains) {
	my $infile  = $self->{"db"}->{"all_align"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all.sto";
	my $outfile = $self->{"db"}->{"all_align"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all.fa";
	my $alnin   = Bio::AlignIO->new( -file => $infile,    -format => 'stockholm' );
	my $seqout  = Bio::SeqIO->new( -file => ">$outfile", -format => 'fasta' );
	my $alnnum  = 0;
	while( my $aln = $alnin->next_aln() ){
	    if( $alnnum > 0 ){
		warn("there is more than one alignment in $alnin, exiting\n");
		exit(0);
	    }
	    for my $seq( $aln->each_seq() ){
		my $id = $seq->display_id();
		$id =~ s/[\/|\-]/\_/g;
		$seq->display_id( $id );
		$seqout->write_seq( $seq );
	    }
	    $alnnum++;
	}
    }
}

=head2 run_align_qc

 Title   : run_align_qc
 Usage   : $project->run_align_qc();
 Function: Cleans up a fasta formatted multiple sequence alignment by removing 
           highly gapped columns (insert states) and very poorly aligning sequences
           (unlikely true homologs).
 Example : $project->run_align_qc();
 Returns : None
 Args    : None

=cut

#Note: there are other qc scripts included in the repository, the ScoreInternalGap coupled with run_align_ColQC is
#recommended for speed.
sub run_align_qc{
    my ( $self, $lencut ) = @_;
    my @domains   = @{ $self->{"domains"} };
    my $set_ct    = 0;
    warn( "Running alignment qc with length cutoff of $lencut\n" );
    foreach my $set ( @domains) {
        my $inaln    = $self->{"db"}->{"all_align"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all.fa";
        my $outaln   = $self->{"db"}->{"qc_align"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all_qc.fa";
        my $gapaln   = $self->{"db"}->{"qc_align"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all_qc.gap";
        my $cmscores = $self->{"db"}->{"cm_scores"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_cmscores.tab";
        #run with flat on to prevent printing of seq coords in headers                                                                                                    
        my @args     = ("-i $inaln", "-o $outaln", "-s $cmscores", "-flat", "-lc $lencut -g $gapaln");
        my $results  = capture( "perl " . $self->{"workdir"} . "align2profile_qc_ScoreInternalGap.pl @args");
        if( $EXITVAL != 0 ){
            warn("Error running align2profile_qc_ScoreInternalGap.pl for $inaln!\n");
            exit(0);
        }
        $set_ct++;
        warn( "Alignment qc for $set is complete.\n" );
    }
    return $set_ct;
}

=head2 run_align_ColQC

 Title   : run_align_ColQC
 Usage   : $project->run_align_ColQC($mincoverage);
 Function: Cleans up a fasta formatted multiple sequence alignment by removing 
           columns with very low coverage
           Code must be compiled as
           g++ align2profile_qc_Col.cpp -o align2profile_qc_Col
 Example : $project->run_align_qc(2);
 Returns : None
 Args    : Minimum number of residues in a column required for the column to be retained

=cut

sub run_align_ColQC{
  my ( $self, $min ) = @_;
  my @domains   = @{ $self->{"domains"} };
  my $set_ct    = 0;
  foreach my $set ( @domains) {
    my $aln      = $self->{"db"}->{"qc_align"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all_qc.fa";
    my $rowaln   = $self->{"db"}->{"qc_align"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all_qc_rows.fa";
    my $gap      = $self->{"db"}->{"qc_align"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all_qc.gap";
    my $rowgap   = $self->{"db"}->{"qc_align"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all_qc_rows.gap";
    rename( $aln, $rowaln);  # This is no longer the best alignment but we'll save it as a different file name
    rename( $gap, $rowgap);
    my $args     = "$rowaln $rowgap $aln $gap $min > " . $self->{"db"}->{"qc_align"} . "log.tmp 2>&1";
    my $EXITVAL  = system($self->{"workdir"} . "align2profile_qc_Col $args");
#    my $args     = "-i $rowaln -j $rowgap -o $aln -p $gap -flat -m $min > " . $self->{"db"}->{"qc_align"} . "log.tmp 2>&1";
#print "RUNNING: perl " . $self->{"workdir"} . "align2profile_qc_Col.pl $args\n";
#    my $EXITVAL  = system("perl " . $self->{"workdir"} . "align2profile_qc_Col.pl $args");
    system( "cat " . $self->{"db"}->{"qc_align"} . "log.tmp"); #Print stdout and stderr of align2profile_qc_Col.pl to stdout (the logfile)
    if( $EXITVAL != 0 ){
      warn("Error running align2profile_qc_Col.pl for $aln, exited with code $EXITVAL!\n");
      exit(0);
    }
    $set_ct++;
    warn( "Alignment column qc for $set is complete.\n" );
  }
  return $set_ct;
}


=head2 align_to_seqs

 Title   : align_to_seqs
 Usage   : $project->align_to_seqs();
 Function: Converts fasta alignment sequences into non-alignment frame sequences by 
           removing all gap characters ( "-", "." );
 Example : $project->run_align_qc();
 Returns : None
 Args    : None

=cut

sub align_to_seqs{
    my ( $self, $formatid )  = @_;
    my @domains   = @{ $self->{"domains"} };
    my $set_ct    = 0;
    foreach my $set ( @domains) {
	my $inaln     = $self->{"db"}->{"qc_align"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all_qc.fa";
	my $outseqs   = $self->{"db"}->{"qc_seqs"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all_qc_seqs.fa";
	my $seqin     = Bio::SeqIO->new( -file => $inaln, -format => 'fasta' );
	my $seqout    = Bio::SeqIO->new( -file => ">$outseqs", -format => 'fasta' );
	while( my $seq = $seqin->next_seq() ){
	    if( $formatid ){
		my $id       = $seq->display_id();
		$id    =~ s/[\||\/|\\|\.|\-]/_/g;
		$seq->display_id( $id );
	    }
	    my $sequence = $seq->seq();
	    $sequence    =~ s/[\-|\.]//g;
	    $seq->seq( $sequence );
	    $seqout->write_seq( $seq );
	}
    }
    return $self;
}

=head2 run_tree

 Title   : run_tree
 Usage   : $project->run_tree();
 Function: Wrapper for FastTree. The high quality multiple sequence alignment for
           each sample set is processed by FastTree using the -dna and -pseudocounts
           options.
 Example : $project->run_tree();
 Returns : None
 Args    : None

=cut

sub run_tree{
    my ( $self ) = @_;
    my @domains   = @{ $self->{"domains"} };
    my $set_ct    = 0;
    foreach my $set ( @domains) {
	my $inaln   = $self->{"db"}->{"qc_align"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all_qc.fa";
	my $outtree = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo.tree";
	my @args    = ( "-nt", "-pseudo", "$inaln", "> $outtree" );
	my $results = capture( "FastTree @args" );
	if( $EXITVAL != 0 ){
	    warn("Error running align2profile_qc.pl for $inaln!\n");
	    exit(0);
	}
    }    
}

#Not only do we prune out reference sequences, we remove problematic characters from sequence ids

sub prune_tips{
    my $self       = shift;
    my $reformatid = shift;
    my @domains    = @{ $self->{"domains"} };
    foreach my $set( @domains ){
	my $intree    = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo.tree";
	my $refalign  = $self->{"db"}->{"ref_align"} . "SSU_" . $set . "_ref.stk";
	my $outtree   = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned_tmp.tree";
	my $outsmooth = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.tree";
	my @seqids    = ();
	my $align     = Bio::AlignIO->new( -file => $refalign, -format => "stockholm" ); 
	while (my $aln = $align->next_aln() ){
	    for my $seq( $aln->each_seq() ){
		my $id = $seq->display_id();
		push @seqids, $id;
	    }
	}
	$align     = ();
	my $treei  = new Bio::TreeIO( -file => $intree, -format => "newick" );
	my $tree   = $treei->next_tree;
	my @leaves = $tree->get_leaf_nodes;
	foreach my $leaf( @leaves ){
	  my $id = $leaf->id();
	  if( $reformatid ){
	    $id    =~ s/[\||\/|\\|\.|\-]/_/g;
	    $leaf->id( $id );
	  }
	  foreach my $refseq( @seqids ){
	    next if ( $self->{"simulation"} );
	    chomp $_;
#	    if( $id =~ m/$refseq/ ){
	    #To ensure reference tree can be built properly for simulation project
	    if( $id eq $refseq ){
	      $tree->remove_Node($leaf);
	      $tree->contract_linear_paths();
	    }
	  }
	}
	my $treeout = new Bio::TreeIO( -file => ">$outtree", -format => "newick" );
	$treeout->write_tree( $tree );
	#Use Bio::Phylo to remove any internal, unbranched nodes, esp at root
	open( TREEIN,  $outtree)       || die "Can't open pruned, rough tree $outtree for read: $!\n";
	open( TREEOUT, ">$outsmooth" ) || die "Can't open pruned, smooth tree $outsmooth for write: $!\n";
	while (<TREEIN>){
	  chomp;
	  my $tree = parse( -format => 'newick', -string => $_ )->first;
	  my $string1 = $tree->to_newick;
	  my $string2 = $tree->remove_unbranched_internals->to_newick;
	  print TREEOUT "$string2\n";
	}
	close TREEIN;
	close TREEOUT;
      }
    return $self;
}

=head2 tree_to_matrix_R

 Title   : tree_to_matrix_R
 Usage   : $project->tree_to_matrix_R();
 Function: Wrapper for tree_to_matrix.R, which is released with this package and
           converts a newick tree to a phylip formatted distance matrix. Execution
           requires R and the ape R library. 
 Example : $project->tree_to_matrix_R();
 Returns : None
 Args    : None

=cut

sub tree_to_matrix_R{
  #Call an R script with the following:
  my ( $self )  = shift;
  my @domains   = @{ $self->{"domains"} };
  my $set_ct    = 0;
  foreach my $set ( @domains) { 
    my $intree    = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.tree";
    my $outmat    = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.Rmat";
#    my $refnames  = $self->{"db"}->{"ref_align"} . "SSU_" . $set . "_refnames.list";
    my @args = ();
    @args = ( "--slave", "--args", "$intree", "$outmat", "< " . $self->{"workdir"} . "tree_to_matrix.R");
    my $results = capture( "R @args" );
    if( $EXITVAL != 0 ){
      warn("Error running tree_to_matrix.R for $intree!\n");
      exit(0);
    }
  }
}

=head2 tree_to_matrix_cpp

 Title   : tree_to_matrix_cpp
 Usage   : $project->tree_to_matrix_cpp(start,end,cutoff);
 Function: Wrapper for tree_to_matrix.cpp, which is released with this package and
           converts a newick tree to a phylip formatted distance matrix. 
           Compilation requires PhyloTree.h and the Boost c++ library. 
           g++ -I /usr/include/boost/ tree_to_matrix.cpp -o tree_to_matrix
 Example : $project->tree_to_matrix_cpp(0, 100, 0.1);
 Returns : None
 Args    : starting (from 0) and ending rows of the matrix to print to file.
           If start=end=0 the entire matrix is printed
           Print only distances < the cutoff

=cut

sub tree_to_matrix_cpp{
  #Call an c++ script with the following:
  my $self;
  my $mstart = 0;
  my $mend   = 0;
  my $format = 'E';
  my $cutoff = 0.1;
  my $do_pruning = 1;
  my $dom = "";
  ( $self, $mstart, $mend, $format, $cutoff, $do_pruning, $dom )  = @_;
  my @domains;
  if( $dom ){
    @domains = ( $dom );
  } else {
    @domains = @{ $self->{"domains"} };
  }
  my $set_ct    = 0;
  my $tag = "";
  if( $mend != 0 ){
    $tag = "_row".$mstart."to".$mend;
  }
  foreach my $set ( @domains) { 
    my $intree       = $self->{"db"}->{"tree"}   . $self->{"sample"}->{"name"} .        "_SSU_" . $set . "_FT_pseudo.tree";
    my $tmptree      = $self->{"db"}->{"tree"}   . $self->{"sample"}->{"name"} .        "_SSU_" . $set . "_FT_pseudo_pruned_tmp.tree";
    my $prunetree    = $self->{"db"}->{"tree"}   . $self->{"sample"}->{"name"} .        "_SSU_" . $set . "_FT_pseudo_pruned.tree";
    my $refalign     = $self->{"db"}->{"ref_align"} . "SSU_" . $set . "_ref.fa";
    my $logfile      = $self->{"db"}->{"tree"} . "log" . $tag . ".tmp";
    my ( $outmat, $frqfile );
    if( $format eq "M" ){
      # Print for mothur
      $outmat    = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . $tag . "_SSU_" . $set . "_FT_pseudo_pruned.phymat";
      $frqfile   = "";
    } elsif( $format eq "E" ){
      # Print for ESPRIT
      $outmat    = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . $tag . "_SSU_" . $set . "_FT_pseudo_pruned.ndist";
      $frqfile   = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"}        . "_SSU_" . $set . "_FT_pseudo_pruned.frq";
    }
    my $args = "$intree $tmptree $prunetree $refalign $outmat $mstart $mend $format $do_pruning $frqfile $cutoff > $logfile 2>&1";
    print "Running ./tree_to_matrix $args\n";
    my $EXITVAL = system( "./tree_to_matrix $args" );
    system( "cat $logfile"); #Print stdout and stderr of tree_to_matrix to stdout (the logfile)
    if( $EXITVAL != 0 ){
      warn("Error running tree_to_matrix for $intree!\n");
      exit(0);
    }
  }
}


sub format_matrix_to_phylip{
    my $self   = shift;
    my @domains   = @{ $self->{"domains"} };
    my $set_ct    = 0;
    foreach my $set ( @domains) {
	my $inmat  = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.Rmat";
	my $outmat = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.phymat";
	my @args     = ("-i $inmat", "-o $outmat");
	my $results  = capture( "perl " . $self->{"workdir"} . "format_Rphylip.pl @args" );
	if( $EXITVAL != 0 ){
	    warn("Error formatting matrix for $inmat!\n");
	    exit(0);
	}
    }    
    return $self;
}


=head2 run_mothur

 Title   : run_mothur
 Usage   : $project->run_mothur();
 Function: Wrapper for mothur, which leverages a phylip distance matrix to
           cluster sequences into OTUs. Accepts a clustering cutoff threshold
           and method clustering method (average, nearest, furthest). Output
           location will be optionable in next mothur update.
 Example : my $clust_cutoff = 0.05;
           my $clust_method = "average";
           $project->run_mothur( $clust_cutoff, $clust_method );
 Returns : None
 Args    : A scalar of clustering threshold and a scalar of clustering method

=cut

#NEED TO RESOLVE OUTPUT SETTINGS!
sub run_mothur{
  #Run mothur
  my ( $self, $cutoff, $method ) = @_;
  my @domains   = @{ $self->{"domains"} };
  my $set_ct    = 0;
  foreach my $set ( @domains) {    
#    my $inmat    = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo.Rmat";
    my $inmat     = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.phymat";
#    my $seqfile  = $self->{"db"}->{"SSU_reads"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . ".fa";
    my $seqfile  = $self->{"db"}->{"qc_seqs"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all_qc_seqs.fa";
    my $outdir   = $self->{"db"}->{"otudir"} . $set . "/";
    unless( -e $outdir ){
	make_path( $outdir );
    }
    #set.dir will be available in the next MOTHUR release
#    my $command  = "#set.dir(output=$outdir); read.dist(phylip=$inmat, cutoff=$cutoff); cluster(method=$method); bin.seqs(fasta=$seqfile)";
#    my $command  = "#read.dist(phylip=$inmat, cutoff=$cutoff); cluster(method=$method); bin.seqs(fasta=$seqfile)";
    #Having problems getting bin.seqs to work inside mothur, so we'll dot it ourselves later
    my $command  = "#set.dir(output=$outdir); read.dist(phylip=$inmat, cutoff=$cutoff); cluster(method=$method, cutoff=$cutoff);";
    print "executing mother using $command\n";
    my $results  = capture( "mothur \"$command\"" );
    if( $EXITVAL != 0 ){
      warn("Error running mothur for $inmat!\n");
      exit(0);
    }  
  }
}

=head2 run_Ecluster

 Title   : run_Ecluster
 Usage   : $project->run_Ecluster();
 Function: Wrapper for ESPRITs hcluster, which leverages a phylip distance matrix to
           cluster sequences into OTUs.
 Example : $project->run_Ecluster( );
 Returns : None
 Args    : 

=cut

sub run_Ecluster{
  #Run Ecluster
  my ( $self, $cutoff, $method ) = @_;
  my @domains   = @{ $self->{"domains"} };
  my $set_ct    = 0;
  foreach my $set ( @domains) {    
    my $inlist   = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.ndist";
    my $sortlist = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.ndist_sort";
    my $frqfile  = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.frq";
    my $logfile  = $self->{"db"}->{"matrix"} . "log.tmp";
    my $command;
#    if( -e $sortlist ){
#      $command   = "-t 100000 -u         $sortlist $frqfile > $logfile 2>&1";
#    } else {
      $command   = "-t 100000 -b 100000000 $inlist $frqfile > $logfile 2>&1";
#    }
    print "executing ESPRIT's hcluser using $command\n";
    my $EXITVAL  = system( "./hcluster $command" );
    system( "cat $logfile" );
    if( $EXITVAL != 0 ){
      warn("Error running Ecluster for $inlist!\n");
      exit(0);
    }  
  }
}


=head2 split_query

 Title   : split_query
 Usage   : $project->split_query($type $num);
 Function: Split an input sample file of reads into several smaller files
           with the names $SampleFile_$FileNumber
 Example : my @subsample_files =$project->split_query("blast", 5);
           # run blast on each input file
           $project->stitch_blast(@subsample_files); # put together all the blast results
           my @subsample_files =$project->split_query("alignQC", 5);
           # run alignQC on each input file
           $project->stitch_alignQC(@subsample_files); # put together all the alignQC results  
 Returns : Array of the names of the files produced
 Args    : Contents of the original file are split into $num files, 
           with an equal number of reads in each (deliminted with a >)

=cut

sub split_query{

    my ( $self, $type, $nsplit) = @_;
    my @domains   = @{ $self->{"domains"} };
    my @subfilenames;   # file names returned, used to submit new jobs
    my ( $name, $path, $suffix ) = fileparse( $self->{"sample"}->{"path"}, qr/\.[^.]*/ );
    my @list = @domains;
    if( $type =~ "blast"){
      # Since I don't use $set for the blast stuff I only want to run everything once
      @list = $list[1];
    }

    foreach my $set ( @list ) {
      $path =~ s/raw\///g;
      my $oldraw    = $path    . "raw/"           . $name    . ".fa";
      my $oldscore  = $path    . "aligns/scores/" . $name    . "_SSU_" . $set . "_cmscores.tab";
      my $infile;
      if( $type =~ "blast"){
	$infile = $self->{"sample"}->{"path"};
      } elsif($type =~ "alignQC") {
	$infile = $self->{"db"}->{"all_align"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all.fa";
      } else {
	print ("Don't recognise type $type, don't know what file to split\n");
	exit(0);
      }

      my $nlines   = capture( "grep -c \">\" \"$infile\"");
      my $nsubfile = ceil( $nlines / $nsplit );
      open( FULLFILE, $infile ) || die "Can't open sample file ".$infile.": $!\n";
      my $lcount =$nlines;
      my $fcount =1;
      print "splitting ".$infile." into ".$nsplit." files with ".$nsubfile." reads each\n";

      while( <FULLFILE> ){
	# Read in the next line of the file
	if( $_ =~ m/>/ ){
	  # This is a new read, so I finished the old one, so increment the # lines
	  $lcount++;
	  
	  # Check if this read needs to go in a new file
	  if( ($lcount>=$nsubfile) && ($fcount<=$nsplit) ){
	    # The last file will have a few extra from the rounding down
	    if( $fcount >= 1 ){ close( SUBFILE ); }
	    # Set-up the new directories and files
	    my $newname = $name . "_" . $fcount;
	    my $newpath = $path;
	    $newpath    =~ s/$name/$newname/g;
	    my $rawfile = $newpath . "raw/" . $newname . ".fa";
	    make_path( $newpath."raw/", $newpath."aligns/raw/", $newpath."aligns/scores/");
	    my $outfile;
	    if( $type =~ "blast"){
	      $outfile = $rawfile;
	    } elsif($type =~ "alignQC") {
	      my $alnfile   = $newpath . "aligns/raw/"    . $newname . "_SSU_" . $set . "_all.fa";
	      my $scorefile = $newpath . "aligns/scores/" . $newname . "_SSU_" . $set . "_cmscores.tab";
	      unlink($scorefile);
	      symlink($oldscore, $scorefile);
	      unlink($rawfile);
	      symlink($oldraw,   $rawfile);
	      $outfile = $alnfile;
	    }
	    open( SUBFILE, ">".$outfile) || die "Can't open new sub-sample file: $!\n";
	    if( $set eq $domains[-1] ){
	      # Only save the files to submit on one version of $set, the parallel code will run BAC and ARC
	      push(@subfilenames, $rawfile);  # Return the raw file names, this is needed to submit
	    }
	    # Reset the line count, increment the file number
	    $lcount = 0;
	    $fcount++;
	  }
	}
	# Print this line to the file
	print SUBFILE $_;
      }
      close( FULLFILE );
      close( SUBFILE );
    }
    return @subfilenames;
}

=head2 split_tree

 Title   : split_tree
 Usage   : $project->split_tree($num);
 Function: Figure out how to split up the tree_to_matrix printing
 Example : my @subsample_files =$project->split_tree(5)
           # run tree_to_matrix each part of the tree
           $project->stitch_matrix(@subsample_files); # put together all the matrices
 Returns : Array of the commands to submit
 Args    : The tree is split into $num parts, with an equal number of leaves in each

=cut

sub split_tree{

  my ( $self, $num ) = @_;
  my @domains   = @{ $self->{"domains"} };
  my @subfilenames;        # commands to run
  foreach my $set ( @domains ){
    my $infile  = $self->{"sample"}->{"path"};
    my $qcaln   = $self->{"db"}->{"qc_align"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all_qc.fa";
    my $refaln  = $self->{"db"}->{"ref_align"} . "SSU_" . $set . "_ref.fa";
    my $nleaves = capture( "grep -c \">\" \"$qcaln\"") - capture( "grep -c \">\" \"$refaln\"");
    my $step    = ceil( $nleaves/$num );
    my $start   = 0;
    my $end     = $step;
    print "splitting ".$infile." tree_to_matrix printing into ".$num." jobs with ".$step." rows each\n";
    while( $start < $nleaves ){
      my $command = $infile . " -s $start -e $end -d $set";
      push(@subfilenames, $command);
      $start = $end   + 1;
      $end   = $start + $step;
    }
  }
  return @subfilenames;

}


=head2 stitch_blast

 Title   : stitch_blast
 Usage   : $project->stitch_blast(@list_of_files);
 Function: Concatenate several blast files together
           Files must have the filename $dir/$base_$num.$suffix
 Example : my @subsample_files = $project->split_query("blast",5);
           # run blast on each input file
           $project->stitch_blast(@subsample_files); # put together all the blast results
 Returns : 
 Args    : 

=cut
    
sub stitch_blast{
 
  my ( $self, @ffiles ) = @_;
  my $blastpath = $self->{"db"}->{"blastout"};
  my $SSUpath   = $self->{"db"}->{"SSU_reads"};
  my $sampdir   = $blastpath;
  $sampdir  =~ s/blastout/samples/g;
  
  foreach my $set( @{ $self->{"domains"} } ){
    my $blastsuffix  = "_SSU_" . $set . ".blast";
    my $SSUsuffix    = "_SSU_" . $set . ".fa";
    
    # I will cat all the blastout and SSU files together into these files
    my $blastoutfile = $blastpath . $self->{"sample"}->{"name"} . $blastsuffix;
    my $SSUoutfile   = $SSUpath   . $self->{"sample"}->{"name"} . $SSUsuffix;
    # Delete them to be sure we start with a clean slate, since I will append these files
    unlink($blastoutfile);
    unlink($SSUoutfile);
    
    foreach my $ffile ( @ffiles ){
      my ( $name, $path, $suffix ) = fileparse( $ffile, qr/\.[^.]*/ );
      my $subsampdir = $sampdir . $name;                       # samples/SubSampleName
	    my $bfile = $blastpath          . $name . $blastsuffix;  # The subsample's blastout file
      my $sfile = $subsampdir."/SSU/" . $name . $SSUsuffix;    # The subsample's SSU file
      
      # Put all the blast results into one file
      # Put all the SSU files together
      capture( "cat $bfile >> $blastoutfile" );
      capture( "cat $sfile >> $SSUoutfile"   );
      
      # Cleanup
      unlink($ffile);                  # Delete the subsample fasta file in the original directory
      unlink($bfile);                  # Delete the blastout file
      my @domains = @{ $self->{"domains"} };
      if( $set eq $domains[-1] ){
	capture( "rm -rf $subsampdir");  # Delete the new sample directory that was made (this includes $sfile)     
      }
    }
  }
}

=head2 stitch_alignQC

 Title   : stitch_alignQC
 Usage   : $project->stitch_alignQC(@list_of_files);
 Function: Concatenate several alignQC files together
           Files must have the filename $dir/$base_$num_SSU_BAC_all_qc.fa
 Example : my @subsample_files = $project->split_query("alignQC", 5);
           # run alignQC on each input file
           $project->stitch_alignQC(@subsample_files); # put together all the alignQC results
 Returns : 
 Args    : 

=cut
    
sub stitch_alignQC{
 
  my ( $self, @ffiles ) = @_;
  
  foreach my $set( @{ $self->{"domains"} } ){
    my $alnsuffix = "_SSU_" . $set . "_all_qc.fa";    
    my $gapsuffix = "_SSU_" . $set . "_all_qc.gap";    
    my $FullalignQC    = $self->{"db"}->{"qc_align"} . $self->{"sample"}->{"name"} . $alnsuffix;
    my $FullalignQCgap = $self->{"db"}->{"qc_align"} . $self->{"sample"}->{"name"} . $gapsuffix;
    # Delete the output files to be sure we start with a clean slate
    unlink($FullalignQC);
    unlink($FullalignQCgap);
    my $numseq = 0;
    my @coverage = ();

    foreach my $ffile ( @ffiles ){
      my ($samplenm, $subsampdir, $suffix) = fileparse( $ffile, qr/\.[^.]*/ );
      $subsampdir =~ s/raw\///;
      my $alnfile =  $ffile;
      $alnfile    =~ s/raw/aligns\/qc/;
      $alnfile    =~ s/$suffix/$alnsuffix/;
      my $gapfile =  $alnfile;
      $gapfile    =~ s/$alnsuffix/$gapsuffix/;
      # cat all the alignment lines together
      capture( "cat $alnfile >> $FullalignQC" );
      # Keep a running total of the average files
      open( GAP, "$gapfile" ) || die "can't open output $gapfile in stitch_alignQC:$!\n";
      my $seqs = <GAP>;
      $numseq  = $numseq + $seqs;
      my $cov_str = <GAP>;
      my @cov = split(" ",$cov_str);
      if( @coverage ){
	my $nres = @coverage;
	for( my $i=0; $i<$nres; $i++ ){
	  $coverage[$i] = $coverage[$i] + $cov[$i]; # Accumulate the coverage
	}
      } else {
	@coverage = @cov;
      }
      close GAP;

      my @domains = @{ $self->{"domains"} };
      if( $set eq $domains[-1] ){
	# Delete the new sample directory that was made
	capture( "rm -rf $subsampdir");  
      }
    }
    open( GAP, ">>$FullalignQCgap" ) || die "can't open output $FullalignQCgap in stitch_alignQC:$!\n";
    # Print the total number of sequences
    print GAP "$numseq\n";
    # Print the coverage
    my $format = ("%1.0f " x @coverage)."\n";
    printf GAP $format, @coverage;
    # Print the % coverage
    my $nres = @coverage;
    for(my $i=0; $i< $nres; $i++ ){
      $coverage[$i] = $coverage[$i]/$numseq;
    }  
    $format = ("%1.5f " x @coverage)."\n";
    printf GAP $format, @coverage;
    close GAP;
  }
}

=head2 stitch_matrix

 Title   : stitch_matrix
 Usage   : $project->stitch_matrix(@list_of_commands);
 Function: Concatenate several matrix files together
           Files must have the filename $dir/matrix/$base_row$Xto$Y_SSU_BAC_FT_pseudo_pruned.phymat
 Example : my @subsample_files = $project->split_tree(5);
           # run tree_to_matrix on each part of the tree
           $project->stitch_matrix(@subsample_files); # put together all the matrices
 Returns :
 Args    : list of command with start/end rows indicated

=cut

sub stitch_matrix{

  my ( $self, @ffiles ) = @_;
  opendir(MATDIR, $self->{"db"}->{"matrix"});
  my @dircontent = readdir(MATDIR);  
  closedir(MATDIR);
  foreach my $format ( ("ndist", "phymat") ){
    foreach my $set( @{ $self->{"domains"} } ){
      my $base       = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"};
      my $suffix     = "_SSU_" . $set . "_FT_pseudo_pruned." . $format;
      my $Fullmatrix = $base . $suffix;
      my @files = grep /$suffix$/, @dircontent;
      print scalar(@files) . " files with suffix $suffix in directory" . $self->{"db"}->{"matrix"} . "\n";
      if( scalar(@files) > 1 ){
	# Only cat togheter a list of files (otherwise this is probably the wrong format)
	# Delete the output files to be sure we start with a clean slate
	unlink($Fullmatrix);
	foreach my $command ( @ffiles ){
	  my( $infile, $s, $start, $e, $end, $d, $dom ) = split(' ', $command);
	  if( $dom =~ $set ){
	    my $tag     = "_row".$start."to".$end;
	    my $subfile = $base . $tag . $suffix;
	    # cat all the alignment lines together
	    capture( "cat $subfile >> $Fullmatrix" );
	    unlink($subfile);
	  }
	}
      }
    }
  }

}


######################################
# SIMULATION SPECIFIC FUNCTIONS HERE #
######################################

#if called, simulation switch is set on
sub set_simulation{
    my $self = shift;
    $self->{"simulation"} = 1;
    return $self;
}

sub build_read2source_tab{
    my $self   = shift;
    my $raw    = $self->{"db"}->{"reads"} . $self->{"sample"}->{"file"};
     my $tab    = $self->{"db"}->{"reads"} . "sim-read2source-lookup.tab";
    my %lookup = ();
    my $seqs = Bio::SeqIO->new( -file => $raw, -format => 'fasta' );
    while( my $seq = $seqs->next_seq() ){
	my $id   = $seq->display_id();
	my $desc = $seq->description();
        #metasim description has: SOURCE_1="NC_008321/1-1970" (
	if( $desc =~ m/SOURCE\_1\=\"(.*)\/\d+\-\d+\"/){
	    my $source = $1;
	    $lookup{$id} = $source;
	}
	else{
	    print ("Can't grab source for read $id in build_read2source_tab!\n");
	    exit(0);
	}
    }
    open( TAB, ">$tab" ) || die "can't open output $tab in build_read2source_tab:$!\n";
    foreach my $id( sort( keys( %lookup ) ) ){
	my $source = $lookup{$id};
	print TAB join( "\t", $id, $source, "\n");
    }
    close TAB;
    return $self;
}

#We're going to remove sequences from the qc alignment file that are full
#length reference or source versions of any read in our raw seq file.
#Since we want to leverage downstream functions in a consistant mannar,
#we'll backup the unedited alignment file and dump our new alignment
#into the old spot.

#use formatid to remove the trailing end of each sequence id (coordinate-coordinate)

#also, last step is replacing the sim read id with its corresponding source id. thus
#each source (or reference) is represented exactly one time in the resulting output

#note that now we're running all reads through pipeline and pruning
#down to 70 at this step. Need to do the following:
#1. open alignment file, get all seq headers
#2. grab all read headers.
#3. use read2source lookup tab to map reads into source bins
#4. for each source, randomly select one read
#5. across all sources, randomly select 70 to keep in alignment
#6. if sim 1, prune back full length source seqs for these 70 
#7. dump all passing seqs to output.

sub prune_alignment{
    my $self      = shift;
    my $tab       = $self->{"db"}->{"reads"} . "sim-read2source-lookup.tab";
    my $simtype   = shift;
    my $formatid  = 1; #shift;
    my $n_sources = 50;
    my @domains   = @{ $self->{"domains"} };
    #Build read to source and source to read lookup tables (redundant but easy)
    my %rlookup    = ();
    my %slookup   = ();
    open( TAB, $tab ) || die "can't open $tab for read at prune_alignment:$!\n";
    while( <TAB> ){
	chomp $_;
	my ( $readid, $source ) = split( "\t" , $_ );
	$rlookup{$readid}       = $source;
	my @reads = ();
	if( defined( $slookup{$source} ) ){
	  @reads = @{ $slookup{$source} };
	  push( @reads, $readid );
	  $slookup{$source} = \@reads;
	}
	else{
	  push @reads, $readid;
	  $slookup{$source} = \@reads;
	}
    }
    close TAB;
    #Now let's move through the alignment, select reads to retain and punt source sequences
    my $set_ct    = 0;
    foreach my $set ( @domains) {
	#backup qc alignment file
	my $align       = $self->{"db"}->{"qc_align"} . $self->{"sample"}->{"name"} . 
	    "_SSU_" . $set . "_all_qc.fa";
	#for sim3
	my $alignout;
	if( $simtype == 3 ){
	  $alignout       = $self->{"db"}->{"qc_align"} . $self->{"sample"}->{"name"} . 
	    "_SSU_" . $set . "_all_qc_sim3.fa";
	}
	my $alnbk       = $self->{"db"}->{"qc_align"} . $self->{"sample"}->{"name"} . 
	    "_SSU_" . $set . "_all_qc.backup.fa";
	#we want to capture which seq ids in our final, pruned alignment are the sources
	#used to build the simulated reads. we made a typo earlier, so these next two 
	#blocks amend the mistake. can be removed next pass through the code
	my $source_list = $self->{"db"}->{"qc_align"} . $self->{"sample"}->{"name"} . 
	  "_SSU_" . $set . "source_list.txt"; 
	if( -e $source_list ){
	  unlink $source_list;
	}
	$source_list = $self->{"db"}->{"qc_align"} . $self->{"sample"}->{"name"} . 
	  "_SSU_" . $set . "_source_list.txt";        
	open( SOURCE_LIST , ">$source_list" ) || 
	  die "Can't open source list $source_list for write in prune alignments:$!\n";
	unless( $simtype ==3 ){
	  if( -e $alnbk ){
	    if( 1 ){
	      unlink( $alnbk );
	    }  
	    else{
	      print "Um, I hate to be a bother, but I don't want to overwrite a previously backed up file. You'd better get in here and check this out:\n$alnbk\n";
	      exit(0);
	    }
	  }
	}
	my $alnin;
	if( $simtype != 3){
	  move( $align, $alnbk );
	  $alnin  = Bio::AlignIO->new( -file => $alnbk, -format => 'fasta' );
	}
	else{
	  $alnin = Bio::AlignIO->new( -file => $align, -format => 'fasta' );
	}
	#Will print alnout to seqio object - handles seq modifications a bit better
	#and since it's fasta, this isn't a problem
	my $alnout;
	if( $simtype == 3){
	  $alnout = Bio::SeqIO->new( -file => ">$alignout", -format => 'fasta' );
	}
	else{
	  $alnout = Bio::SeqIO->new( -file => ">$align", -format => 'fasta' );
	}
	while( my $aln = $alnin->next_aln() ){
	    my @seqs_to_pass = ();
	    if( $formatid ){
		$aln->set_displayname_flat();
	    }
	    #first loop, get all sequence ids. Use id_status to determine if
	    #seq will be retained or dropped. Note this might mean different
	    #things depending on the simulation type
	    my %id_status = ();
	    #initialize each read seq setting as 1, meaning retain
	    my $source_ct = 0;
	    for my $seq( $aln->each_seq() ){
	      my $id = $seq->display_id();
	      $id_status{$id} = 1;	  
	      if( $id =~ m/^r/ ){
		next;
	      }
	      else{
		$source_ct++;
	      }
	    }
	    #check each source's sim read and determine if it passed qc. if so, 
	    #add to the sample box
	    foreach my $source( sort( keys( %slookup ) ) ){
	      next unless( defined( $slookup{$source} ) );
	      my @box = ();
	      my @source_reads = @{ $slookup{$source} };
	      foreach my $read( @source_reads ){ 
		if( exists( $id_status{$read} ) ){
		  push @box, $read;
		}
	      }
	      print "source reads is:\n";
	      print Dumper @source_reads;
	      #now draw at random a single passing read from the source's box
	      #unless simtype 3, where we draw 2
	      if( $simtype == 3 ){
		print "checking $source which has the following reads\n";
		print join( "\t", @box, "\n");
		my @retained_reads = ();
		if( scalar( @box ) == 1 ){
		  print "not enough items in our source's read box for $source. Reads for this " .
		    "are " . join(' ', @{ $slookup{$source} } ) . "\n";
		  next;		  
		}
		elsif( scalar(@box) < 1 ){
		  print "no items in our source's read box for $source. Reads for this " .
		    "are " . @{ $slookup{$source} } . "\n";
		  next;
		}
		else{
		  @retained_reads = @{ get_random_subset( 2, \@box ) };
		  foreach my $item( @retained_reads ){
		    print "drew $item\n";
		  }
		}
		#set reads that aren't drawn to status 0, meaning pass
		foreach my $read ( @box ){
		  my $hashit = 0;
		  foreach my $retained_read( @retained_reads ){
		    if( $read eq $retained_read ){
		      $hashit = 1;
		    }
		  }
		  if( $hashit == 0){
		    $id_status{$read} = 0;
		  }
		  print "considering $read\t";
		  print $id_status{$read} . "\n";		 
	        }
	      }
	      else{
		my $retained_read;
		if( scalar( @box ) == 1 ){
		  $retained_read = $box[0];
		}
		elsif( scalar(@box) < 1 ){
		  print "no items in our source's read box for $source. Reads for this " .
		    "are " . @{ $slookup{$source} } . "\n";
		  next;
		}
		else{
		  $retained_read = @{ get_random_subset( 1, \@box ) }[0];
		}
		#set reads that aren't drawn to status 0, meaning pass
		foreach my $read ( @box ){
		  unless( $read eq $retained_read ){
		    $id_status{$read} = 0;
		  }
		}
	      }
	    }
	    #now that all sources have been processed to exhibit a single, retained read,
	    #sample a subset of n_sources at random to represent in the dataset
	    #first check to make sure there are n_source number of unique, retained, reads
	    my $npass_reads = 0;
	    my @passing_read_ids = ();
	    foreach my $id( sort( keys( %id_status ) ) ){
	      next unless ( $id =~ m/^r/ );
	      if( $id_status{$id} == 1 ){
		$npass_reads++;
		push @passing_read_ids, $id;
	      }
	    }
	    if( $simtype == 3 ){
	      if( $npass_reads < ( 2 * $n_sources ) ){
		warn( "Not enough uniquly sourced reads to continue when using source " .
		      "sampling size of $n_sources. Max available is $npass_reads!\n");
		exit(0);
	      }
	    }
	    else{
	      if( $npass_reads < $n_sources ){
		warn( "Not enough uniquly sourced reads to continue when using source " .
		      "sampling size of $n_sources. Max available is $npass_reads!\n");
		exit(0);
	      }
	    }
	    #now grab all of the passing read source ids and sample randomly
	    my @source_box = ();
	    my %source_proc = (); #hack to deal with simtype 3
	    if( $simtype == 3) {
	      foreach my $readid( @passing_read_ids ){
		my $source = $rlookup{$readid};	      
		$source_proc{$source}++;
		if( $source_proc{$source} == 2){
		  push @source_box, $source;	  
		}
	      }
	    }
	    else{
	      foreach my $readid( @passing_read_ids ){
		my $source = $rlookup{$readid};	      
		push @source_box, $source;
	      }
	    }
	    print "source box is size " . scalar(@source_box) ."\n";
	    print "sources in box are " . join( ' ', sort(@source_box), "\n");
	    my @keep_sources = ();
	    if( scalar( @source_box ) == $n_sources ){
	      @keep_sources = @source_box;
	    }
	    elsif( scalar( @source_box ) < $n_sources ){
	      warn("Not enough unique simulated reads to proceed!\n");
	      exit(0);
	    }
	    else{
	      @keep_sources = @{ get_random_subset( $n_sources, \@source_box ) };
	    }
	    #now set all passing_read_ids without a source in keep_sources to pass
	    READ: foreach my $readid( @passing_read_ids ){
	      my $source = $rlookup{$readid};
	      foreach my $kept ( @keep_sources ){
		next READ if ($kept eq $source);
	      }
	      $id_status{$readid} = 0;
	    }
	    #now set all sources with retained reads to fail (depending on simtype)
	    #now set that source's status to 0, meaning pass
	    #at least for sim 1, might need to alter for sim 2
	    if( $simtype == 1 || $simtype == 3 ){
	      foreach my $source ( @keep_sources ){
		$id_status{$source} = 0;
	      }
	    }
	    else{
	      warn( "Since this is simulation type 2, we'll retain the source sequences..right?\n");
	    }
	    #at this point, we should have a number of retained sequences
	    #equal to the total number of source sequences in original alignment,
	    #with n_sources of those retained being reads. Let's check:
	    my $nreads = 0;
	    my $ntotal = 0;
	    my %source_draws = (); #this is for sim 3 minimum draw checking
	    foreach my $id( sort( keys( %id_status ) ) ) {
	      next if ( $id_status{$id} == 0 );
	      $nreads++ if( $id =~ m/^r/ );	      
	      $ntotal++;
	      if( $id =~ m/^r/ ){
		my $src = $rlookup{$id};
		$source_draws{$src}++;
		print "$src $id " . $source_draws{$src} . "\n";
	      }
	    }
	    if( $simtype == 3){
	      my $src_pass_num = scalar(@source_box);
	      foreach my $src( sort( keys( %source_draws ) ) ){
		my $read_count = $source_draws{$src};
		if( $read_count != 2 ){
		  warn( "Error drawing reads for $src, only got " . 
			$source_draws{$src} . "\n");
		  $src_pass_num = $src_pass_num - 1;
		}
	      }
	      print "$src_pass_num\n";
	    }
	    if( ($simtype == 1 &&
		 ( $nreads != $n_sources ||
		   $ntotal != $source_ct )) ||
	        ($simtype == 2 && 
		 ( $nreads != $n_sources ||
		   $ntotal != $source_ct + $nreads )) ){
	      warn( "Error in setting alignment sequence retention status:\n" .
		    "Original source count: $source_ct\n" .
		    "Total retained sequences: $ntotal\n" .
		    "Number passing reads: $nreads\n" .
		    "Number of desired reads: $n_sources\n" );
	      exit(0);
	    }
	    print( "Alignment pruning statistics:\n" .
		   "Original source count: $source_ct\n" .
		   "Total retained sequences: $ntotal\n" .
		   "Number passing reads: $nreads\n" .
		   "Number of desired reads: $n_sources\n" );
	    
	    #second loop, do our id formatting and seq punting
	    SEQ: for my $seq ( $aln->each_seq() ){
		my $id = $seq->display_id(); 
		next if ( $id_status{$id} == 0 );
		#Find read ids, replace with source id, print
		if( $id =~ m/^r/ ){
		  if( $simtype == 2 ){
		    print SOURCE_LIST $rlookup{$id} . "\n";
		    $id = "r" . $rlookup{$id};
		    $seq->display_id( $id );
		    $alnout->write_seq( $seq );
		    next SEQ;
		  }
		  if( $simtype == 3 ){
		    print SOURCE_LIST $rlookup{$id} . "\n";
		    my $src    = $rlookup{$id};
		    my $src_ct = $source_proc{$src};
		    my $newid  = $src . "_" . $src_ct;
		    $source_proc{$src}--;
		    $seq->display_id( $newid );
		    $alnout->write_seq( $seq );
		    next SEQ;
		  }
		  else{
		    print SOURCE_LIST $rlookup{$id} . "\n";
		    $seq->display_id( $rlookup{$id} );
		    $alnout->write_seq( $seq );
		    next SEQ;
		  }
		}
		#print source ids that don't match read id
		$alnout->write_seq( $seq );
	      }
	}
      }
   return $self;
}

#Here we're going to cut our tree down to versions with fewer tips. the breakdown follows:
#Sim 1 treetype 1 (s1t1): Don't prune any tips on either source or sim tree
#     - this is simply the input files for sim and ref trees
#s1t2: Prune tree and reference tree down to just source/sim tips
#s2t1: Prune tree down to just source and sim tips. Ref tree of just source tips

sub sim_prune_tips{
    my $self       = shift;
    my $tab        = $self->{"db"}->{"reads"} . "sim-read2source-lookup.tab";
    my $simtype    = shift;
    my $reftree    = shift; #path to the reference tree
    my $reformatid = shift;
    my @domains    = @{ $self->{"domains"} };
    #this may be obsolete given the revamp to the sim workflow
    my @source_to_save = ();
    open( TAB, "$tab" ) || 
      die "can't open sim-read2source-lookup.tab $tab:$!\n";
    while(<TAB>){
      chomp $_;
      my ( $read, $source ) = split( "\t", $_ );
      push @source_to_save, $source;
    }
    close TAB;
    #let's process a domain's data
    foreach my $set( @domains ){
	my $intree      = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo.tree";
	my $refalign    = $self->{"db"}->{"ref_align"} . "SSU_" . $set . "_ref.stk";
	my $outtree     = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned_tmp.tree";
	my $outsmooth   = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.tree";
	my $source_list = $self->{"db"}->{"qc_align"} . $self->{"sample"}->{"name"} . 
	  "_SSU_" . $set . "_source_list.txt";        
	my $routtree    = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_REF_pruned_tmp.tree";
	my $routsmooth  = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_REF_pruned.tree";
	#build a reference table that indicates which ids refer to source/sims
	my %sourceids = ();
	open( SOURCE_LIST, $source_list ) || 
	  die "can't open source list $source_list in sim_prune_tips: $!\n";
	while(<SOURCE_LIST>){
	  chomp $_;
	  $sourceids{$_} = 1;
	}

	#get sim and ref trees and start processing them
	my $treei  = new Bio::TreeIO( -file => $intree, -format => "newick" );
	my $tree   = $treei->next_tree;
	my $rtreei = new Bio::TreeIO( -file => $reftree, -format => "newick" );
	my $rtree  = $rtreei->next_tree;
	#start by processing the simulation tree
	#we put code in this if block to properly scope recycled vars (below)
	if( $simtype == 1 || $simtype == 2 ){
	  my @leaves = $tree->get_leaf_nodes;
	  foreach my $leaf( @leaves ){
	    my $id = $leaf->id();
	    if( $reformatid ){
	      $id    =~ s/[\||\/|\\|\.|\-]/_/g;
	      $leaf->id( $id );
	    }
	    #now we prune
	    unless( defined( $sourceids{$id} ) || $id =~ m/^r/ ){
	      $tree->remove_Node($leaf);
	      $tree->contract_linear_paths();
	    }
	  }
	  my $treeout = new Bio::TreeIO( -file => ">$outtree", -format => "newick" );
	  $treeout->write_tree( $tree );
	}
	#now do the same thing for the reftree if simtype is 1
	if( $simtype == 1 || $simtype == 2){
	  my @leaves = $rtree->get_leaf_nodes;
	  foreach my $leaf( @leaves ){
	    my $id = $leaf->id();
	    if( $reformatid ){
	      $id    =~ s/[\||\/|\\|\.|\-]/_/g;
	      $leaf->id( $id );
	    }
	    unless( defined( $sourceids{$id} ) ){ 
	      $rtree->remove_Node($leaf);
	      $rtree->contract_linear_paths();
	    }
	  }	  
	  my $rtreeout = new Bio::TreeIO( -file => ">$routtree", -format => "newick" );
	  $rtreeout->write_tree( $rtree );
	}
    
	#Use Bio::Phylo to remove any internal, unbranched nodes, esp at root
	#start with sim tree
	if( $simtype == 1 || $simtype == 2 ){
	  open( TREEIN,  $outtree)       || die "Can't open pruned, rough tree $outtree for read: $!\n";
	  open( TREEOUT, ">$outsmooth" ) || die "Can't open pruned, smooth tree $outsmooth for write: $!\n";
	  while (<TREEIN>){
	    chomp;
	    my $tree = parse( -format => 'newick', -string => $_ )->first;
	    my $string1 = $tree->to_newick;
	    my $string2 = $tree->remove_unbranched_internals->to_newick;
	    print TREEOUT "$string2\n";
	  }
	  close TREEOUT;
	}
	#now do the ref tree, if simtype 1
	if( $simtype == 1 ){
	  open( TREEIN,  $routtree)       || die "Can't open pruned, rough tree $routtree for read: $!\n";
	  open( TREEOUT, ">$routsmooth" ) || die "Can't open pruned, smooth tree $routsmooth for write: $!\n";
	  while (<TREEIN>){
	    chomp;
	    my $tree = parse( -format => 'newick', -string => $_ )->first;
	    my $string1 = $tree->to_newick;
	    my $string2 = $tree->remove_unbranched_internals->to_newick;
	    print TREEOUT "$string2\n";
	  }
	  close TREEOUT;
	}
      }
    return $self;
}

#This version interfaces with the c++ version of the tree pruning algorithm to prune a tree and produce a 
#distance matrix (phylip). See tree_to_matrix_cpp for more information on the specifics.
#
#Here we're going to cut our tree down to versions with fewer tips. the breakdown follows:
#Sim 1 treetype 1 (s1t1): Don't prune any tips on either source or sim tree
#     - this is simply the input files for sim and ref trees
#s1t2: Prune tree and reference tree down to just source/sim tips
#s2t1: Prune tree down to just source and sim tips. Ref tree of just source tips

sub sim_tree_to_matrix_cpp{
    #initialized params from tree_to_matrix_cpp
    my $mstart = 0;
    my $mend   = 0;
    my $format = 'M';
    my $cutoff = 0.1;
    my $do_pruning = 1;
    my $tag = "";
    #grab run-time params from method call
    my ( $self, $simtype, $reftree );
    ( $self, $simtype, $reftree, $mstart, $mend, $format, $cutoff, $do_pruning )  = @_;
    my $tab        = $self->{"db"}->{"reads"} . "sim-read2source-lookup.tab";
    my @domains    = @{ $self->{"domains"} };
    my $set_ct     = 0;
    if( $mend != 0 ){
      $tag = "_row".$mstart."to".$mend;
    }
    #this may be obsolete given the revamp to the sim workflow
    my @source_to_save = ();
    open( TAB, "$tab" ) || 
      die "can't open sim-read2source-lookup.tab $tab:$!\n";
    while(<TAB>){
      chomp $_;
      my ( $read, $source ) = split( "\t", $_ );
      push @source_to_save, $source;
    }
    close TAB;
    #let's process a domain's data
    foreach my $set( @domains ){
	my $intree      = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo.tree";
	my $refalign    = $self->{"db"}->{"ref_align"} . "SSU_" . $set . "_ref.fa";
	my $source_list = $self->{"db"}->{"qc_align"} . $self->{"sample"}->{"name"} . 
	  "_SSU_" . $set . "_source_list.txt";        
	my $subref_align = $self->{"db"}->{"qc_align"} . $self->{"sample"}->{"name"} . 
	  "_SSU_" . $set . "subref_align.fa";        
	#simulation tree file settings
	my $tmptree     = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . 
	  "_SSU_" . $set . "_FT_pseudo_pruned_tmp.tree";
	my $prunetree   = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . 
	  "_SSU_" . $set . "_FT_pseudo_pruned.tree";
	my $outmat      = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . 
	  $tag . "_SSU_" . $set . "_FT_pseudo_pruned.phymat";
	#source (formerly called reference) file settings. Not to be confused with the reference db sequences!
	my $rtmptree    = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . 
	  "_SSU_" . $set . "_FT_pseudo_REF_pruned_tmp.tree";
	my $rprunetree  = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . 
	  "_SSU_" . $set . "_FT_pseudo_REF_pruned.tree";
	my $routmat      = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . 
	  $tag . "_SSU_" . $set . "_FT_pseudo_REF_pruned.phymat";
	my $logfile     = $self->{"db"}->{"tree"} . "log" . $tag . ".tmp";
	my $rlogfile     = $self->{"db"}->{"tree"} . "REF_log" . $tag . ".tmp";
	#Come back and add ESPRIT format functionality at a later date, Skip for now
	my $frqfile = "";
        #build a reference table that indicates which ids refer to source/sims
	my %sourceids = ();
	open( SOURCE_LIST, $source_list ) || 
	  die "can't open source list $source_list in sim_prune_tips: $!\n";
	while(<SOURCE_LIST>){
	  chomp $_;
	  $sourceids{$_} = 1;
	}
	#Going to open the reference alignnment and print out those sequences that are not sources 
	#into a sub-reference alignment. This is important given the way the the tree_to_matrix_cpp
	#method prunes tips: sequences in this subref align will be pruned from both trees.
	my $refseqin  = Bio::SeqIO->new( -file => $refalign, -format => 'fasta' );
	my $subrefout = Bio::SeqIO->new( -file => ">$subref_align", -format => 'fasta' );
	while( my $refseq = $refseqin->next_seq() ){
	  my $id = $refseq->display_id();
	  $id =~ s/\/.*$//g;
	  next if ( $sourceids{$id} );
	  $subrefout->write_seq( $refseq );
	}
	#get sim and ref trees and start processing them using tree_to_matrix
	#start with the sim tree
	my $args = "$intree $tmptree $prunetree $subref_align $outmat $mstart $mend $format $do_pruning $frqfile $cutoff > $logfile 2>&1";
	print "Running ./tree_to_matrix $args\n";
	my $EXITVAL = system( "./tree_to_matrix $args" );
	system( "cat $logfile"); #Print stdout and stderr of tree_to_matrix to stdout (the logfile)                 
	if( $EXITVAL != 0 ){
	  warn("Error running tree_to_matrix for $intree!\n");
	  exit(0);
	}
	#now do the source tree
	my $rargs = "$reftree $rtmptree $rprunetree $subref_align $routmat $mstart $mend $format $do_pruning $frqfile $cutoff > $logfile 2>&1";
	print "Running ./tree_to_matrix $args\n";
	$EXITVAL = system( "./tree_to_matrix $rargs" );
	system( "cat $rlogfile"); #Print stdout and stderr of tree_to_matrix to stdout (the logfile)                 
	if( $EXITVAL != 0 ){
	  warn("Error running tree_to_matrix for $reftree!\n");
	  exit(0);
	}		
    }
    return $self;
}


sub sim_tree_to_matrix{
  #Call an R script with the following:
  my ( $self )  = shift;
  my @domains   = @{ $self->{"domains"} };
  my $set_ct    = 0;
  foreach my $set ( @domains) { 
    my $intree    = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.tree";
    my $outmat    = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.Rmat";
    my $refintree =  $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_REF_pruned.tree";
    my $refoutmat = $self->{"db"}->{"matrix"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_REF_pruned.Rmat";      
#    my $refnames  = $self->{"db"}->{"ref_align"} . "SSU_" . $set . "_refnames.list";    
    my @args = ();
    @args = ( "--slave", "--args", "$intree", "$outmat", "< " . $self->{"workdir"} . "tree_to_matrix.R");
    my $results = capture( "R @args" );
    if( $EXITVAL != 0 ){
      warn("Error running tree_to_matrix.R for $intree!\n");
      exit(0);
    }
    @args = ();
    @args = ( "--slave", "--args", "$refintree", "$refoutmat", "< " . $self->{"workdir"} . "tree_to_matrix.R");
    $results = capture( "R @args" );
    if( $EXITVAL != 0 ){
      warn("Error running tree_to_matrix.R for $refintree!\n");
      exit(0);
    }
  }
}

sub sim_format_matrix_to_phylip{
    my $self   = shift;
    my @domains   = @{ $self->{"domains"} };
    my $set_ct    = 0;
    foreach my $set ( @domains) {
	my $inmat  = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.Rmat";
	my $outmat = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.phymat";
	my $refinmat  = $self->{"db"}->{"matrix"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_REF_pruned.Rmat";
	my $refoutmat = $self->{"db"}->{"matrix"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_REF_pruned.phymat";
	my @args     = ();
	@args     = ("-i $inmat", "-o $outmat");
	my $results  = capture( "perl " . $self->{"workdir"} . "format_Rphylip.pl @args" );
	if( $EXITVAL != 0 ){
	    warn("Error formatting matrix for $inmat!\n");
	    exit(0);
	}
	@args     = ();
	@args     = ("-i $refinmat", "-o $refoutmat");
	$results  = capture( "perl " . $self->{"workdir"} . "format_Rphylip.pl @args" );
	if( $EXITVAL != 0 ){
	    warn("Error formatting matrix for $refinmat!\n");
	    exit(0);
	}
    }    
    return $self;
}

sub sim_run_mothur{
  #Run mothur
  my ( $self, $cutoff, $method ) = @_;
  my @domains   = @{ $self->{"domains"} };
  my $set_ct    = 0;
  foreach my $set ( @domains) {    
#    my $inmat    = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo.Rmat";
    my $inmat     = $self->{"db"}->{"matrix"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_pruned.phymat";
    my $refinmat    = $self->{"db"}->{"matrix"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_REF_pruned.phymat";
#    my $seqfile  = $self->{"db"}->{"SSU_reads"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . ".fa";
    my $seqfile  = $self->{"db"}->{"qc_seqs"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all_qc_seqs.fa";
    my $outdir   = $self->{"db"}->{"otudir"} . $set . "/";
    unless( -e $outdir ){
	make_path( $outdir );
    }
    #set.dir will be available in the next MOTHUR release
#    my $command  = "#set.dir(output=$outdir); read.dist(phylip=$inmat, cutoff=$cutoff); cluster(method=$method); bin.seqs(fasta=$seqfile)";
#    my $command  = "#read.dist(phylip=$inmat, cutoff=$cutoff); cluster(method=$method); bin.seqs(fasta=$seqfile)";
    #Having problems getting bin.seqs to work inside mothur, so we'll dot it ourselves later
    my $command  = "#set.dir(output=$outdir); read.dist(phylip=$inmat, cutoff=$cutoff); cluster(method=$method);";
    print "executing mother using $command\n";
    my $results  = capture( "mothur \"$command\"" );
    print $results ."\n";
    if( $EXITVAL != 0 ){
      warn("Error running mothur for $inmat!\n");
      exit(0);
    }  
    $command  = "#set.dir(output=$outdir); read.dist(phylip=$refinmat, cutoff=$cutoff); cluster(method=$method);";
    print "executing mother using $command\n";
    $results  = capture( "mothur \"$command\"" );
    print $results . "\n";
    if( $EXITVAL != 0 ){
      warn("Error running mothur for $refinmat!\n");
      exit(0);
    }
    
  }
}

# first arg: positive integer = SIZE of random subset to be returned;
# second arg: reference to array that contains elements of set to be sampled
# returns an array of length SIZE containing a random subset of the elements.
# Courtsey of S. Riesenfeld

sub get_random_subset ($$) {
    my ($size, $array_ref) = @_;
    my $range = scalar(@$array_ref);
    if ($range == $size) {
        my @ret_array=@{$array_ref};
        return \@ret_array;
    }

    unless ($range > $size) {die "In get_random_subset:  Size of array must be at least the desired size of random subset.\n";}
    my @set = (0..($range-1));
    my @chosen = ();
    foreach my $count (0..($size-1)) {
        my $random_num = int(rand($range-$count-1));    
        swap($count, ($random_num+$count+1),\@set);     
        my $index = $set[$count];
        # print "count: $count; at index $index: ". $$array_ref[$index]."\n";
        push(@chosen, $$array_ref[$index]);
    }
    return \@chosen;
}

# first arg: positive integer = index
# second arg: positive integer = index
# third arg: reference to an array
sub swap ($$$) {
    my ($index1, $index2, $array_ref) = @_;
    unless ( ($index1 < scalar(@$array_ref)) and ($index2 < scalar(@$array_ref)) ) {        die "Trying to reference element outside array boundary.\n";}
    # print "Swapping elements at indices $index1 and $index2.\n";
    my $temp = $$array_ref[$index1];
    $$array_ref[$index1]=$$array_ref[$index2];
    $$array_ref[$index2]=$temp;
    # print "at $index1: ". $$array_ref[$index1]."; at $index2: ".
        $$array_ref[$index2]."\n";
}

#Take a reference tree (newick format) and merge with a stockholm alignment in xrate format

sub build_xrate_alignment{
  my $self = shift;
  my @domains   = @{ $self->{"domains"} };
  my $set_ct    = 0;
  foreach my $set ( @domains) { 
    #This should be a masked alignment, but will need to change qc algo to make that work.
    my $align = $self->{"db"}->{"all_align"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all.fa";
    my $fmtalign = $self->{"db"}->{"all_align"} . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_all_fmt.fa"; 
    my $reftree = $self->{"db"}->{"tree"}  . $self->{"sample"}->{"name"} . "_SSU_" . $set . "_FT_pseudo_REF_pruned.tree";

    my $tree = ();
    open( TREE, $reftree ) || die "Can't open $reftree for xrate\n";
    $tree = <TREE>;
    close TREE;

    open ( ALIGN, $align ) || die "Can't open $align for xrate parsing\n";
    open ( OUT, ">$fmtalign" ) || die "Can't open $fmtalign for write\n";
    my $head = <ALIGN>;
    print OUT $head;
    #might need to be careful about the spacing between NH and $tree
    #need to ensure that the alignment names are consistant with tips names in tree
    print OUT "#=GF NH  " . $tree . "\n";
    while( <ALIGN> ){
      my $line = $_;
      print OUT $line;
    }
    close ALIGN;
    close OUT;
  }
}




1;
