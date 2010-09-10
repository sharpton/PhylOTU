#!/usr/bin/perl -w

use strict;
use Bio::AlignIO;
use Getopt::Long;
use Bio::Align::Utilities qw(:all);

#CONSIDERATIONS/RESTRICTIONS
#Currently, all thresholds are inclusive (if you equal threshold, you stay)
#For all gap based cutoffs, only gaps relative to the profile are considered (-). "." are ignored for cutoff purposes.
#Currently, score/p-vals can only be parsed from INFERNAL output.  Will update as needed.
#Alignment file can only contain a single alignment
#Currently, alignment masking is the excition of columns, so not ideal for protein space (should replace with N|?)

#Var initializatoin
my ($input, $output, $seq_scores_tab, $lookup_list, $flat);
my ($seq_len_cut, $aln_seq_score_cut, $aln_seq_p_cut, $aln_seq_score_rat, $gap_ratio_cut, $mm_ratio_cut, $pcnt_cov_cut, $ngaps_internal_cut, $gappy_col_cut);
my ($gap_window, $mm_window, $aln_format);

#Set a few default settings
$seq_len_cut        = 0;
$aln_seq_score_cut  = 0.0;
$aln_seq_p_cut      = 0.0;
$aln_seq_score_rat  = 1.0;
$pcnt_cov_cut       = 15;
$ngaps_internal_cut = 50;
$gap_ratio_cut      = 50;   #Currently not configured
$mm_ratio_cut       = 50;   #Currently not configured
$mm_window          = 50;   #IN BP/AA POSITIONS
$gappy_col_cut      = 0.75; 
$gap_window         = $mm_window;
$aln_format         = "fasta";

#Command line parameter space settings
GetOptions(
	   "i=s"   => \$input,              #Path to alignment object(s)
	   "o=s"   => \$output,             #Output file path
	   "s:s"   => \$seq_scores_tab,     #Formatted for INFERNAL verbose output
	   "l:s"   => \$lookup_list,        #Path to list of sequence ids to prune from alignment
	   "lc:i"  => \$seq_len_cut,        #Minimum sequence length cutoff
	   "sc:i"  => \$aln_seq_score_cut,  #Minimum alignment score for each sequence cutoff
	   "sr:i"  => \$aln_seq_score_rat,  #Minimum bit score per base in sequence cutoff
	   "pc:i"  => \$aln_seq_p_cut,      #Minimum p-val/e-val for each sequence cutoff
	   "grc:i" => \$gap_ratio_cut,      #Minimum number of gaps/window for each sequence cutoff	   "mrc=i" => \$mm_ratio_cut,       #Minimum number of mismatches/window per sequence cutoff
	   "cov:i" => \$pcnt_cov_cut,       #Percentage of sequence positions that must be covered by model/concensus (minimum)
	   "ngc:i" => \$ngaps_internal_cut, #Maximum number of sequence bounded gaps tolerated per sequence (TOTAL GAPS)
	   "gcc:i" => \$gappy_col_cut,      #Maximum number of gaps in a column divided by sequences in an alignment
	   "mmw:i" => \$mm_window,          #Size of mismatch window (Default = 50)
	   "gpw:i" => \$gap_window,         #Size of gap window (Default = mm_window)
	   "f:s"   => \$aln_format,         #Alignment format (Default = FASTA)
	   "flat"  => \$flat,             #use flat to keep coords out of headers
);

print "processing $input with format $aln_format and referencing $seq_scores_tab, will print output to $output\n";

#You're not a bonehead, right?
if(!$input || !$output){
  die "You must specific an input and output file!\n"
}
if(!$seq_scores_tab && ( $aln_seq_score_cut || $aln_seq_p_cut) ) {
  warn "You didn't specify a place to obtain sequence score or pvalue, so the 'sc' and 'pc' cutoffs you set will be ignored!\n"
}

#Initialize the framework
my $in_aln  = Bio::AlignIO->new(-file => "$input",   -format => "$aln_format");
my $out_aln = Bio::AlignIO->new(-file => ">$output", -format => "$aln_format");
my %seq_scores = ();
if (defined($seq_scores_tab)){
  %seq_scores = %{ process_seq_score_tab($seq_scores_tab) };
}
my %init_prunes = ();
if(defined($lookup_list)){
  %init_prunes = %{ process_lookup_list($lookup_list) };
}

#Let's loop through the alignment and look at each sequence, one at a time.  Prune out the bad ones!
while( my $aln = $in_aln->next_aln() ){
  #make sure all sequences are the same length
  if( !$aln->is_flush() ){
    warn "Alignment isn't flush, passing!\n";
    next;
  }
  my $nseqs  = $aln->num_sequences();
  my $alnlen = $aln->length();
  #build a sequence, gap, insert map of the alignment
  #my %pre_alignment_map = %{ build_alignment_map( $aln ) };
  foreach my $seq ($aln->each_seq) {
    my $to_prune   = 0;
    my @reasons    = ();
    my $id         = $seq->display_id();
    my $sequence   = $seq->seq();
    #sequence length shouldn't count gaps/inserts, so copy $sequence and count seq length
    my $seq_seq    = $sequence;
    $seq_seq       =~ s/\-//g;
    $seq_seq       =~ s/\.//g;
    my $seqlen     = length( $seq_seq );
    my $residues   = $sequence;
    $residues      =~ s/(\-|\.)|//g;
    my $n_residues = length($residues);
#    my @seqarray   = split("", $sequence);
#    for(my $i=0; $i < $alnlen; $i++){
#      $pre_alignment_map{ $i+1 }{ $seqarray[$i] }++;
#    }     
    #SEQUENCE IDS A LOOKUP FILE SAYS TO PASS ON
    foreach my $prune_id ( keys( %init_prunes ) ){
     if($id =~ m/$prune_id/){
       $to_prune = 1;
       last;
     }
    }
    #SEQUENCE LENGTH
    if($seqlen < $seq_len_cut){
      $to_prune = 1;
    }
    #INTERNAL GAPS
    if(defined($ngaps_internal_cut)){
      my $cpseq = $sequence;
      #only want internal gaps (bounded by sequence) so let's remove leading and lagging
      $cpseq =~ s/^(\-|\.)+//;
      $cpseq =~ s/(\-|\.)+$//;
      #now count the total number of gaps
      my $ngaps = $cpseq =~ tr/\-/\-/;
      if ($ngaps > $ngaps_internal_cut){
	$to_prune = 1;
      }
    }
    #SCORE AND P-VAL BASED CUTOFF PROCESING
    if(defined($seq_scores_tab)){
      if(defined($seq_scores{$id})){  #IF NOT, ASSUME IT'S REFERENCE SEQUENCE USED TO BUILD MODEL INCLUDED IN ALIGNMENT (MEANS THAT IT'S A GOOD READ)
	if(defined($aln_seq_score_cut) && $seq_scores{$id}{"score"} < $aln_seq_score_cut){
	  $to_prune = 1;
	}
	if(defined($aln_seq_p_cut) && $seq_scores{$id}{"pval"} < $aln_seq_p_cut){
	  $to_prune = 1;
	}
	if(defined($aln_seq_score_rat) ){
	  my $score_per_residue = $seq_scores{$id}{"score"} / $n_residues;
	  if( $score_per_residue < $aln_seq_score_rat ){
	    $to_prune = 1;
	  }
	}
      }
    }
    #DO THE PRUNING OF THE SEQUENCES HERE
    if($to_prune){
      #print "Pruning $id\n";
      $aln->remove_seq($seq);
    }
  }
  #%pre_alignment_map = ();
  
  #POST PRUNING ALIGNMENT CLEANUP/MASKING
  #Now iterate through the sequences that passed and remove bad columns from the alignment.
  #make sure all sequences are the same length
  if( !$aln->is_flush() ){
    warn "Alignment is no longer flush, passing!\n";
    next;
  }
  my $post_nseqs  = $aln->num_sequences();
  my $postalnlen = $aln->length();
  #Build a processed alignment map
  my %post_alignment_map = %{ build_alignment_map( $aln ) };
  #purge gappy columns using %post_alignment_map
  my @gappy_cols = ();
  foreach my $pos( sort { $a <=> $b } ( keys( %post_alignment_map ) ) ){
    if(defined($post_alignment_map{$pos}{"."}) && $post_alignment_map{$pos}{"."}/$post_nseqs > $gappy_col_cut){
      #print "$pos\t$post_alignment_map{$pos}{'.'}\n";
      push(@gappy_cols, $pos);
    }
  }
  my $n_gappy_cols = @gappy_cols;
  my $gap_start = $gappy_cols[0];  #Keep all start and stop vars in aln coordinate space (1 based)
  my $gap_stop  = 1;
  my $aln_start = 1;
  my $aln_stop  = $postalnlen;
  my @aln_slices = ();
  $aln->map_chars('-', 'N');  #Must mask gaps for splicing to work later
  for ( my $n = 0; $n < $n_gappy_cols; $n++ ){
    my $slice;
    if ( exists( $gappy_cols[ $n + 1] ) && $gappy_cols[ $n ] + 1 == $gappy_cols[ $n + 1 ] ){
      #then monotonic increase, move on
      next;
    }
    else{
      #check to see if gaps start at beginning of aln, if so, reset gap_start and move on
      if( $aln_start == $gap_start ){
	$gap_start = $gappy_cols[ $n + 1 ];
	$gap_stop  = $gappy_cols[ $n ];
	next;
      }
      #deal with the lefthand-most aln block
      elsif( !@aln_slices ){
	$slice = $aln->slice( $gap_stop, $gap_start - 1 );
	push( @aln_slices, $slice );
	$gap_start = $gappy_cols[ $n + 1 ];
	$gap_stop  = $gappy_cols[ $n ];
      }
      #deal with the righthand-most aln block
      elsif( $n+1 == $n_gappy_cols && $gappy_cols[$n] < $aln_stop ){
	#need to process the aln blocks on both sides of the gappy column block
	#get the left side:
	$slice = $aln->slice( $gap_stop + 1, $gap_start - 1 );
	push( @aln_slices, $slice );
	$gap_stop = $gappy_cols[$n];
	#get the right side:
	$slice = $aln->slice( $gap_stop + 1, $aln_stop );
	push( @aln_slices, $slice );
      }
      else{
	$slice = $aln->slice( $gap_stop + 1, $gap_start - 1);
	push( @aln_slices, $slice );
	$gap_start = $gappy_cols[ $n + 1];
	$gap_stop  = $gappy_cols[$n];
      }
    }
  }
  %post_alignment_map=();
  #now stitch alignment together
  my $clean_aln = cat(@aln_slices);
  #map back to gap chars
  $aln->map_chars('N','-');
  $clean_aln->map_chars('N','-');
  #Build final alignment map
  my %final_alignment_map = %{ build_alignment_map( $clean_aln ) };
  produce_aln_profile_data( \%final_alignment_map );
  if( $flat ){
    $clean_aln->set_displayname_flat();
  }
  $out_aln->write_aln($clean_aln);
}

####################################
#SUBROUTINES
####################################
sub produce_aln_profile_data{
  my %alignment_map = %{ $_[0] };
  foreach my $pos ( sort  { $a <=> $b } ( keys ( %alignment_map ) ) ){
    my $coverage = 0;
    foreach my $char ( keys( %{ $alignment_map{ $pos } } ) ){
      if ($char eq "-" || $char eq "."){
	next;
      }
      $coverage = $coverage + $alignment_map{ $pos }{ $char };
    }
    print "$pos\t$coverage\n";
  }
}

sub build_alignment_map{
  my $aln = shift;  #An AlignI object;
  my %alignment_map = ();
  my $alnlen = $aln->length();
  foreach my $seq ($aln->each_seq) {
    my $sequence = $seq->seq();
    my @seqarray = split("", $sequence);
    for(my $i=0; $i < $alnlen; $i++){
      $alignment_map{ $i+1 }{ $seqarray[$i] }++;
    }     
  }
  return \%alignment_map
}

sub process_lookup_list{
  my $file = shift;
  my %init_prunes = ();
  open(IN, $file) || die "can't open $file for score parsing: $!\n";
  while(<IN>){
    chomp $_;
    my $id = $_;
    $init_prunes{$id}{"count"}++;
    if($init_prunes{$id}{"count"} > 1){
      warn "Found lookup list id $_ inside lookup list file more than once!\n"
    }
  }
  return \%init_prunes;
}

sub process_seq_score_tab{
  my $file = shift;
  open(SCORES, $file) || die "can't open $file for score parsing: $!\n";
  my %seq_scores = ();
  while(<SCORES>){
    chomp $_;
    if($_ =~ m/^\#/ || $_ =~ m/^$/ ){
      next;
    }
    my @data  = split(" ", $_);
    my $read  = $data[1];
    my $score = $data[3];
    my $pval  = $data[5];
    $seq_scores{$read}{"count"}++;
    if($seq_scores{$read}{"count"} > 1){
      warn "Sequence $read has more than one set of results. Overwriting previous results which are\n\tscore\t" .
	$seq_scores{$read}{"score"} . "\n\tpval\t" .
	$seq_scores{$read}{"pval"} . "\n";
    }
    $seq_scores{$read}{"score"} = $score;
    $seq_scores{$read}{"pval"} = $pval;    
  }
  return \%seq_scores;
  close SCORES;
}
