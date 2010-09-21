#!/usr/bin/perl -w

use strict;
use Bio::AlignIO;
use Getopt::Long;
use Bio::Align::Utilities qw(:all);


#Var initializatoin
my ($input, $ingap, $output, $outgap, $flat, $aln_format, $min);

#Set a few default settings
$min                = 1;
$aln_format         = "fasta";

#Command line parameter space settings
GetOptions(
	   "i=s"   => \$input,              #Path to alignment object(s)
	   "j=s"   => \$ingap,              #Path to gaps file
	   "o=s"   => \$output,             #Output file path
	   "p=s"   => \$outgap,             #Output gaps file
	   "m=i"   => \$min,                #Minimum number of residues a column must contain to be retained
	   "f:s"   => \$aln_format,         #Alignment format (Default = FASTA)
	   "flat"  => \$flat,               #use flat to keep coords out of headers
);

print "processing $input with cutoff $min, will print output to $output\n";

#You're not a bonehead, right?
if(!$input || !$output){ die "You must specific an input and output file!\n" }

#Initialize the framework
my $in_aln  = Bio::AlignIO->new(-file => "$input",   -format => "$aln_format");
my $out_aln = Bio::AlignIO->new(-file => ">$output", -format => "$aln_format");
open( INGAP,  "$ingap"    ) || die "can't open input $ingap in align2profile_qc_Col.pl:$!\n";
unlink($outgap);
open( OUTGAP, ">>$outgap" ) || die "can't open output $outgap in align2profile_qc_Col.pl:$!\n";
my $aln = $in_aln->next_aln();
my $num_seq      = <INGAP>;
my $num_cov_str  = <INGAP>;
my @num_cov      = split(" ", $num_cov_str);
my $perc_cov_str = <INGAP>;
my @perc_cov     = split(" ", $perc_cov_str);
my $size         = @perc_cov;
my @Onum_cov     = ();
my @Operc_cov    = ();
my @remove       = ();
close INGAP;

my $i;
#Find the gappy columns
for( $i=0; $i<$size; $i++ ){
 if( $num_cov[$i] >= $min ){
   #This is a sufficiently filled column, preserve it
   push(@Onum_cov,  $num_cov[$i]);
   push(@Operc_cov, $perc_cov[$i]);
 } else {
   #Mark this column to be removed
   push(@remove, $i+1);
 } 
}
print $aln->length(), "\n";
print "original alignment: $size  reduced by: ",scalar(@remove), "\n";
push(@remove, $i+1);  #One beyond the end marks the last 

#Print the new gaps file
 printf OUTGAP $num_seq;
 my $format = ("%1.0f " x @Onum_cov)."\n";
 printf OUTGAP $format,    @Onum_cov;
 $format    = ("%1.5f " x @Operc_cov)."\n";
 printf OUTGAP $format,    @Operc_cov;
 close OUTGAP;

#Mask character so rows don't get deleted
$aln->map_chars('-','J');
$aln->map_chars('\.','O');

#Collect together only the good columns
my @aln_slices = ();
my $goodstart = 1;
my $goodend   = 1;
my $sum =0;
$i=1;
while ( scalar(@remove) > 0 ){
  my $badcol = shift(@remove);
  $i++;
  if( $goodstart == $badcol ){
    #The next residue is also bad, step forward
    $goodstart++;
  } else {
    $goodend = $badcol - 1;
    if( $goodstart <= $goodend ){
      #print "$goodstart - $goodend\n";
      my $slice = $aln->slice( $goodstart, $goodend );
      push( @aln_slices, $slice );
      $goodstart = $badcol+1;
    }
  }  
}

#Print the good rows to file
my $clean_aln = cat(@aln_slices);
print "new alignment: ", $clean_aln->length, "\n";
#UnMask characters so we print something sensible
$clean_aln->map_chars('J','-');
$clean_aln->map_chars('O','.');

if( $flat ){
  $clean_aln->set_displayname_flat();
}
$out_aln->write_aln($clean_aln);
