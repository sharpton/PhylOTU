#!/usr/bin/perl -w

use strict;
use Getopt::Long;

#Convert a raw Rmatrix to a properly formatted phylip matrix. If -t is
#included, reference a esprit tag to sequence identifier to replace
#matrix tag ids with sequence ids. Note that for MOTHUR's bin_seqs we want the 
#sequence id to perfectly match sequence database headers....

my ($in_raw_matrix, $esprit_tag_2_seqid_lookup, $outmatrix );

GetOptions(
    'i=s' => \$in_raw_matrix,
    't=s' => \$esprit_tag_2_seqid_lookup,
    'o=s' => \$outmatrix,
    );

open( IN, $in_raw_matrix ) || die;
my @inlines = <IN>;
close IN;

open( OUT, ">$outmatrix" ) || die;

my $heads = shift( @inlines); #top line is header, so don't count it
my $ntaxa = @inlines; 
print OUT "\t$ntaxa\n";
#Fast format, no lookup needed
unless( defined( $esprit_tag_2_seqid_lookup ) ){
    foreach my $line ( @inlines ){
	chomp $line;
	print OUT "$line\n";
    }
}
#Convert esprit tag ids to seq ids
else{
    open( LOOK, $esprit_tag_2_seqid_lookup ) || die;
    my %lookup = ();
    while(<LOOK>){
	chomp $_;
	my ($tagid, $seqid, @extra) = split("\t", $_);
	$lookup{$tagid} = $seqid;
    }
    close LOOK;
    foreach my $line( @inlines ){
	chomp $line;
	if( $line =~ m/^\"(.*)\"/ ){
	    my $tag   = $1 - 1; #R is 1 indexed, esprit is 0 indexed
	    my $seqid = $lookup{$tag};
	    $line =~ s/^\"(.*)\"/\"$seqid\"/;
	    print OUT "$line\n";
	}
	else{
	    warn("Can't replace tag id with sequence id!\n");
	    exit(0);
	}
    }
}
close OUT;
