#!/usr/bin/perl -w

#format_Rphylip.pl - Convert an R matrix to Phylip format
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
