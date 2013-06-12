#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %opt;
getopts ( "hf:", \%opt );
die "Usage: $0 -f file_name.txt\n" if $opt{h} or !$opt{f};

open IN, "<$opt{f}";
open OUT, ">tmp.gtf";

while ( my $l1 = <IN> )
{
    my $l2 = <IN>;
    my @split = split /\t/, $l1;
    $split[2] = "mRNA";
    my $entry = join "\t", @split;
    $entry .= $l1;
# plus strand
    if ( $split[6] =~ /\+/ )
    {
        $split[2] = "start_codon";
        my $tmp = $split[4];
        $split[4] = $split[3] + 2;
        my $x = join "\t", @split;
        $entry .= $x;
        $entry .= $l2;
        $split[2] = "stop_codon";
        $split[4] = $tmp;
        $split[3] = $split[4] - 2;
        $x = join "\t", @split;
        $entry .= $x;
        print OUT $entry;
    }
# minus strand
    else
    {
        $split[2] = "stop_codon";
        my $tmp = $split[4];
        $split[4] = $split[3] + 2;
        my $x = join "\t", @split;
        $entry .= $x;
        $entry .= $l2;
        $split[2] = "start_codon";
        $split[4] = $tmp;
        $split[3] = $split[4] - 2;
        $x = join "\t", @split;
        $entry .= $x;
        print OUT $entry;
    }
}
close OUT;
close IN;
