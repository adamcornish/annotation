#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %opt;
getopts ( "hg:s:o:", \%opt );
die "$0 -g <gtf file> -s chromosome_shifts.txt -o output.gtf" if $opt{h};
die "Need both files, moron.\n" unless $opt{g} and $opt{s};
die "$opt{g} doesn't exist. Exiting.\n" unless -e $opt{g};
die "$opt{s} doesn't exist. Exiting.\n" unless -e $opt{s};
die "You didn't specify an output file. Exiting.\n" unless $opt{o};

print "cat shifts\n";
my @shifts = `cat $opt{s}`;
print "cat gtf file\n";
#my @unsorted = `cat $opt{g}`;
my @anno = `cat $opt{g}`;
#print "sort the thing\n";
#my @anno = map { $_->[0] } sort { $a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] } map { [$_, /(^\S+)\s+\S+\s+\S+\s+(\d+)/] } @unsorted;
print "done\n";

for ( my $i = $#shifts; $i > -1; $i-- )
{
    print "working on this: $shifts[$i]";
    my $shift = $shifts[$i];
    $shift =~ s/\s{1,80}/\t/g;
    my (@col) = split /\t/, $shift; 
    for ( my $i = 0; $i < @anno; $i++ )
    {
        if ( $anno[$i] )
        {
            my @split     = split /\t/, $anno[$i];
            my $newline   = "";
            my $printable = 1;
            if ( $col[0] eq $split[0] )
            {
                if ( $col[3] =~ /R/i )
                {
                    $printable = 0 if ( $split[3] < ($col[1] + $col[2]) and $split[3] > $col[1] ) or ( $split[4] < ($col[1] + $col[2]) and $split[4] > $col[1] );
                    $split[3] -= $col[2] if $split[3] > $col[1];
                    $split[4] -= $col[2] if $split[4] > $col[1];
                }
                elsif ( $col[3] =~ /I/i )
                {
                    $split[3] += $col[2] if $split[3] > $col[1];
                    $split[4] += $col[2] if $split[4] > $col[1];
                }
                else
                {
                    die "Danger, Will Robinson!\n";
                }
            }
            $newline = join "\t", @split if $printable;
            $anno[$i] = $newline;
        }
    }
}
open OUT, ">$opt{o}";
print OUT @anno;
close OUT;
