#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %opt;
getopts ( "hf:o:", \%opt );
die "Usage: $0 -f file_name.txt -o output_file.txt\n" if $opt{h} or !$opt{f} or !$opt{o};

print "build tsid hash\n";
my %tsid;

open IN, "<$opt{f}";

while ( my $line = <IN> )
{
    my ($id) = $line =~ /transcript_id "(.+?)"/;
    if ( $line !~ /\tmRNA\t/ )
    {
        $tsid{$id} = "" unless $tsid{$id};
        $tsid{$id} .= $line;
    }
}

print "create mRNA line\n";
open OUT, ">$opt{o}";
while ( my ($k, $v) = each %tsid )
{
    my @lines  = split /\n/, $v;
    my @sorted = map { $_->[0] } sort { $a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] } map { [$_, /(^\S+)\s+\S+\s+\S+\s+(\d+)/] } @lines;
    $v = join "\n", @sorted;
    my @exons = $v =~ /\t(exon)\t/g;
# no exon lines
    if ( $#exons == -1 )
    {
        my ($start, $end) = $v =~ /CDS\s(\d+)\s(\d+)/;
        my ($pre, $post) = $sorted[0] =~ /^(\S+\s\S+)\s\S+\s\d+\s\d+\s\S+\s(.+)/;
        $v = "$pre\tmRNA\t$start\t$end\t.\t$post\n$pre\texon\t$start\t$end\t.\t$post\n$v\n";
        print OUT $v;
    }
# single exon genes
    elsif ( $#exons == 0 )
    {
        my ($start, $end) = $v =~ /exon\s(\d+)\s(\d+)/;
        my ($pre, $post) = $sorted[0] =~ /^(\S+\s\S+)\s\S+\s\d+\s\d+\s\S+\s(.+)/;
        $v = "$pre\tmRNA\t$start\t$end\t.\t$post\n$v\n";
        print OUT $v;
    }
    else
    {
# multiple exon genes
        my ($start, $end) = $v =~ /exon\s(\d+).+exon\s\d+\s(\d+)/s;
        my ($pre, $post) = $sorted[0] =~ /^(\S+\s\S+)\s\S+\s\d+\s\d+\s\S+\s(.+)/;
        print "v:\n$v\n" unless $start;
        $v = "$pre\tmRNA\t$start\t$end\t.\t$post\n$v\n";
        print OUT $v;
    }
}
close OUT;
