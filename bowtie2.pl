#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %opt;
getopts ( "hd:", \%opt );
die "Usage: $0 -d dir\n" if $opt{h} or !$opt{d};

chomp ( my @files = `ls $opt{d}` );

chop $opt{d} if $opt{d} =~ /\/$/;
system ( "mkdir sam" ) unless -d "sam";

my $i = 1;
foreach my $file ( @files )
{
    my ($name) = $file =~ /(.+?)\./;
    my $c = "head -n1 $opt{d}/$file > tmp/$name.fasta && tail -n 100 $opt{d}/$file >> tmp/$name.fasta && nohup bowtie2 --mm --np 0 --very-fast -S sam/$name.sam -x MaSuRCA_Rhesus_Genome_v3_20130528 -f -U tmp/$name.fasta &";
    print "c: $c\n";
    $i++;
    sleep 15 if $i % 100 == 0;
    system ( $c );
}
