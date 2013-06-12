#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %opt;
getopts ( "hf:g:", \%opt );
die "Usage: $0 -f genome.fa -g gtf_file.gtf\n" if $opt{h} or !$opt{f} or !$opt{g};

system ( "mkdir seqs" ) unless -d "seqs";

my %chrom;
my @genes;

print "Building chromosome hash\n";
open GENOME, "<$opt{f}";
my $chr = "";
while ( my $line = <GENOME> )
{
    if ( $line =~ /chr/i )
    {
        ($chr) = $line =~ />(.+)/;
        print "Working on $chr\n";
    }
    else
    {
        $line =~ s/\s//g;
        $chrom{$chr} .= $line;
    }
}
close GENOME;

print "Grabbing genes\n";
chomp ( @genes = `grep -P '\tmRNA\t' $opt{g}` );

print "Printing out sequences";
foreach my $gene ( @genes )
{
    my @split  = split /\t/, $gene;
    my $seq    = substr ( $chrom{$split[0]}, $split[3] - 1, $split[4] - $split[3] + 1 );
    my ($tsid) = $split[8] =~ /transcript_id "(.+?)"/;
    open OUT, ">seqs/$tsid.fasta";
    $seq =~ s/(.{60})/$1\n/g;
    print OUT ">$tsid\n$seq\n";
    close OUT;
}
