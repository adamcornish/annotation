#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %opt;
getopts ( "hf:1:2:", \%opt );
die "Usage: $0 -f v2_annotation.gtf -1 sam_dir_1 -2 sam_dir_2\n" if $opt{h} or !$opt{f} or !$opt{1} or !$opt{2};

print "build gtf hash\n";
my @data = `cat $opt{f}`;
my %anno;
foreach my $line ( @data )
{
    my @split  = split /\t/, $line;
    my ($tsid) = $split[8] =~ /transcript_id "(.+?)"/;
    $anno{$tsid} = "" unless $anno{$tsid};
    $anno{$tsid} .= $line;
}

open OUT, ">zz_new.gtf";
open BAD, ">zz_unmapped.gtf";
my $i = 1;
chomp (my $total = `grep -cP "\\tmRNA\\t" $opt{f}`);
while ( my ($k,$v) = each %anno )
{
    chomp ( my $pos1  = `tail -n1 $opt{1}/$k.sam | cut -f4` );
    chomp ( my $pos2  = `tail -n1 $opt{2}/$k.sam | cut -f4` );
    chomp ( my $rc1   = `tail -n1 $opt{1}/$k.sam | cut -f2` );
    chomp ( my $rc2   = `tail -n1 $opt{2}/$k.sam | cut -f2` );
    my $newv  = "";
    my ($start, $stop) = $v =~ /mRNA\s(\d+)\s(\d+)/;
    if ( $pos1 and $pos2 )
    {
        printf ( "%5d/$total: $k\n", $i++ );
        if ( $pos1 < $pos2 )
        {
            my @lines = split /\n/, $v;
            foreach my $line ( @lines )
            {
                my @cols = split /\t/, $line;
                $cols[3] -= $start - $pos1;
                $cols[4] -= $start - $pos1;
                $newv .= join "\t", @cols;
                $newv .= "\n";
            }
        }
# if R2 is on the left side of R1, then the gene starts at the very end of R1 and heads <- that way until the beginning of R2
# if this is the case, we need to reverse the order of the exons. Exon n becomes exon 1, exon n-1 becomes exon 2, etc.
# flip the exons and then build new CDS? just count the number of exons in it is, yeah?
        elsif ( $pos1 > $pos2 )
        {
            my @exons;
            my @CDS;
            if ( $v =~ /\t\+\t/ ) { $v =~ s/\t\+\t/\t-\t/; }
            else { $v =~ s/\t-\t/\t+\t/; }
            my @lines = split /\n/, $v;

            chomp ( my $seq = `tail -n1 $opt{1}/$k.sam | cut -f10` );
            my $seqlen = length ( $seq );
            $pos1 += $seqlen - 1;
            my ($pre,$post) = $lines[0] =~ /(\S+\s\S+)\s\S+\s\S+\s\S+(.+)/;
            foreach my $line ( @lines ) 
            { 
                my @split = split /\t/, $line;
                if ( $line =~ /(?:\texon\t|\tCDS\t)/ )
                {
                    my $dist_from_start = $split[3] - $start;
                    my $length          = $split[4] - $split[3];
                    push @exons, "$dist_from_start\t$length" if $line =~ /\texon/;
                    push @CDS, "$dist_from_start\t$length"   if $line =~ /\tCDS/;
                }
            }
            my @tmp;
            my @final;
            foreach my $exon ( @exons )
            {
                my @split  = split /\t/, $exon;
                my $coord1 = $pos1 - $split[0] - $split[1];
                my $coord2 = $pos1 - $split[0];
                push @tmp, "$pre\texon\t$coord1\t$coord2$post";
            }
            @final = reverse @tmp;
            foreach my $CDS ( @CDS )
            {
                my @split  = split /\t/, $CDS;
                my $coord1 = $pos1 - $split[0] - $split[1];
                my $coord2 = $pos1 - $split[0];
                push @tmp, "$pre\tCDS\t$coord1\t$coord2$post";
            }
            $newv = join "\n", @tmp;
            $newv .= "\n";
        }
        else
        {
            my @lines = split /\n/, $v;
            foreach my $line ( @lines )
            {
                my @cols = split /\t/, $line;
                $cols[3] -= $start - $pos1;
                $cols[4] -= $start - $pos1;
                $newv .= join "\t", @cols;
                $newv .= "\n";
            }
        }
    }
    else
    {
        $newv = $v;
    }
    $newv =~ s/Chr/chr/g;
    ($rc1 == 4 or $rc2 == 4) ? print BAD $newv : print OUT $newv;
}
close OUT;
