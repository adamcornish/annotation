#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %opt;
getopts ( "hs:f:1:2:", \%opt );
die "Usage: $0 -f v2_annotation.gtf -s seq_dir -1 sam_dir_1 -2 sam_dir_2\n" if $opt{h} or !$opt{f} or !$opt{1} or !$opt{2} or !$opt{s};

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
    chomp ( my $tail1 = `tail -n1 $opt{1}/$k.sam` );
    chomp ( my $tail2 = `tail -n1 $opt{2}/$k.sam` );
    chomp ( my $fasta = `cat $opt{s}/$k.fasta` );
    my @cut1 = split /\t/, $tail1;
    my @cut2 = split /\t/, $tail2;
    my ($rc1,$chr1,$pos1,$tmp) = ($cut1[1],$cut1[2],$cut1[3],$cut1[9]);
    my ($rc2,$chr2,$pos2)      = ($cut2[1],$cut2[2],$cut2[3]);
    my ($seq) = $fasta =~ /.+?\n(.+)/s;
    $seq =~ s/\s//g;
    my $rev  = reverse $tmp;
    my $newv = "";
    my ($start, $stop) = $v =~ /mRNA\s(\d+)\s(\d+)/;
    $rev =~ tr/ACGTacgt/TGCAtgca/;
    if ( $pos1 and $pos2 )
    {
        printf ( "%5d/$total: $k\n", $i++ );
        if ( $pos1 < $pos2 or ( $pos1 == $pos2 and substr($seq,0,100) eq substr($tmp,0,100) ) )
        {
            my @lines = split /\n/, $v;
            foreach my $line ( @lines )
            {
                my @cols = split /\t/, $line;
                $cols[3] -= $start - $pos1;
                $cols[4] -= $start - $pos1;
                $cols[0] = $chr1;
                $newv .= join "\t", @cols;
                $newv .= "\n";
            }
        }
        elsif ( $pos1 > $pos2 or ( $pos1 == $pos2 and substr($seq,0,100) eq substr($rev,0,100) ) )
        {
            my @exons;
            my @CDS;
            if ( $v =~ /\t\+\t/ ) { $v =~ s/\t\+\t/\t-\t/; }
            else { $v =~ s/\t-\t/\t+\t/; }
            my @lines = split /\n/, $v;

            chomp ( my $seq = `tail -n1 $opt{1}/$k.sam | cut -f10` );
            my $seqlen = length ( $seq );
            $pos1 += $seqlen - 1;
            my ($pre,$post) = $lines[0] =~ /\S+(\s\S+)\s\S+\s\S+\s\S+(.+)/;
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
                push @tmp, "$chr1$pre\texon\t$coord1\t$coord2$post";
            }
            @final = reverse @tmp;
            foreach my $CDS ( @CDS )
            {
                my @split  = split /\t/, $CDS;
                my $coord1 = $pos1 - $split[0] - $split[1];
                my $coord2 = $pos1 - $split[0];
                push @tmp, "$chr1$pre\tCDS\t$coord1\t$coord2$post";
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
                $cols[0] = $chr1;
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
