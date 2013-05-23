#!/usr/bin/perl
use warnings;
use strict;

# This scripts assumes the genes are sorted on column 1 and then column 4 (i.e. by chromosome and then coordinate)

print "read in the gtf file\n";
my @lines = `cat $ARGV[0]` unless $#ARGV != 0;

my %hash;
print "populate the hash\n";
foreach my $line ( @lines )
{
    my ($transcript_id) = $line =~ /transcript_id "(.+?)"/;
    $hash{$transcript_id} = "" unless $hash{$transcript_id};
    $hash{$transcript_id} .= $line;
}

$"="\n";

open OUT, ">fixed_frames.gtf";
while ( my ($k, $v) = each %hash )
{
    #check to see if it has CDS
    if ( $v =~ /CDS/ )
    {
        my @lines = split /\n/, $v;
       #plus strand
        if ( $v =~ /\t+\t/ )
        {
            for ( my $i = 0; $i < @lines; $i++ )
            {
                # fix the frames
                if ( $lines[$i] =~ /\tCDS\t/ )
                {
                    my @split = split /\t/, $lines[$i];
                    $split[7] = 0;
                    if ( $lines[$i-2] =~ /\tCDS\t/ )
                    {
                        my @previous_split = split /\t/, $lines[$i-2];
                        $split[7] = ( 3 - ( $previous_split[4] - $previous_split[3] - $previous_split[7] + 1 ) % 3 ) % 3;
                    }
                    $lines[$i] = join ( "\t", @split );
                }
            }
        }
       #minus strand
        else
        {
            for ( my $i = $#lines; $i >= 0; $i-- )
            {
                # fix the frames
                if ( $lines[$i] =~ /\tCDS\t/ )
                {
                    my @split = split /\t/, $lines[$i];
                    $split[7] = 0;
                    if ( $lines[$i+2] =~ /\tCDS\t/ )
                    {
                        my @previous_split = split /\t/, $lines[$i+2];
                        $split[7] = ( 3 - ( $previous_split[4] - $previous_split[3] - $previous_split[7] + 1 ) % 3 ) % 3;
                    }
                    $lines[$i] = join ( "\t", @split );
                }
            }
        }
        $v = join ( "\n", @lines );
        print OUT "$v\n";
    }
}
close OUT;
