#!/usr/bin/perl
use warnings;
use strict;

print "read in the gff file\n";
my @gff = `cat same_number_of_CDS.gff`;

my %hash;
print "populate the hash\n";
foreach my $line ( @gff )
{
    my ($transcript_id) = $line =~ /(NM.+?)(?:;|\s)/;
    $hash{$transcript_id} = "" unless $hash{$transcript_id};
    $hash{$transcript_id} .= $line;
}

#fix the last column so that it has transcript_id,             gene_id, gene_symbol and gene_description
#                                   gene_symbol_transcript_01  col_3,   col_4,          col_5
my %needed;
print "fixing the last column\n";
while ( my ($k, $v) = each %hash )
{
# check to see if it's in our needed list
    if ( `grep $k missing_no_sequence.tab` )
    {
        my @lines = split /\n/, $v;
        foreach my $line ( @lines )
        {
            my @split    = split /\t/, $line;
            my $pop      = pop @split;
            my $pre      = join "\t", @split;
            my $anno     = `grep $k missing_no_sequence.tab`;
            my ($NM, $NP, $gene_id, $gene_symbol, $gene_description) = split /\t/, $anno;
            chop ( $gene_description );
            my $last_col = "transcript_id \"$k\"; gene_id \"$gene_id\"; gene_symbol \"$gene_symbol\"; gene_description \"$gene_description\";";
            $needed{$k} = "" unless $needed{$k};
            $needed{$k} .= "$pre\t$last_col\n";
        }
    }
}
print "done fixing last column\n";
$"="\n";

open OUT, ">same_number_of_CDS_with_needed.gtf";
open BAD_CDS, ">supposedly_bad_CDS.txt";
open BAD_EXONS, ">supposedly_bad_exons.txt";
while ( my ($k, $v) = each %needed )
{
    #check to see if it has CDS
    if ( $v =~ /CDS/ )
    {
        my @lines = split /\n/, $v;

       #plus strand
        if ( $v =~ /mRNA\s\S+\s\S+\s\S+\s\+/ )
        {
            for ( my $i = 1; $i < @lines; $i++ ) # start at 1 because the first line is the mRNA line
            {
                # add a start codon
                if ( $lines[$i] =~ /CDS/ and $lines[$i-1] =~ /exon/ )
                {
                    my @split = split /\t/, $lines[$i];
                    my $start_codon  = "$split[0]\t$split[1]\tstart_codon\t$split[3]\t".($split[3]+2)."\t.\t+\t0\t$split[8]";
                    splice ( @lines, $i, 0, $start_codon );
                    $i++;
                }
                # fix the frames
                if ( $lines[$i] =~ /CDS/ )
                {
                    my @split = split /\t/, $lines[$i];
                    # First CDS in the gene, so its frame is 0
                    if ( $lines[$i-1] !~ /CDS/ )
                    {
                        $lines[$i] =~ s/\.\ttranscr/0\ttranscr/;
                    }
                    else
                    {
                        my @previous_split = split /\t/, $lines[$i-1];
                        $split[7] = ( 3 - ( $previous_split[4] - $previous_split[3] - $previous_split[7] + 1 ) % 3 ) % 3;
                        $lines[$i] = join ( "\t", @split );
                    }
                }
                # add a stop codon
                if ( $i == $#lines )
                {
                    my @split = split /\t/, $lines[$i];
                    my $stop_codon = "$split[0]\t$split[1]\tstop_codon\t".($split[4]-2)."\t$split[4]\t.\t+\t0\t$split[8]";
                    push @lines, $stop_codon;
                    $i++;
                }
            }
        }
       #minus strand
        else
        {
            for ( my $i = $#lines; $i >= 0; $i-- )
            {
                # add a stop codon
                if ( $lines[$i] =~ /CDS/ and $lines[$i-1] =~ /exon/ )
                {
                    my @split = split /\t/, $lines[$i];
                    my $stop_codon  = "$split[0]\t$split[1]\tstop_codon\t".($split[3]-2)."\t$split[3]\t.\t+\t0\t$split[8]";
                    splice ( @lines, $i, 0, $stop_codon );
                    $i++;
                }
                # fix the frames
                if ( $lines[$i] =~ /CDS/ )
                {
                    my @split = split /\t/, $lines[$i];
                    # First CDS in the gene, so its frame is 0
                    if ( $i == $#lines )
                    {
                        $lines[$i] =~ s/\.\ttranscr/0\ttranscr/;
                    }
                    else
                    {
                        my @previous_split = split /\t/, $lines[$i+1];
                        $split[7] = ( 3 - ( $previous_split[4] - $previous_split[3] - $previous_split[7] + 1 ) % 3 ) % 3;
                        $lines[$i] = join ( "\t", @split );
                    }
                }
                # add a start codon
                if ( $i == $#lines )
                {
                    my @split = split /\t/, $lines[$i];
                    my $start_codon  = "$split[0]\t$split[1]\tstart_codon\t".($split[4]+2)."\t$split[4]\t.\t-\t0\t$split[8]";
                    push @lines, $start_codon;
                }
            }
        }
        $v = join ( "\n", @lines );
        @lines              = split /\n/, $v;
        my $big_number      = 500_000;
        my $dont_print_gene = 0;
        my ($start_coord)   = $v =~ /start_codon\s(\d+)/;
        my ($stop_coord)    = $v =~ /stop_codon\s(\d+)/;
        my $gene            = "$lines[0]\n";
        my @start           = split /start_codon/, $v;
        my @stop            = split /stop_codon/, $v;
        # Start at 1 because $lines[0] is the mRNA line
        for ( my $i = 1; $i < @lines; $i++ )
        {
            my @split   = split /\t/, $lines[$i];
            my @n_split = split /\t/, $lines[$i+1] unless $i == $#lines;
            # exon is too far upstream of the start codon and is + strand
            if    ( $split[2] eq "exon" and $split[6] eq "+" and ($start_coord - $split[4] > $big_number) ) { print BAD_EXONS "$lines[$i]\n"; }
            # exon is too far upstream of the start codon and is - strand
            elsif ( $split[2] eq "exon" and $split[6] eq "-" and ($split[3] - $start_coord > $big_number) ) { print BAD_EXONS "$lines[$i]\n"; }
            # exon is too far downstream of the stop codon and is + strand
            elsif ( $split[2] eq "exon" and $split[6] eq "+" and ($split[3] - $stop_coord > $big_number) ) { print BAD_EXONS "$lines[$i]\n"; }
            # exon is too far downstream of the stop codon and is - strand
            elsif ( $split[2] eq "exon" and $split[6] eq "-" and ($stop_coord - $split[4] > $big_number) ) { print BAD_EXONS "$lines[$i]\n"; }
            # CDS is too far away from the next CDS
            elsif ( $split[2] eq "CDS" and $n_split[2] eq "CDS" and ($n_split[3] - $split[4] > $big_number) )
            {
                $dont_print_gene = 1;
                print BAD_CDS $v;
                last;
            }
            else { $gene .= "$lines[$i]\n"; }
        }
        unless ( $dont_print_gene )
        {
            my ($new_start) = $gene =~ /.+?exon\s(\d+)/s;
            my ($new_stop) = $gene =~ /.+exon\s\S+\s(\d+)/s;
            my ($old_start, $old_stop) = $gene =~ /mRNA\s(\d+)\s(\d+)/;
            $gene =~ s/$old_start/$new_start/;
            $gene =~ s/$old_stop/$new_stop/;
           #my ($transcript_id) = $gene =~ /transcript_id "(.+?)"/;
           #$new_hash{$transcript_id} = $gene;
            print OUT $gene;
        }
    }
}
close OUT;
close BAD_CDS;
close BAD_EXONS;
