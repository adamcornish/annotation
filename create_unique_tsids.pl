#!/usr/bin/perl
use warnings;
use strict;

my $file = shift;
open IN, "<$file";
my %hash;
my $tsid;
while ( my $line = <IN> )
{
    my ($pre, $post) = $line =~ /(.+?").+?(".+)/;
    my ($gsym) = $line =~ /gene_symbol "(.+?)"/;
    if ( $line =~ /\tmRNA\t/ )
    {
        $tsid = "";
        for ( my $i = 1; 1; $i++ )
        {
            $tsid = sprintf ( "%s_transcript_%02d", $gsym, $i );
=cut
            if ( $i > 5 )
            {
                print "hashtsid: $hash{$tsid}\n";
                print "line: $line\n";
                print "tsid: $tsid\n";
                my $sin = <STDIN>;
            }
=cut
            if ( !$hash{$tsid} )
            {
                $hash{$tsid} = "$pre$tsid$post\n";
                last;
            }
        }
    }
    else { $hash{$tsid} .= "$pre$tsid$post\n"; }
}

open OUT, ">$file.uniq";
while ( my ($k, $v) = each %hash )
{
    print OUT $v;
}
close OUT;
