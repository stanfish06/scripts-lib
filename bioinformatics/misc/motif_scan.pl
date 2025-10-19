#/usr/bin/perl
use strict;
use warnings;

use Path::Tiny;

my $dir = path(".");
my $file = $dir->child("test.fa");

my $content = $file->slurp_utf8();
my $motif = 'AT(C|A|G)G';
my @matches;
while($content =~ /${motif}/g) {
	push @matches, "match: $&| alt: $1| $-[0]-$+[0]";
}
print join("\n", @matches);
