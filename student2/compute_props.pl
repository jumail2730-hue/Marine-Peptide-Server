#!/usr/bin/env perl
# compute_props.pl - physicochemical property table for a peptide FASTA file (BioPerl).
# Usage:  perl compute_props.pl [--ph 7.0] input.fasta > properties.tsv
# FASTA headers from Student 1 look like  >MP000001|ACP|CancerPPD2
use strict;
use warnings;
use FindBin qw($Bin);
use lib $Bin;
use Getopt::Long;
use Bio::SeqIO;
use PeptideProps;

my $ph = 7.0;
GetOptions('ph=f' => \$ph) or die "Usage: $0 [--ph 7.0] input.fasta > properties.tsv\n";
my $file = shift or die "Usage: $0 [--ph 7.0] input.fasta > properties.tsv\n";

my $pp  = PeptideProps->new(ph => $ph);
my $in  = Bio::SeqIO->new(-file => $file, -format => 'fasta');
my @col = $pp->columns;

print join("\t", qw(pep_id label sequence), @col), "\n";
my ($ok, $skipped) = (0, 0);
while (my $s = $in->next_seq) {
    my ($id, $label) = split /\|/, $s->display_id;
    my $p = $pp->compute($s->seq);
    unless ($p) { $skipped++; print STDERR "skipped $id: non-standard residues\n"; next; }
    print join("\t", $id, $label // '', uc $s->seq, map { $p->{$_} } @col), "\n";
    $ok++;
}
print STDERR "Done: $ok peptides written, $skipped skipped (net charge at pH $ph)\n";
