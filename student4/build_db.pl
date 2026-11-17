#!/usr/bin/env perl
# build_db.pl - turn the marine peptide tables into one SQLite database (as SQL text).
#
# Usage: perl student4/build_db.pl data/final > data/final/marine.sql
#        sqlite3 data/final/marine.db < data/final/marine.sql
#
# Joins on pep_id:  discovery_marine.tsv (organism, group, sources)
#                   discovery_marine_props.tsv (BioPerl properties)
#                   discovery_marine_scores_et.tsv (anticancer probability, call, rank)
#                   marine_acp_benchmark.tsv (known marine ACPs, flagged)
use strict;
use warnings;

my $dir = shift // 'data/final';

sub read_tsv {
    my ($file) = @_;
    open my $fh, '<', $file or die "Cannot open $file: $!\n";
    chomp(my $head = <$fh>);
    my @cols = split /\t/, $head;
    my %rows;
    while (<$fh>) {
        chomp;
        next unless length;
        my %r;
        @r{@cols} = split /\t/, $_, -1;
        $rows{ $r{pep_id} } = \%r;
    }
    close $fh;
    return \%rows;
}

sub sq {                       # SQL string literal
    my $s = defined $_[0] ? $_[0] : '';
    $s =~ s/'/''/g;
    return "'$s'";
}
sub num {                     # SQL number or NULL
    my $n = $_[0];
    return (defined $n && $n =~ /^-?\d+(?:\.\d+)?$/) ? $n : 'NULL';
}

my $info  = read_tsv("$dir/discovery_marine.tsv");
my $props = read_tsv("$dir/discovery_marine_props.tsv");
my $score = read_tsv("$dir/discovery_marine_scores_et.tsv");
my $bench = read_tsv("$dir/marine_acp_benchmark.tsv");

print <<'SQL';
PRAGMA foreign_keys = OFF;
DROP TABLE IF EXISTS peptides;
CREATE TABLE peptides (
    pep_id             TEXT PRIMARY KEY,
    sequence           TEXT NOT NULL,
    length             INTEGER,
    organism           TEXT,
    taxon_group        TEXT,
    sources            TEXT,
    source_ids         TEXT,
    mol_weight         REAL,
    pI                 REAL,
    net_charge         REAL,
    gravy              REAL,
    hydrophobic_pct    REAL,
    cationic_pct       REAL,
    anionic_pct        REAL,
    boman              REAL,
    aliphatic_index    REAL,
    hydrophobic_moment REAL,
    acp_probability    REAL,
    acp_call           TEXT,
    acp_rank           INTEGER,
    is_benchmark       INTEGER DEFAULT 0
);
BEGIN;
SQL

my $n = 0;
for my $id (sort keys %$info) {
    my ($i, $p, $s) = ($info->{$id}, $props->{$id}, $score->{$id});
    unless ($p && $s) { warn "skipped $id: missing from props or scores\n"; next; }
    printf "INSERT INTO peptides VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%d);\n",
        sq($id), sq($i->{sequence}), num($i->{length}), sq($i->{organism}), sq($i->{group}),
        sq($i->{sources}), sq($i->{source_ids}),
        (map { num($p->{$_}) } qw(mol_weight pI net_charge gravy hydrophobic_pct cationic_pct
                                  anionic_pct boman aliphatic_index hydrophobic_moment)),
        num($s->{anticancer_probability}), sq($s->{call}), num($s->{rank}),
        exists $bench->{$id} ? 1 : 0;
    $n++;
}

print <<'SQL';
COMMIT;
CREATE INDEX idx_group ON peptides(taxon_group);
CREATE INDEX idx_prob  ON peptides(acp_probability DESC);
SQL
print STDERR "Wrote $n peptides\n";
