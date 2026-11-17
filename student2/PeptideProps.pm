package PeptideProps;
# PeptideProps.pm - physicochemical properties of a peptide, using BioPerl.
#   use PeptideProps;
#   my $pp = PeptideProps->new(ph => 7.0);
#   my $p  = $pp->compute('GIGKFLHSAKKFGKAFVGEIMNS');   # hashref of properties
# Properties: length, mol_weight (Da), pI (EMBOSS pK set), net_charge at pH,
#   gravy (Kyte-Doolittle), hydrophobic_pct (A C F I L M V W), cationic_pct (K R),
#   anionic_pct (D E), boman = -(mean transfer free energy, Radzicka & Wolfenden;
#   APD prints the opposite sign), aliphatic_index (Ikai 1980),
#   hydrophobic_moment (Eisenberg, 100 degrees).
use strict;
use warnings;
use Bio::Seq;
use Bio::Tools::SeqStats;
use Bio::Tools::pICalculator;

our @COLUMNS = qw(length mol_weight pI net_charge gravy hydrophobic_pct
                  cationic_pct anionic_pct boman aliphatic_index hydrophobic_moment);

my %KD = (A=>1.8, R=>-4.5, N=>-3.5, D=>-3.5, C=>2.5, Q=>-3.5, E=>-3.5, G=>-0.4, H=>-3.2, I=>4.5,
          L=>3.8, K=>-3.9, M=>1.9, F=>2.8, P=>-1.6, S=>-0.8, T=>-0.7, W=>-0.9, Y=>-1.3, V=>4.2);

my %EISENBERG = (A=>0.62, R=>-2.53, N=>-0.78, D=>-0.90, C=>0.29, Q=>-0.85, E=>-0.74, G=>0.48,
                 H=>-0.40, I=>1.38, L=>1.06, K=>-1.50, M=>0.64, F=>1.19, P=>0.12, S=>-0.18,
                 T=>-0.05, W=>0.81, Y=>0.26, V=>1.08);

my %BOMAN = (L=>4.92, I=>4.92, V=>4.04, F=>2.98, M=>2.35, W=>2.33, A=>1.81, C=>1.28, G=>0.94,
             Y=>-0.14, T=>-2.57, S=>-3.40, H=>-4.66, Q=>-5.54, K=>-5.55, N=>-6.64, E=>-6.81,
             D=>-8.72, R=>-14.92, P=>0.0);

sub new {
    my ($class, %o) = @_;
    my $self = bless { ph => $o{ph} // 7.0, pkset => $o{pkset} // 'EMBOSS' }, $class;
    $self->{calc} = Bio::Tools::pICalculator->new(-places => 2, -pKset => $self->{pkset});
    return $self;
}

sub columns { return @COLUMNS }

sub count {
    my ($seq, $set) = @_;
    my $n = () = $seq =~ /[$set]/g;
    return $n;
}

# Returns a hashref of properties, or undef if the sequence is not 20 standard residues.
sub compute {
    my ($self, $seq) = @_;
    $seq = uc($seq // '');
    $seq =~ s/\s+//g;
    return undef unless $seq =~ /^[ACDEFGHIKLMNPQRSTVWY]+$/;
    my $n = length $seq;
    my %p = (length => $n);

    my $bio = Bio::Seq->new(-seq => $seq, -alphabet => 'protein');
    my $mw  = Bio::Tools::SeqStats->get_mol_wt($bio);
    $p{mol_weight} = sprintf '%.2f', $mw->[0];

    $self->{calc}->seq($bio);
    $p{pI}         = sprintf '%.2f', $self->{calc}->iep;
    $p{net_charge} = sprintf '%.2f', $self->{calc}->charge_at_pH($self->{ph});

    my ($kd, $bo) = (0, 0);
    for my $aa (split //, $seq) { $kd += $KD{$aa}; $bo += $BOMAN{$aa}; }
    $p{gravy} = sprintf '%.3f', $kd / $n;
    $p{boman} = sprintf '%.3f', -$bo / $n;

    $p{hydrophobic_pct} = sprintf '%.1f', 100 * count($seq, 'ACFILMVW') / $n;
    $p{cationic_pct}    = sprintf '%.1f', 100 * count($seq, 'KR') / $n;
    $p{anionic_pct}     = sprintf '%.1f', 100 * count($seq, 'DE') / $n;

    my ($a, $v, $il) = (count($seq, 'A'), count($seq, 'V'), count($seq, 'IL'));
    $p{aliphatic_index} = sprintf '%.2f', 100 * ($a + 2.9 * $v + 3.9 * $il) / $n;

    my ($sx, $sy) = (0, 0);
    my $pi = 4 * atan2(1, 1);
    my @aa = split //, $seq;
    for my $i (0 .. $#aa) {
        my $ang = $i * 100 * $pi / 180;
        $sx += $EISENBERG{ $aa[$i] } * cos($ang);
        $sy += $EISENBERG{ $aa[$i] } * sin($ang);
    }
    $p{hydrophobic_moment} = sprintf '%.3f', sqrt($sx * $sx + $sy * $sy) / $n;

    return \%p;
}

1;
