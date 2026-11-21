package Charts;
# Charts.pm - server-side SVG charts for the marine peptide web server.
# Pure Perl (core modules only): no JavaScript library and no external files.
use strict;
use warnings;
use List::Util qw(max sum);
use POSIX qw(ceil);

our @AA = split //, 'ACDEFGHIKLMNPQRSTVWY';
our @KEYPROPS = (
    [ net_charge         => 'Net charge (pH 7)' ],
    [ pI                 => 'Isoelectric point' ],
    [ gravy              => 'GRAVY (hydrophobicity)' ],
    [ hydrophobic_moment => 'Hydrophobic moment' ],
    [ cationic_pct       => 'Cationic residues (%)' ],
    [ length             => 'Length (residues)' ],
);

my ($BLUE, $GREEN, $GREY, $GRID, $TEXT) = ('#0b5394', '#2e9e5b', '#9aa5b1', '#d9e2ec', '#334e68');
my $FONT = 'font-family="system-ui,Arial,sans-serif" font-size="12"';

sub esc {
    my $s = defined $_[0] ? $_[0] : '';
    $s =~ s/&/&amp;/g; $s =~ s/</&lt;/g; $s =~ s/>/&gt;/g; $s =~ s/"/&quot;/g;
    return $s;
}

sub shorten {
    my ($s, $n) = @_;
    $s = '' unless defined $s;
    return length($s) > $n ? substr($s, 0, $n - 3) . '...' : $s;
}

sub nice_axis {                 # -> (top, step) with a tidy step size
    my ($max, $min_top) = @_;
    $max = $min_top if $max < $min_top;
    for my $step (1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500) {
        return (ceil($max / $step) * $step, $step) if $max / $step <= 5;
    }
    return (ceil($max / 500) * 500, 500);
}

sub text {                      # text element
    my ($x, $y, $s, %o) = @_;
    my $anchor = $o{anchor} // 'start';
    my $fill   = $o{fill}   // $TEXT;
    my $extra  = $o{bold} ? ' font-weight="600"' : '';
    return sprintf('<text x="%.1f" y="%.1f" text-anchor="%s" fill="%s" %s%s>%s</text>',
        $x, $y, $anchor, $fill, $FONT, $extra, esc($s));
}

sub line {
    my ($x1, $y1, $x2, $y2, $color, $extra) = @_;
    return sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"%s/>',
        $x1, $y1, $x2, $y2, $color, $extra // '');
}

sub open_svg {
    my ($w, $h, $label) = @_;
    return qq{<svg class="chart" viewBox="0 0 $w $h" width="100%" style="max-width:${w}px;height:auto" }
         . qq{role="img" aria-label="} . esc($label) . qq{" xmlns="http://www.w3.org/2000/svg">};
}

sub legend {                    # items: [color, label, shape], starts at (x,y)
    my ($x, $y, @items) = @_;
    my $out = '';
    for my $it (@items) {
        my ($color, $label, $shape) = @$it;
        $shape //= 'box';
        if ($shape eq 'dot') {
            $out .= sprintf('<circle cx="%.1f" cy="%.1f" r="5" fill="%s"/>', $x + 5, $y - 4, $color);
        } elsif ($shape eq 'diamond') {
            $out .= sprintf('<path d="M%.1f %.1f l5 5 l-5 5 l-5 -5 z" fill="%s"/>', $x + 5, $y - 9, $color);
        } else {
            $out .= sprintf('<rect x="%.1f" y="%.1f" width="11" height="11" fill="%s"/>', $x, $y - 10, $color);
        }
        $out .= text($x + 16, $y, $label);
        $x += 22 + 7.2 * length($label);
    }
    return $out;
}

# ---- reference statistics from the training set ----------------------------
sub load_reference {
    my ($file) = @_;
    open my $fh, '<', $file or return undef;
    chomp(my $head = <$fh>);
    my @cols = split /\t/, $head;
    my (%n, %aac, %vals);
    while (<$fh>) {
        chomp;
        next unless length;
        my %r;
        @r{@cols} = split /\t/, $_, -1;
        my $cl = $r{label} // '';
        next unless $cl eq 'ACP' || $cl eq 'nonACP';
        my $seq = $r{sequence} // '';
        my $len = length $seq;
        next unless $len;
        $n{$cl}++;
        my %cnt;
        $cnt{$_}++ for split //, $seq;
        $aac{$cl}{$_} += ($cnt{$_} // 0) / $len for @AA;
        for my $kp (@KEYPROPS) {
            my $v = $r{ $kp->[0] };
            push @{ $vals{ $kp->[0] }{$cl} }, $v if defined $v && $v =~ /^-?\d+(?:\.\d+)?$/;
        }
    }
    close $fh;
    return undef unless $n{ACP} && $n{nonACP};

    my %ref;
    for my $cl (qw(ACP nonACP)) {
        $ref{aac}{$cl}{$_} = $aac{$cl}{$_} / $n{$cl} for @AA;
    }
    for my $kp (@KEYPROPS) {
        my $p = $kp->[0];
        my @all = sort { $a <=> $b } map { @{ $vals{$p}{$_} || [] } } qw(ACP nonACP);
        next unless @all;
        my %m;
        for my $cl (qw(ACP nonACP)) {
            my @v = @{ $vals{$p}{$cl} || [] };
            $m{$cl} = @v ? sum(@v) / @v : undef;
        }
        $ref{prop}{$p} = { lo => $all[ int(0.05 * $#all) ], hi => $all[ int(0.95 * $#all) ],
                           ACP => $m{ACP}, nonACP => $m{nonACP}, n => { %n } };
    }
    $ref{n} = { %n };
    return \%ref;
}

# ---- chart 1: probability per submitted peptide ----------------------------
sub prob_chart {
    my ($res) = @_;             # [ { name, prob, call } ]
    my $n = @$res or return '';
    my ($W, $bx, $bw, $rowh, $top) = (760, 230, 440, 28, 30);
    my $H = $top + $n * $rowh + 30;
    my $svg = open_svg($W, $H, 'Anticancer probability for each peptide');
    for my $t (0, 0.25, 0.5, 0.75, 1) {
        my $x = $bx + $t * $bw;
        $svg .= line($x, $top - 6, $x, $top + $n * $rowh, $GRID);
        $svg .= text($x, $top + $n * $rowh + 16, sprintf('%.2f', $t), anchor => 'middle');
    }
    my $x5 = $bx + 0.5 * $bw;
    $svg .= line($x5, $top - 6, $x5, $top + $n * $rowh, '#b42318', ' stroke-dasharray="5 4" stroke-width="1.5"');
    $svg .= text($x5, $top - 12, 'threshold 0.5', anchor => 'middle', fill => '#b42318');
    my $i = 0;
    for my $r (@$res) {
        my $y = $top + $i * $rowh;
        my $p = $r->{prob} + 0;
        my $color = $p >= 0.5 ? $GREEN : $GREY;
        $svg .= text($bx - 10, $y + 18, shorten($r->{name}, 30), anchor => 'end');
        $svg .= sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="18" fill="%s" rx="2"/>',
            $bx, $y + 5, $p * $bw, $color);
        $svg .= text($bx + $p * $bw + 6, $y + 19, sprintf('%.3f', $p), bold => 1);
        $i++;
    }
    $svg .= legend($bx, $H - 2 - 14 + 14, [ $GREEN, 'predicted ACP (p >= 0.5)' ], [ $GREY, 'predicted non-ACP' ]);
    return $svg . '</svg>';
}

# ---- chart 2: amino-acid composition vs. training classes ------------------
sub aac_chart {
    my ($seq, $ref) = @_;
    return '' unless $ref && defined $seq && length $seq;
    my $len = length $seq;
    my %cnt;
    $cnt{$_}++ for split //, $seq;
    my @q = map { 100 * ($cnt{$_} // 0) / $len } @AA;
    my @a = map { 100 * $ref->{aac}{ACP}{$_} } @AA;
    my @b = map { 100 * $ref->{aac}{nonACP}{$_} } @AA;
    my ($ytop, $ystep) = nice_axis(max(@q, @a, @b), 10);

    my ($W, $H, $px0, $px1, $py0, $py1) = (760, 300, 50, 750, 44, 250);
    my $gw = ($px1 - $px0) / @AA;
    my $bwid = 9;
    my $svg = open_svg($W, $H, 'Amino acid composition of the peptide compared with ACP and non-ACP training peptides');
    $svg .= legend(50, 20, [ $BLUE, 'this peptide' ], [ $GREEN, 'mean of training ACPs' ], [ $GREY, 'mean of training non-ACPs' ]);
    for (my $v = 0; $v <= $ytop + 1e-9; $v += $ystep) {
        my $y = $py1 - ($py1 - $py0) * $v / $ytop;
        $svg .= line($px0, $y, $px1, $y, $GRID);
        $svg .= text($px0 - 6, $y + 4, sprintf('%g', $v), anchor => 'end');
    }
    $svg .= sprintf('<text transform="translate(13 %.1f) rotate(-90)" text-anchor="middle" fill="%s" %s>%% of residues</text>',
        ($py0 + $py1) / 2, $TEXT, $FONT);
    for my $i (0 .. $#AA) {
        my $gx = $px0 + $i * $gw + ($gw - 3 * $bwid - 2) / 2;
        my $j = 0;
        for my $pair ([ $q[$i], $BLUE ], [ $a[$i], $GREEN ], [ $b[$i], $GREY ]) {
            my $h = ($py1 - $py0) * $pair->[0] / $ytop;
            $svg .= sprintf('<rect x="%.1f" y="%.1f" width="%d" height="%.1f" fill="%s"/>',
                $gx + $j * ($bwid + 1), $py1 - $h, $bwid, $h, $pair->[1]);
            $j++;
        }
        $svg .= text($px0 + ($i + 0.5) * $gw, $py1 + 16, $AA[$i], anchor => 'middle', bold => 1);
    }
    $svg .= text(($px0 + $px1) / 2, $H - 12, 'amino acid (one-letter code)', anchor => 'middle');
    return $svg . '</svg>';
}

# ---- chart 3: histogram of marine scores with the query marked --------------
sub hist_chart {
    my ($probs, $marks) = @_;   # \@probabilities, [ [label, prob], ... ]
    return '' unless $probs && @$probs;
    my $nb = 20;
    my @c = (0) x $nb;
    for my $p (@$probs) {
        my $b = int($p * $nb);
        $b = $nb - 1 if $b >= $nb;
        $b = 0 if $b < 0;
        $c[$b]++;
    }
    my ($ctop, $cstep) = nice_axis(max(@c), 5);
    my ($W, $H, $px0, $px1, $py0, $py1) = (760, 300, 56, 740, 52, 250);
    my $bw = ($px1 - $px0) / $nb;
    my $svg = open_svg($W, $H, 'Distribution of anticancer probabilities across the marine peptides');
    for (my $v = 0; $v <= $ctop + 1e-9; $v += $cstep) {
        my $y = $py1 - ($py1 - $py0) * $v / $ctop;
        $svg .= line($px0, $y, $px1, $y, $GRID);
        $svg .= text($px0 - 6, $y + 4, sprintf('%g', $v), anchor => 'end');
    }
    $svg .= sprintf('<text transform="translate(14 %.1f) rotate(-90)" text-anchor="middle" fill="%s" %s>number of marine peptides</text>',
        ($py0 + $py1) / 2, $TEXT, $FONT);
    for my $i (0 .. $nb - 1) {
        my $h = ($py1 - $py0) * $c[$i] / $ctop;
        $svg .= sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>',
            $px0 + $i * $bw + 1, $py1 - $h, $bw - 2, $h, ($i + 0.5) / $nb >= 0.5 ? $GREEN : $GREY);
    }
    for my $t (0, 0.25, 0.5, 0.75, 1) {
        $svg .= text($px0 + $t * ($px1 - $px0), $py1 + 16, sprintf('%.2f', $t), anchor => 'middle');
    }
    my $x5 = $px0 + 0.5 * ($px1 - $px0);
    $svg .= line($x5, $py0 - 4, $x5, $py1, '#b42318', ' stroke-dasharray="5 4" stroke-width="1.5"');
    my $i = 0;
    for my $m (@{ $marks || [] }) {
        my $x = $px0 + $m->[1] * ($px1 - $px0);
        my $ty = $py0 - 8 - 15 * ($i % 3);
        $svg .= line($x, $ty + 3, $x, $py1, $BLUE, ' stroke-width="2"');
        $svg .= sprintf('<circle cx="%.1f" cy="%.1f" r="4" fill="%s"/>', $x, $ty, $BLUE);
        my $anchor = $x > $px1 - 140 ? 'end' : 'start';
        $svg .= text($x + ($anchor eq 'end' ? -7 : 7), $ty + 4, shorten($m->[0], 22) . sprintf(' (%.2f)', $m->[1]), anchor => $anchor, fill => $BLUE, bold => 1);
        $i++;
    }
    $svg .= text(($px0 + $px1) / 2, $H - 10, 'predicted anticancer probability', anchor => 'middle');
    $svg .= legend($px0, 16, [ $GREEN, 'p >= 0.5 (ACP)' ], [ $GREY, 'p < 0.5' ]);
    return $svg . '</svg>';
}

# ---- chart 4: where the peptide sits on key properties ----------------------
sub prop_strips {
    my ($row, $ref) = @_;
    return '' unless $ref && $row;
    my @items = grep { $ref->{prop}{ $_->[0] } && defined $row->{ $_->[0] } && $row->{ $_->[0] } =~ /^-?\d+(?:\.\d+)?$/ } @KEYPROPS;
    return '' unless @items;
    my ($W, $x1, $x2, $rowh, $top) = (760, 250, 640, 44, 44);
    my $H = $top + @items * $rowh + 8;
    my $svg = open_svg($W, $H, 'Position of the peptide on key physicochemical properties');
    $svg .= legend(10, 20, [ $BLUE, 'this peptide', 'dot' ], [ $GREEN, 'training ACP mean', 'diamond' ], [ $GREY, 'training non-ACP mean', 'diamond' ]);
    my $i = 0;
    for my $it (@items) {
        my ($p, $label) = @$it;
        my $st = $ref->{prop}{$p};
        my ($lo, $hi) = ($st->{lo}, $st->{hi});
        $hi = $lo + 1 if $hi <= $lo;
        my $y = $top + $i * $rowh + 18;
        my $sx = sub { my $v = shift; $v = $lo if $v < $lo; $v = $hi if $v > $hi; return $x1 + ($v - $lo) / ($hi - $lo) * ($x2 - $x1); };
        $svg .= text(10, $y + 4, $label);
        $svg .= line($x1, $y, $x2, $y, $GRID, ' stroke-width="6" stroke-linecap="round"');
        $svg .= text($x1, $y + 20, sprintf('%g', $lo), anchor => 'middle');
        $svg .= text($x2, $y + 20, sprintf('%g', $hi), anchor => 'middle');
        for my $pair ([ nonACP => $GREY ], [ ACP => $GREEN ]) {
            next unless defined $st->{ $pair->[0] };
            my $x = $sx->($st->{ $pair->[0] });
            $svg .= sprintf('<path d="M%.1f %.1f l6 6 l-6 6 l-6 -6 z" fill="%s"/>', $x, $y - 6, $pair->[1]);
        }
        my $v = $row->{$p} + 0;
        my $out = ($v < $lo || $v > $hi) ? ' stroke="#b42318" stroke-width="2"' : ' stroke="#ffffff" stroke-width="1.5"';
        $svg .= sprintf('<circle cx="%.1f" cy="%.1f" r="7" fill="%s"%s/>', $sx->($v), $y, $BLUE, $out);
        $svg .= text($x2 + 22, $y + 4, sprintf('%g', $v), bold => 1, fill => $BLUE);
        $i++;
    }
    return $svg . '</svg>';
}

# ---- everything the predict page needs, in one call --------------------------
sub predict_charts {
    my ($res, $probs, $ref, $limit) = @_;
    $limit //= 10;
    return ('', '') unless $res && @$res;
    my @few = @$res[ 0 .. ($#$res < $limit - 1 ? $#$res : $limit - 1) ];
    for my $r (@few) {
        $r->{aac_svg}  = aac_chart($r->{sequence}, $ref);
        $r->{prop_svg} = prop_strips($r, $ref);
    }
    return (prob_chart($res), hist_chart($probs, [ map { [ $_->{name}, $_->{prob} ] } @few ]));
}

1;
