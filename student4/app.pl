#!/usr/bin/env perl
# app.pl - web interface for the marine bioactive peptide database (Mojolicious::Lite).
#
# Run:  perl student4/app.pl daemon -l http://127.0.0.1:3000
# Open: http://127.0.0.1:3000
use Mojolicious::Lite -signatures;
use DBI;
use FindBin qw($Bin);
use lib $Bin;
use Charts;

my $dbfile = $ENV{MARINE_DB} // "$Bin/../data/final/marine.db";
die "Database not found: $dbfile (run build_db.pl first)\n" unless -e $dbfile;

# reference statistics of the training set, used by the composition and property charts
my $REF = Charts::load_reference($ENV{MARINE_TRAIN} // "$Bin/../data/final/training_balanced_props.tsv");

helper marine_probs => sub ($c) {
    return $c->db->selectcol_arrayref('SELECT acp_probability FROM peptides');
};

helper db => sub {
    state $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '',
        { RaiseError => 1, sqlite_unicode => 1 });
    return $dbh;
};

# ---- home: summary ---------------------------------------------------------
get '/' => sub ($c) {
    my $db = $c->db;
    my ($total, $acp) = $db->selectrow_array(
        q{SELECT COUNT(*), SUM(acp_call = 'ACP') FROM peptides});
    my $groups = $db->selectall_arrayref(
        q{SELECT taxon_group, COUNT(*) n, SUM(acp_call = 'ACP') acp, ROUND(AVG(acp_probability),3) avg_p
          FROM peptides GROUP BY taxon_group ORDER BY n DESC}, { Slice => {} });
    my $top = $db->selectall_arrayref(
        q{SELECT pep_id, sequence, organism, taxon_group, acp_probability
          FROM peptides ORDER BY acp_probability DESC LIMIT 10}, { Slice => {} });
    $c->render(template => 'home', total => $total, acp => $acp, groups => $groups, top => $top,
        hist_svg => Charts::hist_chart($c->marine_probs, []));
};

# ---- browse / search -------------------------------------------------------
get '/peptides' => sub ($c) {
    my $db    = $c->db;
    my $q     = $c->param('q')     // '';
    my $group = $c->param('group') // '';
    my $minp  = $c->param('min_prob') // '';
    my $sort  = $c->param('sort')  // 'acp_probability';
    my $dir   = ($c->param('dir') // 'desc') eq 'asc' ? 'ASC' : 'DESC';
    my $page  = int($c->param('page') // 1);
    $page = 1 if $page < 1;
    my $per = 25;

    my %ok = map { $_ => 1 } qw(acp_probability length net_charge pI gravy mol_weight acp_rank pep_id);
    $sort = 'acp_probability' unless $ok{$sort};

    my (@where, @bind);
    if (length $q) {
        push @where, '(sequence LIKE ? OR organism LIKE ? OR pep_id LIKE ?)';
        push @bind, ("%$q%") x 3;
    }
    if (length $group) { push @where, 'taxon_group = ?'; push @bind, $group; }
    if (length $minp && $minp =~ /^\d*\.?\d+$/) { push @where, 'acp_probability >= ?'; push @bind, $minp; }
    my $w = @where ? 'WHERE ' . join(' AND ', @where) : '';

    my ($count) = $db->selectrow_array("SELECT COUNT(*) FROM peptides $w", undef, @bind);
    my $rows = $db->selectall_arrayref(
        "SELECT pep_id, sequence, length, organism, taxon_group, net_charge, pI, acp_probability, acp_call, is_benchmark
         FROM peptides $w ORDER BY $sort $dir LIMIT $per OFFSET " . (($page - 1) * $per),
        { Slice => {} }, @bind);
    my $groups = $db->selectcol_arrayref(
        q{SELECT DISTINCT taxon_group FROM peptides ORDER BY taxon_group});
    my $pages = int(($count + $per - 1) / $per) || 1;

    $c->render(template => 'peptides', rows => $rows, count => $count, groups => $groups,
        q => $q, group => $group, minp => $minp, sort => $sort, dir => $dir,
        page => $page, pages => $pages);
};

# ---- one peptide -----------------------------------------------------------
get '/peptide/:id' => sub ($c) {
    my $row = $c->db->selectrow_hashref('SELECT * FROM peptides WHERE pep_id = ?', undef, $c->stash('id'));
    return $c->reply->not_found unless $row;
    $c->render(template => 'peptide', p => $row,
        aac_svg  => Charts::aac_chart($row->{sequence}, $REF),
        prop_svg => Charts::prop_strips($row, $REF),
        hist_svg => Charts::hist_chart($c->marine_probs, [ [ $row->{pep_id}, $row->{acp_probability} ] ]));
};

# ---- JSON API --------------------------------------------------------------
get '/api/peptide/:id' => sub ($c) {
    my $row = $c->db->selectrow_hashref('SELECT * FROM peptides WHERE pep_id = ?', undef, $c->stash('id'));
    return $c->render(json => { error => 'not found' }, status => 404) unless $row;
    $c->render(json => $row);
};

# ---- predict: paste your own sequences -------------------------------------
use File::Temp qw(tempdir);

sub capture {                  # run a command (no shell), return its stdout
    my @cmd = @_;
    open(my $fh, '-|', @cmd) or die "cannot run $cmd[0]: $!\n";
    local $/;
    my $out = <$fh> // '';
    close $fh;
    return ($out, $?);
}

sub read_table {               # TSV text -> list of hashes
    my ($text) = @_;
    my @lines = grep { length } split /\n/, $text;
    my @cols  = split /\t/, shift @lines;
    return map { my %r; @r{@cols} = split /\t/, $_, -1; \%r } @lines;
}

sub run_predict {              # returns (\@results, $error)
    my ($input) = @_;
    my (@q, @bad);
    if ($input =~ /^\s*>/) {                       # FASTA
        for my $rec (split /^>/m, $input) {
            next unless $rec =~ /\S/;
            my ($head, @body) = split /\n/, $rec;
            push @q, [ ($head // '') =~ s/^\s+|\s+$//gr, join('', @body) ];
        }
    } else {                                       # one sequence per line
        my $i = 0;
        for my $line (split /\n/, $input) {
            next unless $line =~ /\S/;
            $i++;
            push @q, [ "sequence $i", $line ];
        }
    }
    return (undef, 'Please paste at least one sequence.') unless @q;
    return (undef, 'Please submit at most 50 sequences at a time.') if @q > 50;

    my @ok;
    for my $x (@q) {
        (my $s = uc $x->[1]) =~ s/\s+//g;
        my $name = length $x->[0] ? substr($x->[0], 0, 60) : 'unnamed';
        if ($s !~ /^[ACDEFGHIKLMNPQRSTVWY]+$/) { push @bad, "$name: only the 20 standard amino acid letters are allowed"; next; }
        if (length($s) < 5 || length($s) > 50)  { push @bad, "$name: length " . length($s) . " is outside the 5-50 residue range of the model"; next; }
        push @ok, [ $name, $s ];
    }
    return (undef, join('; ', @bad)) unless @ok;

    my $dir = tempdir(CLEANUP => 1);
    my $fa  = "$dir/in.fasta";
    open(my $o, '>', $fa) or die "tmp: $!";
    my $n = 0;
    print $o ">Q" . ++$n . "|query\n$_->[1]\n" for @ok;
    close $o;

    my ($props, $st1) = capture($^X, "$Bin/../student2/compute_props.pl", $fa);
    return (undef, 'Property calculation failed.') if $st1 != 0 || $props !~ /\S/;
    open(my $p, '>', "$dir/props.tsv") or die "tmp: $!";
    print $p $props;
    close $p;

    my $model = $ENV{MARINE_MODEL} // "$Bin/../data/final/model/model_et.joblib";
    my ($sc, $st2) = capture('python3', "$Bin/../student3/predict_et.py", $model, "$dir/props.tsv");
    return (undef, 'Scoring failed (is the model file present?).') if $st2 != 0 || $sc !~ /\S/;

    my %prop  = map { $_->{pep_id} => $_ } read_table($props);
    my %score = map { $_->{pep_id} => $_ } read_table($sc);
    my @res;
    $n = 0;
    for my $x (@ok) {
        my $id = 'Q' . ++$n;
        next unless $prop{$id} && $score{$id};
        push @res, { name => $x->[0], %{ $prop{$id} },
                     prob => $score{$id}{anticancer_probability}, call => $score{$id}{call} };
    }
    return (\@res, @bad ? 'Skipped: ' . join('; ', @bad) : undef);
}

get '/predict' => sub ($c) {
    $c->render(template => 'predict', results => undef, error => undef, input => '',
        prob_svg => '', hist_svg => '');
};

post '/predict' => sub ($c) {
    my $input = $c->param('seqs') // '';
    my ($res, $err) = run_predict($input);
    my ($prob_svg, $hist_svg) = Charts::predict_charts($res, $c->marine_probs, $REF);
    $c->render(template => 'predict', results => $res, error => $err, input => $input,
        prob_svg => $prob_svg, hist_svg => $hist_svg);
};

app->start;

__DATA__

@@ layouts/default.html.ep
<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title><%= title %> - Marine Peptide Server</title>
<style>
  body { font-family: system-ui, Arial, sans-serif; margin: 0; color: #1f2933; background: #fff; }
  header { border-bottom: 1px solid #d9e2ec; padding: 14px 24px; display: flex; gap: 24px; align-items: baseline; flex-wrap: wrap; }
  header b { font-size: 18px; color: #0b5394; }
  header a { color: #334e68; text-decoration: none; }
  main { max-width: 1100px; margin: 0 auto; padding: 24px; }
  table { border-collapse: collapse; width: 100%; font-size: 14px; }
  th, td { text-align: left; padding: 6px 10px; border-bottom: 1px solid #e4e7eb; vertical-align: top; }
  th { background: #f5f7fa; }
  .seq { font-family: ui-monospace, Consolas, monospace; word-break: break-all; }
  .acp { color: #0b7a3b; font-weight: 600; }
  .muted { color: #829ab1; }
  .bar { background: #e4e7eb; height: 10px; width: 140px; display: inline-block; vertical-align: middle; }
  .bar i { display: block; height: 10px; background: #0b5394; }
  form.f { display: flex; gap: 10px; flex-wrap: wrap; margin-bottom: 16px; }
  input, select, button { padding: 6px 8px; font-size: 14px; }
  .cards { display: flex; gap: 16px; flex-wrap: wrap; margin: 16px 0; }
  .card { border: 1px solid #d9e2ec; border-radius: 6px; padding: 12px 18px; }
  .card b { font-size: 22px; display: block; }
  .pg a { margin-right: 10px; }
  .chart { display: block; margin: 8px 0 20px; }
  details { margin: 14px 0; border-top: 1px solid #e4e7eb; padding-top: 8px; }
  summary { cursor: pointer; }
</style></head>
<body>
<header><b>Marine Peptide Server</b>
  <a href="/">Home</a> <a href="/peptides">Browse peptides</a> <a href="/predict">Predict</a></header>
<main><%= content %></main>
</body></html>

@@ home.html.ep
% layout 'default'; title 'Home';
<h2>Marine bioactive peptides - anticancer predictions</h2>
<p>Candidate peptides from marine organisms (UniProt), with BioPerl physicochemical properties
and an Extra Trees anticancer score.</p>
<div class="cards">
  <div class="card"><b><%= $total %></b>marine peptides</div>
  <div class="card"><b><%= $acp %></b>predicted anticancer (p &ge; 0.5)</div>
  <div class="card"><b><%= scalar @$groups %></b>taxonomic groups</div>
</div>
<h3>Score distribution of all marine peptides</h3>
<p class="muted">Each bar counts marine peptides in a probability range. Green bars are predicted anticancer (p &ge; 0.5).</p>
<%== $hist_svg %>
<h3>By taxonomic group</h3>
<table><tr><th>Group</th><th>Peptides</th><th>Predicted ACP</th><th>Mean probability</th></tr>
% for my $g (@$groups) {
<tr><td><a href="/peptides?group=<%= $g->{taxon_group} %>"><%= $g->{taxon_group} %></a></td>
<td><%= $g->{n} %></td><td><%= $g->{acp} %></td><td><%= $g->{avg_p} %></td></tr>
% }
</table>
<h3>Top 10 candidates</h3>
<table><tr><th>ID</th><th>Sequence</th><th>Organism</th><th>Probability</th></tr>
% for my $r (@$top) {
<tr><td><a href="/peptide/<%= $r->{pep_id} %>"><%= $r->{pep_id} %></a></td>
<td class="seq"><%= $r->{sequence} %></td><td><%= $r->{organism} %></td>
<td><%= sprintf('%.3f', $r->{acp_probability}) %></td></tr>
% }
</table>

@@ peptides.html.ep
% layout 'default'; title 'Browse';
<h2>Browse peptides</h2>
<form class="f" method="get" action="/peptides">
  <input name="q" value="<%= $q %>" placeholder="sequence, organism or ID">
  <select name="group"><option value="">All groups</option>
% for my $g (@$groups) {
    <option value="<%= $g %>"<%= $g eq $group ? ' selected' : '' %>><%= $g %></option>
% }
  </select>
  <input name="min_prob" value="<%= $minp %>" placeholder="min probability" size="14">
  <select name="sort">
% for my $s (qw(acp_probability length net_charge pI gravy mol_weight)) {
    <option value="<%= $s %>"<%= $s eq $sort ? ' selected' : '' %>>sort: <%= $s %></option>
% }
  </select>
  <select name="dir"><option value="desc"<%= $dir eq 'DESC' ? ' selected' : '' %>>high to low</option>
    <option value="asc"<%= $dir eq 'ASC' ? ' selected' : '' %>>low to high</option></select>
  <button>Search</button>
</form>
<p class="muted"><%= $count %> peptides - page <%= $page %> of <%= $pages %></p>
<table><tr><th>ID</th><th>Sequence</th><th>Len</th><th>Organism</th><th>Group</th><th>Charge</th><th>Probability</th></tr>
% for my $r (@$rows) {
<tr><td><a href="/peptide/<%= $r->{pep_id} %>"><%= $r->{pep_id} %></a><%= $r->{is_benchmark} ? ' *' : '' %></td>
<td class="seq"><%= $r->{sequence} %></td><td><%= $r->{length} %></td>
<td><%= $r->{organism} %></td><td><%= $r->{taxon_group} %></td><td><%= $r->{net_charge} %></td>
<td class="<%= $r->{acp_call} eq 'ACP' ? 'acp' : '' %>"><%= sprintf('%.3f', $r->{acp_probability}) %></td></tr>
% }
</table>
<p class="pg">
% if ($page > 1) {
<a href="<%= url_with->query([page => $page - 1]) %>">&laquo; previous</a>
% }
% if ($page < $pages) {
<a href="<%= url_with->query([page => $page + 1]) %>">next &raquo;</a>
% }
</p>
<p class="muted">* known marine anticancer peptide (benchmark set)</p>

@@ peptide.html.ep
% layout 'default'; title $p->{pep_id};
<h2><%= $p->{pep_id} %> <span class="muted">- <%= $p->{taxon_group} %></span></h2>
<p><b>Organism:</b> <%= $p->{organism} %><br>
<b>Source:</b> <%= $p->{sources} %> <%= $p->{source_ids} %></p>
<p class="seq" style="font-size:18px"><%= $p->{sequence} %></p>
<p><b>Anticancer probability:</b> <%= sprintf('%.3f', $p->{acp_probability}) %>
 <span class="bar"><i style="width:<%= int($p->{acp_probability} * 100) %>%"></i></span>
 <span class="<%= $p->{acp_call} eq 'ACP' ? 'acp' : '' %>"><%= $p->{acp_call} %></span>
 - rank <%= $p->{acp_rank} %> of 593</p>
% if ($p->{is_benchmark}) {
<p class="muted">This is a known marine anticancer peptide. Its score comes from a model that
saw it in training, so it is not a blind prediction.</p>
% }
<h3>Physicochemical properties (BioPerl)</h3>
<table>
% for my $k (qw(length mol_weight pI net_charge gravy hydrophobic_pct cationic_pct anionic_pct boman aliphatic_index hydrophobic_moment)) {
<tr><th><%= $k %></th><td><%= $p->{$k} %></td></tr>
% }
</table>
<h3>Charts</h3>
<h4>Score compared with all marine peptides</h4>
<%== $hist_svg %>
<h4>Amino acid composition</h4>
<%== $aac_svg %>
<h4>Key properties compared with training peptides</h4>
<%== $prop_svg %>
<p><a href="/api/peptide/<%= $p->{pep_id} %>">JSON</a></p>

@@ predict.html.ep
% layout 'default'; title 'Predict';
<h2>Predict anticancer activity</h2>
<p>Paste one peptide per line, or several in FASTA format (5-50 standard amino acids, up to 50 peptides).
Properties are calculated with BioPerl and scored with the Extra Trees model.</p>
<form method="post" action="/predict">
  <textarea name="seqs" rows="7" style="width:100%;font-family:monospace" placeholder="GIGKFLHSAKKFGKAFVGEIMNS"><%= $input %></textarea>
  <p><button>Predict</button></p>
</form>
% if ($error) {
<p style="color:#b42318"><%= $error %></p>
% }
% if ($results && @$results) {
<table><tr><th>Name</th><th>Sequence</th><th>Len</th><th>MW</th><th>pI</th><th>Charge</th><th>GRAVY</th><th>Probability</th><th>Call</th></tr>
%   for my $r (@$results) {
<tr><td><%= $r->{name} %></td><td class="seq"><%= $r->{sequence} %></td><td><%= $r->{length} %></td>
<td><%= $r->{mol_weight} %></td><td><%= $r->{pI} %></td><td><%= $r->{net_charge} %></td><td><%= $r->{gravy} %></td>
<td><%= sprintf('%.3f', $r->{prob}) %> <span class="bar"><i style="width:<%= int($r->{prob} * 100) %>%"></i></span></td>
<td class="<%= $r->{call} eq 'ACP' ? 'acp' : '' %>"><%= $r->{call} %></td></tr>
%   }
</table>
<h3>Charts</h3>
<h4>Anticancer probability of each peptide</h4>
<%== $prob_svg %>
<h4>Where your peptides fall among the marine peptides</h4>
<p class="muted">Histogram of the scores of all 593 database peptides; blue lines mark your peptides.</p>
<%== $hist_svg %>
% my $shown = 0;
% for my $r (@$results) {
%   next unless $r->{aac_svg};
<details<%= $shown++ ? '' : ' open' %>><summary><b><%= $r->{name} %></b> - composition and property profile</summary>
<h4>Amino acid composition</h4>
<%== $r->{aac_svg} %>
<h4>Key properties compared with training peptides</h4>
<%== $r->{prop_svg} %>
</details>
% }
<p class="muted">The probability is a model score from sequence composition, not a measured activity.
Validation AUC is about 0.76, so use it to prioritise peptides for lab testing.</p>
% }
