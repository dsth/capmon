#!/homes/dsth/dev/localperl/bin/perl
use strict;
use Storable qw(store retrieve freeze thaw dclone);
use Data::Dumper;
my $gff = $ARGV[0];
$gff || die qq{\nusage: perl $0 <gff_validation_2_prob_file};
my $cap_mode = $ARGV[1] eq 'cap' ? 1 : 0;
my $lib = '/homes/dsth/dev/NewCap/backend/monitor:/net/isilon3/production/panda/ensemblgenomes/development/dsth/imports/:/homes/dsth/ensembl_src/bioperl-1.2.3:/homes/dsth/dev/EG14/ensembl/modules';

$ENV{PERL5LIB} = $lib;

print qq{\n$lib\n} if (!$cap_mode);

print qq{\nchecking file $gff\n} if (!$cap_mode);
open (my $f, '<', $gff) || die qq{problem opening temporary file};

my @h;
my %p2id;
while (my $l=<$f>) {
    if($l=~/[^\t]+\t[^\t]+\tgene\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t.*?ID=([^;]+?);/) {
        push @h, $1;
    }
    my ($id,$p) = (q{},q{});
    if ($l=~/ID=([^;]+?);/) { $id = $1; }
    if ($l=~/Parent=([^;]+?);/) { $p = $1; }
    if ($id && $p) { push @{$p2id{$p}}, $id; }
}

close($f);

my $c;
print qq{\nExhaustively checking genes...};
my $tmp = '_tmp.gff';

for my $g (@h) {

    my @listt;
    my @listl;
    print qq{\nchecking $g} if (!$cap_mode);
    push @listt, @{$p2id{$g}}; 
    for my $t (@listt) {
        push @listl, @{$p2id{$t}}; 
    }
    my %allids;
    @allids{$g,@listt,@listl} = ();
    open (my $f, '<', $gff) || die qq{couldn't open temporary file};
    open (my $of, '>', $tmp) || die qq{couldn't open temporary file};
    while (my $l=<$f>) {
        my $id = q{};
        if ($l=~/ID=([^;]+?);/) { 
            if (exists $allids{$1}) {
                print $of $l;
            }
        } elsif ($l=~/Parent=([^;]+?);/) { 
            if (exists $allids{$1}) {
                print $of $l;
            }
        }
    }

    close($f);
    close ($of);
    my $cmd = q{/homes/dsth/dev/localperl/bin/perl }
      .q{/homes/dsth/dev/NewCap/backend/monitor/GffDoc.pl }
      .q{-type contig=ignore -type match=ignore -type match_part=ignore -type pseudogenic_tRNA=mRNA:pseudogenic_tRNA }
      .q{-type ncRNA=mRNA:ncRNA   -type tRNA=mRNA:tRNA -type miRNA=mRNA:miRNA -type rRNA=mRNA:rRNA -mRNA_biotype_dominant } 
      .q{-non_protein_coding_types pseudogenic_tRNA -coordsystem scaffold -non_coding_cds -non_protein_coding -leafnonunique }
      .q{--dbname dummy --host mysql-eg-devel-3.ebi.ac.uk --user ensrw --port 4208 }
      .q{--pass scr1b3d3 --validate --file _tmp.gff};
    open (my $ph, qq{$cmd 2>&1 |}) || die;
    my $out = do { local( $/ ) ; <$ph> } ;
    if ($out !~ /> Saved 0 genes\./) {
        print qq{\n> There seems to be a problem with $g};
        if ($cap_mode) {
            unlink '_tmp.gff';
        } else { 
            rename('_tmp.gff', '_'.$g.'.gff');
            print qq{\nstored as _${g}.gff\n};
        }
        my $cap_msg;
        if ($cap_mode && $out =~ /(.*EXCEPTION.*\n.*\n.*\n.*\n.*\n.*\n)/xm) {
            $cap_msg = $1;
            $cap_msg =~ s/Exception::GffDoc/Exception::GffDoc::ReferentialIntegrity/;
            $cap_msg =~ s/Types::IllegalCombination/Undefined::CDS/;
            $cap_msg =~ s/GffDoc::Preparser/GffDoc::Gff3::mRNA/;
            $cap_msg =~ s/ - Read the above message!/\./;
            print qq{\n$cap_msg};
            print '-' x 51;
        }
    } else {
        unlink('_tmp.gff');
    }
    $c=$cmd;
}
print qq{\nPERL5LIB=$lib \\\n} if (!$cap_mode);
print substr($c,0,length($c)-8).qq{\\\n\\\n\n\n} if (!$cap_mode);
