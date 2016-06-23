#!/usr/bin/perl -w
# -*-Perl-*-
# Last changed Time-stamp: <2016-06-23 09:03:38 sven>

# input barriers .bar file and a treekin trajectory
# output coarse grained dynamic obtained from merging lmins

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use strict;
use warnings;

my $minh=3;
my $conv=1;
pod2usage(-verbose => 0)
    unless GetOptions("minh=f"  => \$minh,
		      "conv=f" =>  \$conv,
		      "man"   => sub{pod2usage(-verbobose => 2)},
		      "help"  => sub{pod2usage(-verbose => 1)});

my @regs = ("^.{25}\\({11}\\.{6}\\){11}[\\.\\(\\)]{11}|^.{26}\\({10}\\.{6}\\){10}[\\.\\(\\)]{11}", # ?25(((((((((((......))))))))))) | ?26((((((((((......))))))))))
	    "^.{2}\\({3}\\.{3}\\({8}\\.{5}\\){5}\\.{3}\\){3}\\.{3}\\){3}|^.{2}\\({3}\\.{3}\\({8}\\.{5}\\){4}\\.{3}\\){4}\\.{3}\\){3}" # ?2(((...((((((((.....)))))...)))...))) | ?2(((...((((((((.....))))...))))...)))"
	   );  # list of regexes that mark special features, only lmins with the same tags can be merged.
my @tags = ("");  # foreach lmin list which regex in @regs matched

my ($seq, @lmin) = read_bar();

foreach (@lmin[1..$#lmin]) { # [0] is empty
  my $str = $_->[1];
  my $t = "";
  foreach my $r (@regs) {
    $t .= ($str =~ /$r/) ? '1' : '0';
  }
  push @tags, $t;
}

my @sbar = sort {$b->[4] <=> $a->[4]} @lmin[1..$#lmin];

my @merge = ();
$merge[$_] = $_ for (0..$#lmin);

foreach my $s (@sbar) {
  next if $s->[4]>=$minh;
  next if $s->[0]==1; # just in case: never merge the ground state;
  my $f = $s->[3];
  next unless $tags[$s->[0]] eq $tags[$f];
  $f = $merge[$f] while ($merge[$f] != $f);
  $merge[$s->[0]] = $f;
}

my $newbar = $ARGV;
$newbar =~ s/\.bar$//; $newbar .= "_$minh" . ".bar";
open(my $fh, ">", $newbar) or die "can't open > $newbar";

my @new = print_bar($fh, $seq, \@lmin, \@merge);  # @new contains the new lmin labels
close($fh);

my $newtkin = $ARGV[-1];
$newtkin =~ s/\.tkin//; $newtkin .= "_$minh" . ".tkin";
open($fh, ">", $newtkin) or die "can't open > $newtkin";
process_treekin($fh, \@merge, $conv);
close($fh);

#---
sub read_bar {
  $_ = <>;
  warn "no seq in bar file" unless /^\s+(\S+)/;
  my $seq = $1;
  my @lmin;
  my $nn = 1;
  while (<>) {
    my @F = split;
    next if @F < 5; # shouldn't happen
    #    splice(@F,2,1) if ($F[2] eq '('); # got 2 fields from e.g. "( -1.00)"
    #    $F[2] =~ s/[()]//g;               # remove brackets from energy field
    $lmin[$nn++] = \@F;
    last if eof;
  }
  return $seq, @lmin;
}

sub print_bar {
  my ($fh, $seq, $lmin, $merge) = @_;
  print $fh "    $seq\n";
  my $n=1;
  my @new = (0);
  foreach my $l (@{$lmin}[1..$#$lmin]) {
    my ($num, $str, $e, $f, $b) = @{$l};
    #    next if $f==0; # skip unconnected
    next unless $num == $merge->[$num];
    $f = $merge->[$f] while ($merge->[$f] < $f);
    printf $fh "%3d %s %6.2f %4d %6.2f\n", $n, $str, $e, $new[$f], $b;
    $new[$num] = $n++;
  }
  return @new;
}

sub process_treekin {
  my ($fh, $merge, $conv) = @_;
  my $head;
  while (<>) {
    print $fh $_, next if /^#/;
    print $fh "# coarse grained with min height ",$minh," and regexes ",join("\t", @regs),"\n" unless $head;
    $head=1;
    my @l = split;
    my @ll = @l;

    printf $fh "%e ", ($l[0]*$conv); # print time
    for (my $i=$#l; $i>0; $i--) { # from right to left since we always have merge[l]<l;
      next unless(defined $merge->[$i]); # if I use absorbance I get more states than minima
      next if $merge->[$i] == $i;
      $l[$merge->[$i]] += $l[$i];
    }
    my ($s1,$s2) = (0,0);
    for my $i (1..$#l) {
      $s2 += $ll[$i];
      next if (defined $merge->[$i] and $merge->[$i] != $i); # if I use absorbance I get more states than minima
      printf $fh "%e ", $l[$i];
      $s1 += $l[$i];
    }
    print $fh "\n";
    # print STDERR $s1-1, " ", $s2-1, " ", 1-$s1/$s2, "\n";
  }
}
