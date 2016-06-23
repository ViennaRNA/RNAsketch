#!/usr/bin/perl -w
# -*-Perl-*-
# Last changed Time-stamp: <2016-06-20 15:12:27 sven>
# $Id$


use Getopt::Long;
use RNA;
#adjust modules folder
use lib "/home/mescalin/sven/bin/MODULES";
#load used modules
#use modulename

#adapt to each program
my $usage= << "JUS";
  usage:   
  options: -d    print debug
  results: 
  example: 
JUS

#adapt to used variable(s)
use strict;
#version 1 global variable(s)
use vars qw/$opt_debug $opt_bar $opt_treekin $opt_seq $opt_k $opt_s1 $opt_s2/;
#version 2 local (only in the programm known) variable(s)

#a=s a must be String
#b=i b must be Integer
#d   d = 1 =>  if 'd' is defined, esle d = undef 
usage() unless GetOptions("d|debug" => \$opt_debug,
                          "s1|structure1=s" => \$opt_s1,
			  "s2|structure2=s" => \$opt_s2,
			  "seq|sequence=s" => \$opt_seq,
			  "bi|barmap=s" => \$opt_bar,
			  "ti|treekin=s" => \$opt_treekin,
			  "k|k-neighbor=s" => \$opt_k
);

#debug message will be printet on standard-error if opt_debug is defined
if($opt_debug){
  print STDERR $usage;
}
$opt_k = 1 unless($opt_k);
$opt_s1 = "xxxxxxxxxxxxxxxxxxxxxxxxx\\(\\(\\(\\(\\(\\(\\(\\(\\(\\(\\(......\\)\\)\\)\\)\\)\\)\\)\\)\\)\\)\\)xxxxxxxxxxx" unless($opt_s1);
$opt_s2 = "\\(\\(\\(\\(\\(\\.\\.\\.\\(\\(\\(\\(\\(\\(\\(\\(.....\\)\\)\\)\\)\\)\\.\\.\\.\\)\\)\\)\\.\\.\\.\\)\\)\\)\\)\\)..xxxxxxxxxxxxxxxxxxxxxx" unless($opt_s2);

&RNAboooooor($opt_s1);
&RNAboooooor($opt_s2);

sub RNAboooooor{
  my $input = shift(@_);
  my ($subSeq, $subStr) = ();
  my %kNeighbors;
  if($input =~ /^([x]*)([\\\(\)\.]+)([x]*)$/){
    print length($1),"\t",$2,"\t",length($3),"\n";
    my $start = length($1);
    $subStr = $2;
    $subStr =~ s/\\//gi;
    $subSeq = substr($opt_seq, $start, length($subStr));
    print $subSeq,"\n",$subStr,"\n";

    my $bef = "echo ".$subSeq." | RNAfold";
    my @a = split(/\n/, `$bef`);
    my ($str,$e)=split(/\s+/, $a[1]);
    $e =~ s/[\(\)]//gi;
    
    $bef = "echo ".$subSeq." | RNAsubopt -e ".($e*-1+5);
    @a = split(/\n/, `$bef`);
    foreach my $l (@a){
      my @line = split(/\s+/, $l);
      if($line[0] =~ /^[ACGTU]+$/){
	next;
      }
      elsif($line[0] =~ /^[\(\.\)]+$/){
	my ($inboth,$dist) = &calculate_bp_distance($line[0],$subStr);
	push(@{$kNeighbors{$dist}}, \@line) if($dist<=$opt_k);
      }
      else{
	print STDERR "ignored subopt line\n>",$l;
      }
    }
  }

  foreach my $k (sort{$a <=> $b}keys %kNeighbors){
    print "k = ",$k,"\t #neighbors = ",scalar(@{$kNeighbors{$k}}),"\n";
    foreach my $N (sort {$a->[1] <=> $b->[1]} @{$kNeighbors{$k}}){
      print $N->[0],"\t",$N->[1],"\n";
    }
  }
}

# Compare two structures of unequal length and determine the base pair distance 
# d_BP(P_a,P_b) = |BP_a| + |BP_b| - 2|BP_a \cap  BP_b| Note that we take only the 
# sub string of the longer sequence equal to the length of the shorter one.
# :param s1: First structure to be compared
# :param s2: Second structure to be compared
sub calculate_bp_distance{
    my ($s1, $s2) = @_;
    #re-order and make s1 the shorter string
    if(length($s1)>length($s2)){
      my $t=$s1;
      $s1=$s2;
      $s2=$s1;
    }

    my @bpt1=@{&create_bp_table($s1)};
    my @bpt2=@{&create_bp_table($s2)};
    my $inboth=0;
    for(my $i=0; $i<scalar(@bpt1); $i++){
      $inboth++ if($bpt1[$i] != -1 and $bpt1[$i] == $bpt2[$i]);
    }
    my $dist=($s1 =~ tr/\)//) + (substr($s2, 0, length($s1)) =~ tr/\)//) - 2*$inboth;
    return ($inboth, $dist);
}

# Takes a structure in dot bracket notation and returns a base pair table.
# Unpaired positions are -1, otherwise the index of the adjacent bracket is listed
# :param structure: string with dot-bracket notation of the strcture
# :return bpt: base pair table
sub create_bp_table(){
  my $structure=shift(@_);
  my @bpo=();
  my @bpt = (-1) x length($structure);
  for(my $i=0; $i<length($structure); $i++){
    my $substr=substr($structure, $i, 1);
    if($substr eq "("){
      push(@bpo, $i);
    }
    elsif($substr eq ")"){
      if(scalar(@bpo)){
	$bpt[pop(@bpo)] = $i;
      }
      else{
	print STDERR "Unbalanced brackets: too few opening brackets";
	exit;
      }
    }
  }
  if(scalar(@bpo)){
    print STDERR "Unbalanced brackets: too few closing brackets";
    exit;
  }
  return \@bpt;
}


exit;

$opt_s1=~tr/\x/./;
$opt_s2=~tr/\x/./;

my $sequence="";
my %states=();
open(FI, $opt_bar) or die "can't open ",$opt_bar,"\n";
while(<FI>){
  my @line = split(/\s+/, $_);
  shift(@line) if($line[0] eq ""); #remove empty element
  #found sequence
  if(scalar(@line) == 1 and $line[0] =~ /^[AUGC]+$/){
    $sequence = $line[0];
  }
  #found bar entry
  elsif($line[0] =~ /\d+/ and $line[1] =~ /^[\(\.\)]+$/){
    if($line[1] =~ /^$opt_s1$/ and $line[1] =~ /^$opt_s2$/){
      print STDERR "Structure 1 and 2 should be exclusive\n",$line[1],"\n",$opt_s1,"\n",$opt_s2,"\n";
      
    }
    elsif($line[1] =~ /^$opt_s1$/){
      $states{$line[0]}=[$line[1],"s1"];
    }
    elsif($line[1] =~ /^$opt_s2$/){
      $states{$line[0]}=[$line[1],"s2"];
    }
    else{
      $states{$line[0]}=[$line[1],"s3"];
    }
  }
  else{
    print STDERR "Don't know what to do with this line. Ignored line.\n";
  }
}
close(FI);

foreach my $k (sort {$a <=> $b} keys %states){
  print "# ",$k,"\t",join("\t", @{$states{$k}}),"\n";
}

open(FI, $opt_treekin) or die "can't open ",$opt_treekin,"\n";
while(<FI>){
  if($_ =~ /^#/){ #print comments
    print $_; next; 
  }
  my @line = split(/\s+/, $_);
  shift(@line) if($line[0] eq ""); #remove empty element
  my @curves = ();
  for(my $i=1; $i<scalar(@line); $i++){
    if(exists $states{$i}){
      if($states{$i}->[1] eq "s1"){
	$curves[0]+=$line[$i];
      }
      elsif($states{$i}->[1] eq "s2"){
	$curves[1]+=$line[$i];
      }
      elsif($states{$i}->[1] eq "s3"){
	$curves[2]+=$line[$i];
      }
      else{
	print STDERR "Should never happen1\n";
      }
    }
    else{
      $curves[1]+=$line[$i];
    }
  }
  print $line[0]," ",join(" ", @curves),"\n";
}
close(FI);

__END__
