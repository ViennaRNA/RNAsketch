#!/usr/bin/perl -w
# -*-Perl-*-
# Last changed Time-stamp: <2016-06-20 10:32:57 sven>
# $Id$


use Getopt::Long;
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
use vars qw/$opt_debug $opt_bar $opt_treekin/;
#version 2 local (only in the programm known) variable(s)
my ($opt_s1, $opt_s2 ) = (".........................\\(\\(\\(\\(\\(\\(\\(\\(\\(\\(\\(......\\)\\)\\)\\)\\)\\)\\)\\)\\)\\)\\)...........", "\\(\\(\\(\\(\\(\\.\\.\\.\\(\\(\\(\\(\\(\\(\\(\\(.....\\)\\)\\)\\)\\)\\.\\.\\.\\)\\)\\)\\.\\.\\.\\)\\)\\)\\)\\)........................");
#a=s a must be String
#b=i b must be Integer
#d   d = 1 =>  if 'd' is defined, esle d = undef 
usage() unless GetOptions("d|debug" => \$opt_debug,
                          "s1|structure1=s" => \$opt_s1,
			  "s2|structure2=s" => \$opt_s2,
			  "bi|barmap=s" => \$opt_bar,
			  "ti|treekin=s" => \$opt_treekin
);

#debug message will be printet on standard-error if opt_debug is defined
if($opt_debug){
  print STDERR $usage;
}

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
