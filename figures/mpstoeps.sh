#!/usr/bin/env perl

###########################################################################
# This script converts MetaPost output files (with .mps file extensions ) #
# into stand-alone (self-contained) Encapsulated PostScript files (with   #
# all necessary fonts included.  Although MetaPost figures' final         #
# destination are often inclusion into TeX documents, it is often useful  #
# to view the figures during their development process.  To this end,     #
# stand-alone EPS figures are often desired.                              #
###########################################################################

##################################################
# Author          : Troy Henderson               #
# Created On      : Tue Jun 29 14:27:08 CDT 2004 #
# Last Modified By: Troy Henderson               #
# Last Modified On: Thu Dec 22 04:30:59 EST 2005 #
# Update Count    : 3                            #
##################################################

use POSIX qw(ceil);
use File::Temp;
use File::Copy;
use Cwd;

########################################
# Make sure input argument are correct #
########################################
$0 =~/.*[\/\\](.*)/;$scriptname = $1;

if (@ARGV != 1 & @ARGV != 2) {
	&usage;
}
if (@ARGV == 1 ) {
	$MPSfile=$ARGV[0];
} else {
	&usage;
}

$MPSfile =~ /(.*)\.mps$/; $basename = $1;
if ( $MPSfile ne "$basename\.mps" ) {
	&usage;
}

########################################
# Start the process of embedding fonts #
########################################
$cwd = getcwd;
$dir = File::Temp::tempnam(".","");
mkdir $dir;
copy("$MPSfile","$dir");
chdir($dir);

open(MYFILE,"< $MPSfile");
foreach $line (<MYFILE>) {
   if($line =~ "\%\%BoundingBox: ") {
  	   $line =~ s/\%\%BoundingBox: //;
     	$line =~ s/\n//;
      @line = split(/ /,$line);
  	   $width=$line[2]-$line[0];
     	$height=$line[3]-$line[1];
      $width=$line[2]-$line[0];
  	   $height=$line[3]-$line[1];
     	last;
   }
}
close(MYFILE);

open MYFILE, "> $basename.tex";
#$width+=24;
#$height+=24;
$width++;
$height++;
print MYFILE "\\documentclass[12pt]{article}\n";
print MYFILE "\\usepackage{graphicx}\n";
print MYFILE "\\usepackage[papersize={$width"."bp,$height"."bp},margin={0pt,0pt},noheadfoot,nomarginpar]{geometry}\n";
print MYFILE "\\begin{document}\n";
print MYFILE "\\thispagestyle{empty}\n";
print MYFILE "\\noindent\\includegraphics{$MPSfile}\n";
print MYFILE "\\end{document}";
close(MYFILE);
$command="latex --interaction=nonstopmode $basename.tex > null;
dvips -T $width"."bp,$height"."bp -E -q -o $basename.eps $basename.dvi > null";
close MYFILE;
system($command);

######################################################
# Strip out unwanted lines for the newly created EPS #
######################################################
open MYFILE, "$basename.eps";
open MYFILE2,"> $cwd/$basename.eps";
$BBOX = 0;
foreach $line (<MYFILE>) {
	if ($BBOX == 0  &  $line =~ /^%%BoundingBox:/ ) {
		@nums = $line =~ /(\d+) (\d+) (\d+) (\d+)/;
#		$nums[0]++; $nums[1]++; $nums[2]--; $nums[3]--;
		$nums[0]--; $nums[1]++; $nums[2]--; $nums[3]--;
		$BBOX = 1;
		print MYFILE2 ("%%BoundingBox: $nums[0] $nums[1] $nums[2] $nums[3]\n");
	} else {
		print MYFILE2 $line unless $line =~ /^%DVIPS|^%%Creator|^%%Title/;
	}
}
close MYFILE2;
close MYFILE;

#################
# Clean up mess #
#################
system("rm -f *");
chdir $cwd;

##############################
# Remove temporary directory #
##############################
rmdir($dir);

sub usage {
	print("Usage: $scriptname <MetaPost Output file>\n\n");
	print("Note: MetaPost Output filename MUST have a \.mps extension\n");
	exit;
}
