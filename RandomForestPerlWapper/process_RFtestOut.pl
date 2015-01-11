# This is a program to transform the RF_output format into 
# Each line: pnScore 	label 
# 	

use strict; 
die "Usage: command input_RFoutput out_ScoreLabel_file\n" if scalar(@ARGV) < 2;
my ($input_RFoutput, $out_ScoreLabel_file) = @ARGV;

open(IN, $input_RFoutput) || die(" Can not open file(\"$input_RFoutput\").\n"); 
open(OUT, "> $out_ScoreLabel_file") || die(" Can not open file(\"$out_ScoreLabel_file\").\n");
my $count = 0; 
while (<IN>)	
{
	chomp $_;	
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines
		
	my @cur_line = split('\s+', $_); 
	
	# so we use pos_class_prop - rand_class_prop as score 	
	my $score = $cur_line[ $#cur_line - 1 ] - $cur_line[ $#cur_line ];
	my $label = $cur_line[1];
	print OUT "$score	$label\n"; 
	$count = $count + 1; 
}
print "\nProcess RFoutput File: $input_RFoutput -> $out_ScoreLabel_file\n"; 
print "Total $count lines processed ! \n"; 
close(IN);
close(OUT); 
