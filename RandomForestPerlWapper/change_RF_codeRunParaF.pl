# If you use this code in your paper , please cite : 
# 
# Y. Qi, HK. Dhiman, et al, Z. Bar-Joseph, J. Klein-Seetharaman,(2009) 
# "Systematic prediction of human membrane receptor interactions"
# PROTEOMICS 2009, 9, 5243-5255
#
# This program is a wrapper to change the RF code's training size and compile, then run the training process 
#
# Here we suppose that there is one program MinGW and we could use the command "g77" to compile F77 code 
# 
# We assume that the input feature file containing all real features and the features having no missing values


use strict;
die "Usage: command inputParaFile  \n" if scalar(@ARGV) < 1;
my ( $inputParaFile ) = @ARGV;


#--------- read in the input para file ----------------

# ( $oldRFcodeFile, $inputTrainFile, $inputTestFile, $outRFcodeFile, $outRFexeFile, $outRFtree, $outRFimpFast ) 
my %inputPara = (); 
my $count = 0; 
open(INP, $inputParaFile) || die(" Can not open file(\"$inputParaFile\").\n"); 
while (<INP>)	
{
	chomp; 
	chop; 
	my $per_line = $_; 
	my @items = split(' = ', $per_line); 	
	$inputPara{"$items[0]"} = $items[1]; 
	print "$per_line\n"; 	
	$count = $count + 1; 	
}
close(INP); 
print "\n==>$count input parameters specified !!! \n";


die "Not enough parameters specified ! \n"  if ($count < 10 ); 


my $oldRFcodeFile = $inputPara{"oldRFcodeFile"}; 
my $inputTrainFile = $inputPara{"inputTrainFile"}; 
my $inputTestFile  = $inputPara{"inputTestFile"}; 
my $outRFcodeFile  = $inputPara{"outRFcodeFile"}; 
my $outRFexeFile  = $inputPara{"outRFexeFile"};
my $outRFtrees =  $inputPara{"outRFtrees"}; 
my $outRFimpFast = $inputPara{"outRFimpFast"}; 
my $mTry = $inputPara{"mTry"}; 
my $numTrees = $inputPara{"numTrees"}; 
my $posCostFactor = $inputPara{"posCostFactor"}; 



#---------------------  count trainSize  -------------------------------

my $numFea = 0; 

open(IN, $inputTrainFile) || die(" Can not open file(\"$inputTrainFile\").\n"); 
my $inputTrainSize = 0; 
while (<IN>)	
{
	chomp;
	my $per_line = $_; 
	my @items = split(' ', $per_line);
	$numFea = $#items;
	$inputTrainSize = $inputTrainSize +1; 		
}
close(IN); 
print "\nInput train data: $inputTrainSize examples. $numFea features. \n"; 



#---------------------  count testSize  -------------------------------

open(IN, $inputTestFile) || die(" Can not open file(\"$inputTestFile\").\n"); 
my $inputTestSize = 0; 
while (<IN>)	
{
	chomp;
	my $per_line = $_; 
	my @items = split(' ', $per_line);
	my $numTestFea = $#items;
	$inputTestSize = $inputTestSize +1; 	
	die "Test feature file line $inputTestSize: numFeature $inputTestSize wrong ! \n"  if ( $numTestFea != $numFea ); 
}
close(IN); 
print "\nInput test data: $inputTestSize examples. \n"; 




#---------------------  change RF code  -------------------------------

open(IN, $oldRFcodeFile) || die(" Can not open file(\"$oldRFcodeFile\").\n"); 
open(OUT, "> $outRFcodeFile") || die(" Can not open file(\"$outRFcodeFile\").\n"); 

my $line_num = 0; 
while (<IN>)	
{
	my $per_line = $_; 
	$line_num = $line_num +1; 	
	
	if ( $per_line =~ m/     1	mdim= / )
	{
		print "Original: ".$per_line;
		my @items = split(/,/, $per_line); 
		
		my $mdimPart = $items[0];
		my @mdimarray = split(/=/, $mdimPart); 
		$mdimarray[1] = $numFea ; 
		$items[0] = join("= ", @mdimarray); 

		my $nsample = $items[1];
		my @samplearray = split(/=/, $nsample); 
		$samplearray[1] = $inputTrainSize ; 
		$items[1] = join("= ", @samplearray); 
		
		my $changeLine = join(",", @items);	
		print "Changed: ".$changeLine."\n";
		print OUT $changeLine; 
	}
	elsif ( $per_line =~ m/     1	maxcat=3,ntest=/ )
	{
		print "Original: ".$per_line;
		my @items = split(/,/, $per_line); 
		
		my $nsample = $items[1];
		my @samplearray = split(/=/, $nsample); 
		$samplearray[1] = $inputTestSize ; 
		$items[1] = join("= ", @samplearray); 
		
		my $changeLine = join(",", @items);	
		print "Changed: ".$changeLine."\n";
		print OUT $changeLine; 
	}
	else {
		if ( $per_line =~ m/     2	jbt=/ )
		{
			print "Original: ".$per_line;
			my @items = split(/,/, $per_line); 
		
			my $mdimPart = $items[0];
			my @mdimarray = split(/=/, $mdimPart); 
			$mdimarray[1] = $numTrees ; 
			$items[0] = join("= ", @mdimarray); 

			my $nsample = $items[1];
			my @samplearray = split(/=/, $nsample); 
			$samplearray[1] = $mTry ; 
			$items[1] = join("= ", @samplearray); 
		
			my $changeLine = join(",", @items);	
			print "Changed: ".$changeLine."\n";
			print OUT $changeLine; 
		}
		elsif ( $per_line =~ m/classwt\(1\)=/ )
		{
			print "Original: ".$per_line;
			my @items = split(/=/, $per_line); 
			$items[1] = "$posCostFactor\n" ; 
			my $changeLine = join("= ", @items);	
			print "Changed: ".$changeLine."\n";
			print OUT $changeLine; 
		}	
		else {
			print OUT $per_line;
		}
	}
}
close(IN); 
close(OUT); 



#---------------------  compile RF train  -------------------------------

my $g77cmd = 'g77 '; 
my $compile = "$g77cmd  $outRFcodeFile -o $outRFexeFile "; 
print "$compile \n"; 
system($compile); 



#---------------------  run RF train  -------------------------------


my $traincmd = "$outRFexeFile $inputTrainFile $outRFtrees  $outRFimpFast $inputTestFile $inputTestFile.out "; 
print "$traincmd \n"; 
system($traincmd); 
