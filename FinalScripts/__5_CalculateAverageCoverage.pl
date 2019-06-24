#######################################################################################
###						    Calculate Average Coverage Exome			    		###
#######################################################################################

use strict;
use warnings;
use Getopt::Long;

#######################################################################################
###				 				  Initialize and Declare 			    		    ###
#######################################################################################

my $AnalysisDirMasking 				= "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/exomeBWsamples/masking/analysis/";
my $AnalysisDirNoMasking 			= "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/exomeBWsamples/nomasking/analysis/";
my $ABCC6_LengthTarget 				= 1725;
my $NonPseudoExomeLengthTarget		= 46998255;
my $NonPseudoRegionsFilePath 		= "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/bed/NonPseudo.bed";
my $ABCC6Target						= "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/bed/ABCC6.probes.bed";
my $GenerateBed						= "no";
my %AverageCoverage					= ();

#######################################################################################
###				 				  	 Get Options 			    		     		###
#######################################################################################

GetOptions ("GenerateBed=s"			=>	\$GenerateBed);

#######################################################################################
###				 				  	    Analysis			    		     		###
#######################################################################################

if ($GenerateBed eq "yes"){

	#######################################################################################
	### 	  	 				  			LOOP NOMASKING								###
	#######################################################################################
		
	opendir (NOMASK, $AnalysisDirNoMasking) or die ("Can't open $AnalysisDirNoMasking!\n");
	while (my $BamFile = readdir(NOMASK)){

		if ($BamFile =~ m/^(.*)\.sorted\.nonpseudo\.clipped\.unique\.bam$/){
				
			my 	$Sample 					= $1;
			my 	$NONPSEUDO_BAM_FP 			= $AnalysisDirNoMasking . $BamFile;
			my 	$ABCC6_BAM_FP				= $AnalysisDirNoMasking . $Sample . ".sorted.abcc6.clipped.unique.bam";
			
			my 	$ShellScript				= $AnalysisDirNoMasking . "/$Sample\.cov.sh";
			
			
			open SH, ">$ShellScript" or die ("Can't open $ShellScript\n");
		
			print SH "#PBS -l nodes=1:ppn=1\n";
			print SH "#PBS -l walltime=24:00:00\n";
			print SH "module load BamUtil/1.0.14-intel-2018a\n";
			print SH "module load BEDTools/2.27.1-intel-2018a\n";
			print SH "cd $AnalysisDirNoMasking\n";
			
			print SH "bedtools genomecov -ibam $NONPSEUDO_BAM_FP -bg > $NONPSEUDO_BAM_FP\.cov.bed\n";
		
			print SH "bedtools intersect -a $NONPSEUDO_BAM_FP\.cov.bed -b $NonPseudoRegionsFilePath > $Sample\.sorted.nonpseudo.clipped.unique.cov.bed\n";
			
			print SH "rm $NONPSEUDO_BAM_FP\.cov.bed\n";
					
			print SH "bedtools genomecov -ibam $ABCC6_BAM_FP -bg > $ABCC6_BAM_FP\.cov.bed\n";
			
			print SH "bedtools intersect -a $ABCC6_BAM_FP\.cov.bed -b $ABCC6Target > $Sample\.sorted.abcc6.clipped.unique.cov.bed\n";
			
			print SH "rm $ABCC6_BAM_FP\.cov.bed\n";
			
			close SH;
			
			system("qsub $ShellScript");
		}
	}
	closedir(NOMASK);

	#######################################################################################
	### 	  	 				  	 	 	 LOOP MASKING 								###
	#######################################################################################

	opendir (MASK, $AnalysisDirMasking) or die ("Can't open $AnalysisDirMasking!\n");
	while (my $BamFile = readdir(MASK)){
		
		if ($BamFile =~ m/^(.*)\.sorted\.nonpseudo\.clipped\.unique\.bam$/){
			
			my 	$Sample 					= $1;
			my 	$NONPSEUDO_BAM_FP 			= $AnalysisDirMasking . $BamFile;
			my 	$ABCC6_BAM_FP				= $AnalysisDirMasking . $Sample . ".sorted.abcc6.clipped.unique.bam";
			
			my 	$ShellScript				= $AnalysisDirMasking . "/$Sample\.cov.sh";
			
			open SH, ">$ShellScript" or die ("Can't open $ShellScript\n");
		
			print SH "#PBS -l nodes=1:ppn=1\n";
			print SH "#PBS -l walltime=24:00:00\n";
			print SH "module load BamUtil/1.0.14-intel-2018a\n";
			print SH "module load BEDTools/2.27.1-intel-2018a\n";
			print SH "cd $AnalysisDirMasking\n";
			
			print SH "bedtools genomecov -ibam $NONPSEUDO_BAM_FP -bg > $NONPSEUDO_BAM_FP\.cov.bed\n";
		
			print SH "bedtools intersect -a $NONPSEUDO_BAM_FP\.cov.bed -b $NonPseudoRegionsFilePath > $Sample\.sorted.nonpseudo.clipped.unique.cov.bed\n";
			
			print SH "rm $NONPSEUDO_BAM_FP\.cov.bed\n";
			
			print SH "bedtools genomecov -ibam $ABCC6_BAM_FP -bg > $ABCC6_BAM_FP\.cov.bed\n";
			
			print SH "bedtools intersect -a $ABCC6_BAM_FP\.cov.bed -b $ABCC6Target > $Sample\.sorted.abcc6.clipped.unique.cov.bed\n";
			
			print SH "rm $ABCC6_BAM_FP\.cov.bed\n";
			
			close SH;
			
			system("qsub $ShellScript");
		}
	}
	closedir(MASK);
}
else {

	#######################################################################################
	### 	  	 				  			LOOP NOMASKING								###
	#######################################################################################
		
	opendir (NOMASK, $AnalysisDirNoMasking) or die ("Can't open $AnalysisDirNoMasking!\n");
	while (my $BedFile = readdir(NOMASK)){

		if ($BedFile =~ m/^(.*)\.sorted\.nonpseudo\.clipped\.unique\.cov\.bed$/){
		
			print $BedFile . "\n";
			
			my 	$BedFilePath 										= $AnalysisDirNoMasking . $BedFile;
			my 	$Sample 											= $1;
			
			my 	$AverageCoveragePerBase 							= CalculateAverageCoverage ($BedFilePath, $NonPseudoExomeLengthTarget);
				$AverageCoveragePerBase 							=~ s/\./,/g;
					
			$AverageCoverage{$Sample}{"nonpseudo"}{"nomasking"} 	=  $AverageCoveragePerBase;
		}
		elsif ($BedFile =~ m/^(.*)\.sorted\.abcc6\.clipped\.unique\.cov\.bed$/){
		
			my 	$BedFilePath 										= $AnalysisDirNoMasking . $BedFile;
			my 	$Sample 											= $1;
			
			my 	$AverageCoveragePerBase 							= CalculateAverageCoverage ($BedFilePath, $ABCC6_LengthTarget);
				$AverageCoveragePerBase 							=~ s/\./,/g;
					
			$AverageCoverage{$Sample}{"abcc6"}{"nomasking"} 		=  $AverageCoveragePerBase;
		}
	}
	closedir(NOMASK);

	#######################################################################################
	### 	  	 				  	 	 	 LOOP MASKING 								###
	#######################################################################################

	opendir (MASK, $AnalysisDirMasking) or die ("Can't open $AnalysisDirMasking!\n");
	while (my $BedFile = readdir(MASK)){
		
		if ($BedFile =~ m/^(.*)\.sorted\.nonpseudo\.clipped\.unique\.cov\.bed$/){
		
			print $BedFile . "\n";
			
			my 	$BedFilePath 										= $AnalysisDirMasking . $BedFile;
			my 	$Sample 											= $1;
			
			my 	$AverageCoveragePerBase 							= CalculateAverageCoverage ($BedFilePath, $NonPseudoExomeLengthTarget);
				$AverageCoveragePerBase 							=~ s/\./,/g;
					
			$AverageCoverage{$Sample}{"nonpseudo"}{"masking"} 		=  $AverageCoveragePerBase;
		}
		elsif ($BedFile =~ m/^(.*)\.sorted\.abcc6\.clipped\.unique\.cov\.bed$/){
		
			my 	$BedFilePath 										= $AnalysisDirMasking . $BedFile;
			my 	$Sample 											= $1;
			
			my 	$AverageCoveragePerBase 							= CalculateAverageCoverage ($BedFilePath, $ABCC6_LengthTarget);
				$AverageCoveragePerBase 							=~ s/\./,/g;
					
			$AverageCoverage{$Sample}{"abcc6"}{"masking"} 			=  $AverageCoveragePerBase;
		}
	}
	closedir(MASK);

}

#######################################################################################
### 	  	 				  	    	 Print   	 								###
#######################################################################################

my $OutputFilePath = "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/outputscripts/AverageCoverage.txt";

open RAT, ">$OutputFilePath" or die ("Can't open $OutputFilePath\n");
print RAT "Sample\tExomeNonPseudo_NonMasked\tExomeNonPseudo_Masked\tABCC6_NonMasked\tABCC6_Masked\n";

for my $DNA (sort keys %AverageCoverage){

	print RAT "$DNA\t$AverageCoverage{$DNA}{nonpseudo}{nomasking}\t$AverageCoverage{$DNA}{nonpseudo}{masking}\t$AverageCoverage{$DNA}{abcc6}{nomasking}\t$AverageCoverage{$DNA}{abcc6}{masking}\n";
}

close RAT;

#######################################################################################
### 	  	 				  	    	   SUB   	 								###
#######################################################################################

sub CalculateAverageCoverage {

	(my $CoverageBedFilePath,
	 my $NrOfPositions)			= @_;
	 my $TotalCoverage			= 0;
	
	open BED, "$CoverageBedFilePath" or die ("Can't open $CoverageBedFilePath\n");
	while (<BED>){
		
		my 	$Line = $_;
			$Line =~ s/\n//g;
			
		(undef, my $Start, my $End, my $CovPerBase) = split ("\t", $Line);
		
		$TotalCoverage += $CovPerBase*($End-$Start);
	}
	close BED;
	
	my $AverageCoverage = $TotalCoverage/$NrOfPositions;
	
	return $AverageCoverage;
}