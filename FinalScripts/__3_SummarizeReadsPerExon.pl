# ml load SAMtools/1.8-intel-2018a
#######################################################################################
###								Summarize Counts Per Exon 			    			###
#######################################################################################

use warnings;

#######################################################################################
###				 				  Initialize and Declare 			    			###
#######################################################################################

my @AnalysisDirs 
							= 	(		
									"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/nomasking/analysis/",
									"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/masking/analysis/",
									"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/exomeBWsamples/nomasking/analysis/",
									"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/exomeBWsamples/masking/analysis/"
								); 

my @MaskingAnalysisDirs 
							= 	(	
									"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/masking/analysis/",
									"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/exomeBWsamples/masking/analysis/"
								); 

my @NoMaskingAnalysisDirs 
							= 	(	
									"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/nomasking/analysis/",
									"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/exomeBWsamples/nomasking/analysis/"
								);

my $OutputDir 				= 		"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/outputscripts/";

my %StatFiles 				= (); 
my %CoveragePerExon 		= ();
my %TotalNrOfReads 			= ();
my %TotalReads 				= ();
my %UniqueReads				= ();
my %MultiReads				= ();
my %FailedDNAs				= ();

my @FailedDNAs 				= qw (	Pool-50-Run148_S50 
									Pool-50-Run148_S50 
									Pool-51-Run148_S51 
									Pool-51-Run148_S51 
									Pool-56-Run148_S56 
									Pool-56-Run148_S56 
									Pool-57-Run148_S57 
									Pool-57-Run148_S57 
									Pool-58-Run148_S58 
									Pool-58-Run148_S58 
									Pool-59-Run148_S59 
									Pool-59-Run148_S59);

###################################################################################
###                    				 Failed DNAs                       			###
###################################################################################

foreach my $FailedDNA (@FailedDNAs){

	$FailedDNAs{$FailedDNA} = undef;
}

###################################################################################
###                read in number of ABCC6 reads after masked mapping           ###
###################################################################################

# print "\nread in number of ABCC6 reads after masked mapping\n";

# foreach my $MaskingAnalysisDir (@MaskingAnalysisDirs){

	# $MaskingAnalysisDir =~ m/\/kyukon\/scratch\/gent\/gvo000\/gvo00082\/research\/wst\/abcc6\/(.*)\/masking\/analysis\//;
	
	# my $Enrichment 	= $1;

	# opendir(DIR, $MaskingAnalysisDir) or die $!;
	# while (my $StatFile = readdir(DIR)){

		# if (    $StatFile =~ m/^.*\.sorted\.abcc6\.clipped\.abcc6\.r.*\.stat$/){

				# $StatFile =~ m/^(.*)\.sorted\.abcc6\.clipped\.abcc6\.r(.*)\.stat$/;

			# my $Sample = $1;
			# my $Exon   = $2;
			
			# next if (exists $FailedDNAs{$Sample});
			
			# if ($Exon =~ m/r4/){
				# $Exon =~ s/r4/4/g;
			# }

			# my $StatPath = $MaskingAnalysisDir . $StatFile;

			# open STAT, "$StatPath" or die ("$!");
			# while (<STAT>){

				# next if ($_ !~ m/^chr16/);

				# my $Line = $_;
				   # $Line =~ s/\n//g;

				# (undef, undef, my $NrOfReads, undef) = split ("\t", $Line);
								
				# $TotalNrOfReads{$Enrichment}{$Sample}{$Exon} = $NrOfReads;
			# }
			# close STAT;
		# }
	# }
# }

# print "\nread in number of ABCC6 reads after masked mapping\n";

# foreach my $NoMaskingAnalysisDir (@NoMaskingAnalysisDirs){

	# $NoMaskingAnalysisDir =~ m/\/kyukon\/scratch\/gent\/gvo000\/gvo00082\/research\/wst\/abcc6\/(.*)\/nomasking\/analysis\//;
	
	# my $Enrichment 	= $1;
	
	# opendir(DIR, $NoMaskingAnalysisDir) or die $!;
	# while (my $StatFile = readdir(DIR)){

		# if (    $StatFile =~ m/^.*\.sorted\.abcc6\.clipped\.*\.r.*\.stat$/){

				# $StatFile =~ m/^(.*)\.sorted\.abcc6\.clipped\(.*)\.r(.*)\.stat$/;

			# my $Sample 	= $1;
			# my $Gene	= $2;
			# my $Exon   	= $3;
			
			# next if (exists $FailedDNAs{$Sample});
			
			# if ($Exon =~ m/r4/){
				# $Exon =~ s/r4/4/g;
			# }

			# my $StatPath = $MaskingAnalysisDir . $StatFile;

			# open STAT, "$StatPath" or die ("$!");
			# while (<STAT>){

				# next if ($_ !~ m/^chr16/);

				# my $Line = $_;
				   # $Line =~ s/\n//g;

				# (undef, undef, my $NrOfReads, undef) = split ("\t", $Line);
				
				# if (!exists $TotalNrOfReads{$Enrichment}{$Sample}{$Exon}){
				
					# $TotalNrOfReads{$Enrichment}{$Sample}{$Exon} = $NrOfReads;
				# }
				# else {
				
					# $TotalNrOfReads{$Enrichment}{$Sample}{$Exon} += $NrOfReads;
				# }
			# }
			# close STAT;
		# }
	# }
# }

# for exon 8, in the latest dataset, Ford sequencing by Shana, it appeared to be that the sum of reads that uniquely map to ABCC6 & ABCC6P1 >  nr. of reads that uniquely map on the masked ABCC6 locus. This is caused by the fact that duplicate reads are treated different between a masked and unmasked alignment, i.e. there are more duplicate reads in a masked alignment and so more reads are removed in this. For this reason we check the non uniquely mapped reads in no masking and count how many unique's there are per exon.

print "\nDetermine total number of reads per exon\n\n";

foreach my $NoMaskingAnalysisDir (@NoMaskingAnalysisDirs){

	$NoMaskingAnalysisDir =~ m/\/kyukon\/scratch\/gent\/gvo000\/gvo00082\/research\/wst\/abcc6\/(.*)\/nomasking\/analysis\//;
	
	my $Enrichment 	= $1;

	opendir(DIR, $NoMaskingAnalysisDir) or die $!;
	while (my $BamFile = readdir(DIR)){
		
		if ($BamFile =~ m/^(.*)\.sorted\.abcc6\.clipped\.(.*)\.r(.*)\.bam$/ && $BamFile !~ m/unique/){
		
			my 	$Sample 	= $1;
			my 	$Gene		= $2;
			my 	$Exon   	= $3;
			
			if ($Exon =~ m/r4/){
				$Exon =~ s/r4/4/g;
			}
			
			# print $BamFile . "\n";
			# print $Enrichment . "\t" . $Sample . "\t" . $Gene . "\t" . $Exon . "\n";
		
			my $BamFilePath = $NoMaskingAnalysisDir . $BamFile;
			my $SamFilePath = $BamFilePath;
			   $SamFilePath =~ s/\.bam$/\.sam/g;
			my %ReadNames	= ();
		
			# convert bam to sam #
			
			system ("samtools view -h -o $SamFilePath $BamFilePath");
			
			# open sam and store readnames in hash #
			
			open SAM, "$SamFilePath" or die ("Can't open $SamFilePath\n");
			while (<SAM>){
			
				next if ($_ =~ m/^@/);
				
				my @LineValues 	= split (m/\t/, $_);
				my $ReadName 	= $LineValues[0];
				
				if 		(!exists $TotalReads{$Enrichment}{$Sample}{$Exon}{$ReadName}){
					
					$TotalReads{$Enrichment}{$Sample}{$Exon}{$ReadName} = 1;
				}
				elsif 	($TotalReads{$Enrichment}{$Sample}{$Exon}{$ReadName} == 1){
					
					$TotalReads{$Enrichment}{$Sample}{$Exon}{$ReadName}++;
				}			
			}
			close SAM;
			
			system ("rm $SamFilePath");
		}
	}
}

# Summarize hash

foreach my $Enrichment (keys %TotalReads){
	foreach my $Sample (keys %{$TotalReads{$Enrichment}}){
		foreach my $Exon (keys %{$TotalReads{$Enrichment}{$Sample}}){
			
			my $Total = 0;
		
			foreach my $ReadName (keys %{$TotalReads{$Enrichment}{$Sample}{$Exon}}){
				
				$Total += $TotalReads{$Enrichment}{$Sample}{$Exon}{$ReadName};
			}
				
			$TotalNrOfReads{$Enrichment}{$Sample}{$Exon} = $Total;
		}
	}
}

#######################################################################################
###				       			 Read in all stat files						   		###
#######################################################################################

print "\nRead in files\n";

foreach my $AnalysisDir (@AnalysisDirs){
	
	$AnalysisDir =~ m/\/kyukon\/scratch\/gent\/gvo000\/gvo00082\/research\/wst\/abcc6\/(.*)\/(.*)\/analysis\//;
	
	print "\t" . $AnalysisDir . "\n";
	
	my $Enrichment 	= $1;
	my $Analysis 	= $2;
	my $Type 		= "";
	my $Gene 		= "";
	my $Exon 		= "";

	opendir(DIR, $AnalysisDir) or die $!;
	while (my $StatFile = readdir(DIR)){
	
		if($StatFile =~ m/\.stat$/){
		
			my $StatPath = $AnalysisDir . $StatFile;
			my $DNA 	 = $StatFile;
			   $DNA 	 =~ s/\..*//g;
			   
			next if (exists $FailedDNAs{$DNA});
		
			$StatFiles{$Enrichment}{$DNA} = undef;

			if ($StatPath =~ m/^.*\.sorted\.abcc6\.clipped\.(.*)\.r(.*)\.unique\.stat$/){
				
				$Gene = $1;
				$Exon = $2;
				$Exon =~ s/r//g;
				$Type = "unique";
			}
			elsif ($StatPath =~ m/^.*\.sorted\.abcc6\.clipped\.(.*)\.r(.*)\.stat$/){
				
				 $Gene = $1;
				 $Exon = $2;
				 $Exon =~ s/r//g;
				 $Type = "all";
			}
			else {
				 next;
			}
		
			open STAT, "$StatPath" or die $!;
			while (<STAT>){
				
				next if ($_ !~ m/^chr16/);
				
				my $Line = $_;
				   $Line =~ s/\n//g;
				
				(undef, undef, my $NrOfReads, undef) = split ("\t", $Line);
				
				$CoveragePerExon{$Enrichment}{$Analysis}{$Type}{$Gene}{$Exon}{$DNA} = $NrOfReads;
				
				# Count unique reads for a given exon #
				# Sum for a given exon over ABCC6, ABCC6P1 and ABCC6P2_Exons
				# The difference between all reads and this = multireads
				
				if 		(exists $UniqueReads{$Enrichment}{$Analysis}{$Type}{$Exon}{$DNA} && $Type eq "unique"){
				
					$UniqueReads{$Enrichment}{$Analysis}{$Type}{$Exon}{$DNA} += $NrOfReads;
					
				}
				elsif 	($Type eq "unique") {
				
					$UniqueReads{$Enrichment}{$Analysis}{$Type}{$Exon}{$DNA} = $NrOfReads;
				
				}				
			}
			close STAT;
		}
	}
}

#######################################################################################
###				       				Build Multiread Hash				   		    ###
#######################################################################################

foreach my $Enrichment (keys %UniqueReads){
	foreach my $Analysis (keys %{$UniqueReads{$Enrichment}}){
		foreach my $Type (keys %{$UniqueReads{$Enrichment}{$Analysis}}){
			foreach my $Exon (keys %{$UniqueReads{$Enrichment}{$Analysis}{$Type}}){
				foreach my $DNA (keys %{$UniqueReads{$Enrichment}{$Analysis}{$Type}{$Exon}}){
				
					if (!exists $TotalNrOfReads{$Enrichment}{$DNA}{$Exon}){
						
						$TotalNrOfReads{$Enrichment}{$DNA}{$Exon} = 0;
					}
				
					$MultiReads{$Enrichment}{$Analysis}{$Type}{$Exon}{$DNA} = 	$TotalNrOfReads{$Enrichment}{$DNA}{$Exon} - 
																				$UniqueReads{$Enrichment}{$Analysis}{$Type}{$Exon}{$DNA};
				}
			}
		}
	}
}

%UniqueReads = ();

#######################################################################################
###				       			Extend CoveragePerExon Hash				   		    ###
#######################################################################################

# $CoveragePerExon{$Enrichment}{$Analysis}{$Type}{$Gene}{$Exon}{$DNA} = $NrOfReads;
# $MultiReads{$Enrichment}{$Analysis}{$Type}{$Exon}{$DNA}

foreach my $Enrichment (keys %CoveragePerExon){
	foreach my $Analysis (keys %{$CoveragePerExon{$Enrichment}}){
		foreach my $Type (keys %{$CoveragePerExon{$Enrichment}{$Analysis}}){
			
			foreach my $Exon (keys %{$MultiReads{$Enrichment}{$Analysis}{$Type}}){
				foreach my $DNA (keys %{$MultiReads{$Enrichment}{$Analysis}{$Type}{$Exon}}){
				
					$CoveragePerExon{$Enrichment}{$Analysis}{$Type}{"multi"}{$Exon}{$DNA} = $MultiReads{$Enrichment}{$Analysis}{$Type}{$Exon}{$DNA};
				}
			}
		}
	}
}

#######################################################################################
###				       				  Write Out files					   		    ###
#######################################################################################

my @ABCC6_Exons 	= (1, 2, "3~4", 5, 6, 7, 8, 9);
my @ABCC6P1_Exons 	= (1, 2, "3~4", 5, 6, 7, 8, 9);
my @ABCC6P2_Exons 	= (1, 2, "3~4");

print "\nWrite out files\n";

my $CoverageOutputNormDat  	= $OutputDir . "exoncoverage.norm.dat";

open COVNORMDAT, 	">$CoverageOutputNormDat" 	or die $!;
print COVNORMDAT 	"Gene~Exon\tDNA\tReadsNorm\tEnrichment\tAnalysis\tType\n";


# $CoveragePerExon{$Enrichment}{$Analysis}{$Type}{$Gene}{$Exon}{$DNA} = $NrOfReads;

foreach my $Enrichment (sort keys %CoveragePerExon){
	foreach my $Analysis (sort keys %{$CoveragePerExon{$Enrichment}}){
		foreach my $Type (sort keys %{$CoveragePerExon{$Enrichment}{$Analysis}}){
			
			my $CoverageOutput 			= $OutputDir . $Enrichment . "." . $Analysis . "." . $Type . ".exoncoverage";
			my $CoverageOutputNorm  	= $OutputDir . $Enrichment . "." . $Analysis . "." . $Type . ".exoncoverage.norm";
			
			open COV, 			">$CoverageOutput" 			or die $!;
			open COVNORM, 		">$CoverageOutputNorm" 		or die $!;

			print COV 		"Gene~Exon\t";
			print COVNORM 	"Gene~Exon\t";		
			
			foreach my $DNA (sort keys %{$StatFiles{$Enrichment}}){
				print COV 		"$DNA\t";
				print COVNORM 	"$DNA\t";
			}
			print COV 		"\n";			
			print COVNORM 	"\n";		

			foreach my $Gene (sort keys %{$CoveragePerExon{$Enrichment}{$Analysis}{$Type}}){
			
				my @Exons = ();
				
				if ($Gene eq "abcc6"){
				
					@Exons = @ABCC6_Exons;
				}
				elsif ($Gene eq "abcc6p1"){
				
					@Exons = @ABCC6P1_Exons;
				}
				elsif ($Gene eq "abcc6p2"){
					
					@Exons = @ABCC6P2_Exons;
				}
				elsif ($Gene eq "multi"){
				
					@Exons = @ABCC6_Exons;	
				}
				else {
					die ($Gene);
				}
			
				foreach my $Exon (@Exons){
				 
					# print COV 		"$Gene\~$Exon\t";
					# print COVNORM 	"$Gene\~$Exon\t";
					
					print COV 		"$Exon\t";
					print COVNORM 	"$Exon\t";
					
					if (exists $CoveragePerExon{$Enrichment}{$Analysis}{$Type}{$Gene}{$Exon}){
					
						foreach my $DNA (sort keys %{$CoveragePerExon{$Enrichment}{$Analysis}{$Type}{$Gene}{$Exon}}){
							
							my 	$ReadsNorm			= "NA";
							my 	$PrintEnrichment 	= "";
							
							if ($TotalNrOfReads{$Enrichment}{$DNA}{$Exon} && $TotalNrOfReads{$Enrichment}{$DNA}{$Exon} >= 20){
								
								$ReadsNorm = $CoveragePerExon{$Enrichment}{$Analysis}{$Type}{$Gene}{$Exon}{$DNA}/$TotalNrOfReads{$Enrichment}{$DNA}{$Exon};
							}

							if ($Enrichment =~ m/exome/){
								
								$PrintEnrichment = "WES";
							}
							elsif ($Enrichment =~ m/ford/){
							
								$PrintEnrichment = "GSEA";
							}
							else {
								die ($!);
							}
							
							print COV 		"$CoveragePerExon{$Enrichment}{$Analysis}{$Type}{$Gene}{$Exon}{$DNA}\t";
							print COVNORM 	"$ReadsNorm\t";
							
							print COVNORMDAT "$Gene\~$Exon\t$DNA\t$ReadsNorm\t$PrintEnrichment\t$Analysis\t$Type\n";
						}
						print COV "\n";
						print COVNORM "\n";
					}
					else {
						print "exon not in coverage hash\n";
					}
				}
			}
			close COV;
			close COVNORM;
		}
	}
}

close COVNORMDAT;
