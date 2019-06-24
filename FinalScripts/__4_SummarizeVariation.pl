###################################################################################
### 	 				 Summarize Sequence Differences ABCC6 					###
###################################################################################

use warnings;
use Data::Dumper;

###################################################################################
### 			    		   INITIALIZE & DECLARE 							###
###################################################################################

my @AnalysisDirs
							= 	(		
									"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/masking/analysis/",
									"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/nomasking/analysis/",
									"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/exomeBWsamples/masking/analysis/",
									"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/exomeBWsamples/nomasking/analysis/"
								); 

my @MaskingAnalysisDirs 
							= 	(	
									"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/masking/analysis/",
									"/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/exomeBWsamples/masking/analysis/"
								);
my %SDsInAmp
							= 	( 	"chr16_16317594_16317594_T/G" 	=>	"1\tP2",
									"chr16_16317510_16317510_T/G"	=>	"1\tP1",
									"chr16_16317430_16317430_G/A"	=>	"1\tP2",
									"chr16_16317402_16317402_G/A"	=>	"1\tP1",
									"chr16_16315752_16315752_A/G"	=>	"2\tP1~P2",
									"chr16_16315751_16315751_C/T"	=>	"2\tP2",
									"chr16_16315608_16315608_T/C"	=>	"2\tP2",
									"chr16_16315534_16315534_C/T"	=>	"2\tP2",
									"chr16_16315529_16315528_-/A"	=>	"2\tP1",
									"chr16_16313865_16313865_G/A"	=>	"3\tP2",
									"chr16_16313792_16313792_C/T"	=>	"3\tP2",
									"chr16_16313667_16313667_A/G"	=>	"3\tP1",
									"chr16_16313196_16313196_C/G"	=>	"4\tP1~P2",
									"chr16_16313217_16313217_T/G"	=>	"4\tP1~P2",
									"chr16_16313224_16313224_A/C"	=>	"4\tP2",
									"chr16_16313368_16313368_G/A"	=>	"4\tP1",
									"chr16_16313412_16313412_G/A"	=>	"4\tP2",
									"chr16_16313512_16313512_C/T"	=>	"4\tP1",
									"chr16_16305928_16305928_A/G"	=>	"6\tP1",
									"chr16_16302586_16302586_T/C"	=>	"7\tP1",
									"chr16_16297310_16297310_T/C"	=>	"8\tP1",
									"chr16_16297410_16297410_G/A"	=>	"8\tP1",
									"chr16_16297424_16297424_T/C"	=>	"8\tP1",
									"chr16_16295799_16295799_C/T"	=>	"9\tP1",
									"chr16_16295893_16295893_A/G"	=>	"9\tP1",
									"chr16_16295902_16295902_G/A"	=>	"9\tP1",
									"chr16_16295957_16295957_T/C"	=>	"9\tP1",
									"chr16_16295707_16295707_C/A"	=>	"9\tP1"
								); 

my %SDsAll					= ();
my %Regions					= ();
my %Coverage 				= ();
my %VcfFiles 				= ();
my %SD_Variants 			= ();
my %AllVariants				= ();
my %TotalNrOfReads			= ();
my %SamplesPerVar			= ();
my %FailedDNAs				= ();
my %Histograms				= ();

my $OutputDir 				= "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/outputscripts/";
my $BedFileExons			= "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/bed/ABCC6_Exons.bed";
my $ClustalO_FP				= "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/CLustalO";

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
###                		   determine all sequence differences          			###
###################################################################################

open CLUST, "$ClustalO_FP" or die ("Can't open $ClustalO_FP\n");
while (<CLUST>){
	
	my 	$Line = $_;
		$Line =~ s/\n//g;
		
	(my $ABCC6_Pos, 
	 my $ABCC6_Base, 
	 my $ABCC6P1_Pos, 
	 my $ABCC6P1_Base, 
	 my $ABCC6P2_Pos, 
	 my $ABCC6P2_Base, 
	 my $Comp)
	= split (m/\t/, $Line);
	
	$ABCC6_Base 	=~ s/ //g;
	$ABCC6P1_Base 	=~ s/ //g;
	$ABCC6P2_Base 	=~ s/ //g;
	$ABCC6_Base		=~ tr/ACGTacgt/TGCAtgca/;
	$ABCC6P1_Base	=~ tr/ACGTacgt/TGCAtgca/;
	$ABCC6P2_Base	=~ tr/ACGTacgt/TGCAtgca/;
	
	my $SequenceDiff = "/";
	
	if ($ABCC6_Base 	ne $ABCC6P1_Base 	&& 
		$ABCC6_Base 	=~ m/[A|C|G|T]/ 	&& 
		$ABCC6P1_Base 	=~ m/[A|C|G|T]/ 	&& 
		$ABCC6P2_Base 	eq "-"){ 					# ABCC6P2 is absent.
		
		$SequenceDiff 				= "chr16_$ABCC6_Pos\_$ABCC6_Pos\_$ABCC6_Base\/$ABCC6P1_Base";
		$SDsAll{$SequenceDiff} 		= "0\tP1";
		
	}
	
	elsif (	$ABCC6_Base 	=~ m/[A|C|G|T]/ 	&& 
			$ABCC6P1_Base 	=~ m/[A|C|G|T]/ 	&&
			$ABCC6P2_Base 	=~ m/[A|C|G|T]/){		# ABCC6P2 is present.
		
		if ($ABCC6_Base ne $ABCC6P1_Base){
		
			
			if 		($ABCC6P1_Base eq $ABCC6P2_Base){
				$SequenceDiff 				= "chr16_$ABCC6_Pos\_$ABCC6_Pos\_$ABCC6_Base\/$ABCC6P1_Base";
				$SDsAll{$SequenceDiff} 		= "0\tP1~P2";
			}
			elsif 	($ABCC6_Base eq $ABCC6P2_Base) {
				
				$SequenceDiff 				= "chr16_$ABCC6_Pos\_$ABCC6_Pos\_$ABCC6_Base\/$ABCC6P1_Base";
				$SDsAll{$SequenceDiff} 		= "0\tP1";
			}
			else {
			
				$SequenceDiff 				= "chr16_$ABCC6_Pos\_$ABCC6_Pos\_$ABCC6_Base\/$ABCC6P1_Base";
				$SDsAll{$SequenceDiff} 		= "0\tP1";
				
				$SequenceDiff 				= "chr16_$ABCC6_Pos\_$ABCC6_Pos\_$ABCC6_Base\/$ABCC6P2_Base";
				$SDsAll{$SequenceDiff} 		= "0\tP2";
			}
		}
		elsif ($ABCC6_Base ne $ABCC6P2_Base){
			
			$SequenceDiff 				= "chr16_$ABCC6_Pos\_$ABCC6_Pos\_$ABCC6_Base\/$ABCC6P2_Base";
			$SDsAll{$SequenceDiff} 		= "0\tP2";
		}
	}	
}
close CLUST;

# We add 1 SD from Amp to be compatible, i.e. insertion

$SDsAll{"chr16_16315529_16315528_-/A"} = "0\tP1";

###################################################################################
###                read in number of ABCC6 reads after masked mapping           ###
###################################################################################

# print "\nread in number of ABCC6 reads after masked mapping\n";

# foreach my $MaskingAnalysisDir (@MaskingAnalysisDirs){

	# $MaskingAnalysisDir =~ m/\/kyukon\/scratch\/gent\/gvo000\/gvo00082\/research\/wst\/abcc6\/(.*)\/masking\/analysis\//;
	
	# my $Enrichment 	= $1;

	# opendir(DIR, $MaskingAnalysisDir) or die $!;
	# while (my $StatFile = readdir(DIR)){
	
		# # D1614657.sorted.abcc6.clipped.abcc6p1.r1.stat 
		# # D1614657.sorted.abcc6.clipped.abcc6p1.r1.unique.stat

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
				
				
				# if ($Exon eq "3~4"){
					
					# $TotalNrOfReads{$Enrichment}{$Sample}{3} = $NrOfReads;
					# $TotalNrOfReads{$Enrichment}{$Sample}{4} = $NrOfReads;
				
				# }
				# else {
					# $TotalNrOfReads{$Enrichment}{$Sample}{$Exon} = $NrOfReads;
				# }
			# }
			# close STAT;
		# }
	# }
# }
				
###################################################################################
###                 		  Read In Bed File Exons                      		###
###################################################################################

print "Read in bedfile\n";

open BED, "$BedFileExons" or die ("Can't open $BedFileExons\n");
while (<BED>){
	
	 my	$Line = $_;
		$Line =~ s/\n//g;
		
	(my $Chromosome, my $Start, my $End, my $Region) = split ("\t", $Line);
		
		$Region = lc($Region);
		$Region =~ s/_/\./g;
		$Region =~ s/exon//g;
		
		my 	$Exon = $Region;
			$Exon =~ s/^.*\.//g;
	 
	 $Regions{$Start} = "$Chromosome\t$Start\t$End\t$Exon";
}
close BED;
				
###################################################################################
### 			    			  Read In Variants								###
###################################################################################

print "\nRead in vcf\n"; 

foreach my $AnalysisDir (@AnalysisDirs){
	
	$AnalysisDir =~ m/\/kyukon\/scratch\/gent\/gvo000\/gvo00082\/research\/wst\/abcc6\/(.*)\/(.*)\/analysis\//;
	
	print "\t" . $AnalysisDir . "\n";
	
	my $Enrichment 	= $1;
	my $Analysis	= $2;
	my $Reads		= "";
	my $Caller		= "";
	
	opendir(DIR, $AnalysisDir) or die $!;
	while (my $VcfFile = readdir(DIR)){
		
		if($VcfFile =~ m/\.vcf$/){
				
			if 	($VcfFile =~ m/\.unique\./)	{$Reads = "unique";}
			else	 						{$Reads = "all";}
			
			if 		($VcfFile =~ m/\.gatk35\.default\.vcf$/)		{$Caller = "GATK_Def";}
			elsif 	($VcfFile =~ m/gatk35/)	 						{$Caller = "GATK_mmq0";}
			else 													{$Caller = "LoFreq";}
			
			my 	$DNA = $VcfFile;
				$DNA =~ s/\..*//g;
				
			next if (exists $FailedDNAs{$DNA});
				
			$VcfFiles{$Enrichment}{$Analysis}{$Reads}{$Caller}{$DNA} = $VcfFile;
			
			my $VcfPath = $AnalysisDir . $VcfFile;
			
			open VCF, "$VcfPath" or die $!;
			while(<VCF>){
			
				next if ($_ =~ m/^#/);
				
				my 	$Line = $_;
					$Line =~ s/\n//g;
					
				(my $VariantsRef, my $AFsRef, my $Depth) = DetermineVariantInVCF_Line ($Line, $Caller);
				
				for (my $I = 0; $I < scalar @$VariantsRef; $I++){
				
					my $Variant = $VariantsRef->[$I];
					my $AF		= $AFsRef->[$I];
					
					if ($AF	ne "/"){
					
						if ($SDsInAmp{$Variant} || $SDsAll{$Variant}){
					
							$SD_Variants{$Enrichment}{$Analysis}{$Reads}{$Caller}{$Variant}{$DNA} = "$AF\_$Depth";
						}
						else {
							
							$AllVariants{$Enrichment}{$Analysis}{$Reads}{$Caller}{$Variant}{$DNA} = "$AF\_$Depth";
							$SamplesPerVar{$Variant}{$DNA} = undef;
						}
					}
				}
			}
			close VCF;
		}
	}
	closedir(DIR);
}

###################################################################################
### 			    			  Read In Coverage								###
###################################################################################

print "\nRead in coverage\n"; 

foreach my $AnalysisDir (@AnalysisDirs){
	
	$AnalysisDir =~ m/\/kyukon\/scratch\/gent\/gvo000\/gvo00082\/research\/wst\/abcc6\/(.*)\/(.*)\/analysis\//;
	
	print "\t$AnalysisDir\n";
	
	my $Enrichment  = $1;
	my $Analysis	= $2;
	my $Reads		= "";
	
	opendir(DIR, $AnalysisDir) or die $!;
	while (my $CovFile = readdir(DIR)){
	
		if ($CovFile =~ m/abcc6/ && $CovFile =~ m/\.cov$/){
		
			my $CovPath = $AnalysisDir . $CovFile;
			
			if 	($CovFile =~ m/unique/)	{$Reads = "unique";}
			else 						{$Reads = "all";}
			
			my 	$DNA = $CovFile;
				$DNA =~ s/\..*//g;
				
			next if (exists $FailedDNAs{$DNA});
			
			open COV, "$CovPath" or die $!;
			while (<COV>){
				
			    my 	$Line = $_;
					$Line =~ s/\n//g;
			   
			   (my $Chr, my $Pos, my $Cov) = split("\t", $Line);
					   
				$Coverage{$Enrichment}{$Analysis}{$Reads}{$Pos}{$DNA} = $Cov;
			}
			close COV;
		}
	}
	closedir(DIR);
}

###################################################################################
### 			    			   Write Out SDs								###
###################################################################################

print "\nWrite out SDs in Amplicons\n";

my $SDsOutPathNormDat = $OutputDir . "SDs.dat";

open SDDAT,  ">$SDsOutPathNormDat" 	or die $!;
print SDDAT  "Variation\tDNA\tVAF\tEnrichment\tAnalysis\tType\tCaller\tAmp\n";

foreach my $Enrichment (sort keys %SD_Variants){
	foreach my $Analysis (sort keys %{$SD_Variants{$Enrichment}}){
		foreach my $Reads (sort keys %{$SD_Variants{$Enrichment}{$Analysis}}){
			foreach my $Caller (sort keys %{$SD_Variants{$Enrichment}{$Analysis}{$Reads}}){
			
				# Write separate file per analysis for manual control #
				
				my $OutPathAmp 	= $OutputDir . $Enrichment . "." . $Analysis . "." . $Reads. "." . $Caller . "." . "SDsAmp";
				my $OutPathAll 	= $OutputDir . $Enrichment . "." . $Analysis . "." . $Reads. "." . $Caller . "." . "SDsAll";
				
				### Amp ###
				
				open OUTAMP, 		">$OutPathAmp" or die $!;
				print OUTAMP 		"Variant~Exon~Pi\t";
				
				foreach my $DNA (sort keys %{$VcfFiles{$Enrichment}{$Analysis}{$Reads}{$Caller}}){ # sorted patients
					
					print OUTAMP  "$DNA\t";
				}			
				print OUTAMP "\n";
				
				### All ###
				
				open OUTALL, 		">$OutPathAll" or die $!;
				print OUTALL 		"Variant~Exon~Pi\t";
				
				foreach my $DNA (sort keys %{$VcfFiles{$Enrichment}{$Analysis}{$Reads}{$Caller}}){ # sorted patients
					
					print OUTALL  "$DNA\t";
				}			
				print OUTALL "\n";	
				
				
				# Loop over SDsInAmp #
				
				foreach my $Variant (sort {$b cmp $a} keys %SDsAll){
				
					(my $Pos) 			= $Variant =~ m/^chr\d+_(\d+)_\d+_.*$/;
					(my $Exon, my $Pi) 	= ("/", "/");
					 my $Amp 			= "no";
					 my $Obs			= 0;					
					
					foreach my $DNA (sort keys %{$VcfFiles{$Enrichment}{$Analysis}{$Reads}{$Caller}}){
										
						my $AF 		= "/";
						my $Depth 	= "/";
					
						if 		(exists $SD_Variants{$Enrichment}{$Analysis}{$Reads}{$Caller}{$Variant}{$DNA}){
							
							$Obs++;
						}
						elsif 	(exists $Coverage{$Enrichment}{$Analysis}{$Reads}{$Pos}{$DNA} && 
										$Coverage{$Enrichment}{$Analysis}{$Reads}{$Pos}{$DNA} >= 20){
							$Obs++;
						}
					}
					
					if ($SDsInAmp{$Variant}){
						
						($Exon, $Pi) 	= split ("\t", $SDsInAmp{$Variant});
						 $Amp			= "yes";
						
						print OUTAMP "$Variant\~$Exon\~$Pi\t";
					}
					else {
						($Exon, $Pi) = split ("\t", $SDsAll{$Variant});
					}

					if ($Obs >= 1){
						print OUTALL "$Variant\~$Exon\~$Pi\t";
					}
					
					if ($Obs >= 1){
					
						foreach my $DNA (sort keys %{$VcfFiles{$Enrichment}{$Analysis}{$Reads}{$Caller}}){
											
							my $AF 		= "/";
							my $Depth 	= "/";
						
							if 		(exists $SD_Variants{$Enrichment}{$Analysis}{$Reads}{$Caller}{$Variant}{$DNA}){
								
								($AF, $Depth) = split(m/\_/, $SD_Variants{$Enrichment}{$Analysis}{$Reads}{$Caller}{$Variant}{$DNA});
								
								$AF =~ s/AF\=//g;
							}
							elsif 	(exists $Coverage{$Enrichment}{$Analysis}{$Reads}{$Pos}{$DNA} && $Coverage{$Enrichment}{$Analysis}{$Reads}{$Pos}{$DNA} >= 20){
								
								$AF 	= 0;
								$Depth 	= $Coverage{$Enrichment}{$Analysis}{$Reads}{$Pos}{$DNA};
							}
							else {
								$AF 	= "NA";
							}
							
							if ($SDsInAmp{$Variant}){
								print OUTAMP "$AF\t";
							}
							print OUTALL "$AF\t";
							
							
							my $EnrichmentPrint = "";
							
							if 		($Enrichment =~ m/ford/){
							
								$EnrichmentPrint = "GSEA";
							}
							elsif 	($Enrichment =~ m/exome/){
							
								$EnrichmentPrint = "WES";
							}
							else {die ($!);}
																		
							print SDDAT "$Variant\~$Exon\~$Pi\t$DNA\t$AF\t$EnrichmentPrint\t$Analysis\t$Reads\t$Caller\t$Amp\n";
						}
						
						if ($SDsInAmp{$Variant}){
							print OUTAMP "\n";
						}
						print OUTALL "\n";
					}
				}	
				close OUTAMP;
				close OUTALL;
			}
		}
	}
}
close SDDAT;


###################################################################################
### 			    			 Write Out Other Vars							###
###################################################################################

print "\nWrite out other variants\n";

my $VarsOutPathNormDat = $OutputDir . "OtherVars.dat";

open VARDAT,  ">$VarsOutPathNormDat" 	or die $!;
print VARDAT  "Variation\tDNA\tVAF\tDepth\tEnrichment\tAnalysis\tType\tCaller\n";

foreach my $Enrichment (sort keys %AllVariants){
	
	# Write also separate files per enrichment #
				
	my $OutPathAmp 	= $OutputDir . $Enrichment . ".VarComp";
	
	open VARPERENR, ">$OutPathAmp" or die ("Can't open $OutPathAmp\n");
	print VARPERENR "ID\tVariantDNA\tType\tCaller\tVAF\tDepth\n";
	
	foreach my $Analysis (sort keys %{$AllVariants{$Enrichment}}){
		foreach my $Reads (sort keys %{$AllVariants{$Enrichment}{$Analysis}}){
			foreach my $Caller (sort keys %{$AllVariants{$Enrichment}{$Analysis}{$Reads}}){
				foreach my $Variant (sort keys %{$AllVariants{$Enrichment}{$Analysis}{$Reads}{$Caller}}){
					
					(my $Pos) 	= $Variant =~ m/^chr\d+_(\d+)_\d+_.*$/;
					
					foreach my $DNA (sort keys %{$SamplesPerVar{$Variant}}){
										
						my $AF 		= "/";
						my $Depth 	= "/";
					
						if 		(exists $AllVariants{$Enrichment}{$Analysis}{$Reads}{$Caller}{$Variant}{$DNA}){
							
							($AF, $Depth) = split(m/\_/, $AllVariants{$Enrichment}{$Analysis}{$Reads}{$Caller}{$Variant}{$DNA});
							
							$AF =~ s/AF\=//g;
						}
						elsif 	(exists $Coverage{$Enrichment}{$Analysis}{$Reads}{$Pos}{$DNA} && $Coverage{$Enrichment}{$Analysis}{$Reads}{$Pos}{$DNA} >= 20) {
							
							$AF 	= 0;
							$Depth 	= $Coverage{$Enrichment}{$Analysis}{$Reads}{$Pos}{$DNA};
						}
						else {
							$AF 	= "NA";
						}
						
						print VARPERENR "$Variant\|$DNA\|$Reads\|$Caller\t$Variant\|$DNA\|$Caller\t$Reads\t$Caller\t$AF\t$Depth\n";
												
						
						my $EnrichmentPrint = "";
						
						if 		($Enrichment =~ m/ford/){
						
							$EnrichmentPrint = "GSEA";
						}
						elsif 	($Enrichment =~ m/exome/){
						
							$EnrichmentPrint = "WES";
						}
						else {die ($!);}
																	
						print VARDAT "$Variant\t$DNA\t$AF\t$Depth\t$EnrichmentPrint\t$Analysis\t$Reads\t$Caller\n";
					}
				}
			}
		}
	}
	
	close VARPERENR;
}
close VARDAT;

###################################################################################
### 			    			   Write Out Hists								###
###################################################################################

# print "\nWrite out histograms\n";

# my $HistogramOutPath = $OutputDir . "Hist.txt";
# open HIST, ">$HistogramOutPath" or die ("Can't open $HistogramOutPath\n");

# foreach my $Variant (sort keys %{$SD_Variants{"exomeBWsamples"}{"masking"}{"unique"}{"GATK_Def"}}){

	# print $Variant . "\n";

	# (my $Pos) 	= $Variant =~ m/^chr\d+_(\d+)_\d+_.*$/;

	# print HIST $Variant . "\t";
	
	# foreach my $DNA (sort keys %{$VcfFiles{"exomeBWsamples"}{"masking"}{"unique"}{"GATK_Def"}}){
	
		# if 		(exists $SD_Variants{"exomeBWsamples"}{"masking"}{"unique"}{"GATK_Def"}{$Variant}{$DNA}){
		
			# (my $AF, my $Depth) = split (m/_/, $SD_Variants{"exomeBWsamples"}{"masking"}{"unique"}{"GATK_Def"}{$Variant}{$DNA});
			
			# if ($Depth > 20){
				
				# print HIST $AF . "\t";
			# }
			# else {
				# print HIST "NA\t";
			# }
		# }
		
		# elsif 	(exists $Coverage{"exomeBWsamples"}{"masking"}{"unique"}{$Pos}{$DNA} && 
						# $Coverage{"exomeBWsamples"}{"masking"}{"unique"}{$Pos}{$DNA} >= 20){
			
			# print HIST "0\t";
		# }
		# else {
		
			# print HIST "NA\t";
		# }
	# }
	# print HIST "\n";
# }

# close HIST;

###################################################################################
### 			    				Subroutines 								###
###################################################################################

sub FormatGATKCoordinatesOfVariant{
	
	(my $ReferenceAllele, my $VariantAllele, my $Position) = @_;
	
	do {
	
		if (substr($ReferenceAllele, -1, 1) eq substr($VariantAllele, -1, 1)){
			$ReferenceAllele = substr($ReferenceAllele, 0, (length($ReferenceAllele)-1));
			$VariantAllele = substr($VariantAllele, 0,(length($VariantAllele)-1));
		}
		elsif (substr($ReferenceAllele, 0, 1) eq substr($VariantAllele, 0, 1)){
			$ReferenceAllele = substr($ReferenceAllele, 1, (length($ReferenceAllele)-1));
			$VariantAllele = substr($VariantAllele, 1,(length($VariantAllele)-1));
			$Position++;
		}
		else {
			print $ReferenceAllele . "\t" . $VariantAllele . "\t" . $Position . "\n";
			last;
		}
	}
	until ($ReferenceAllele eq "" || $VariantAllele eq "");
	
	if ($ReferenceAllele eq "") {$ReferenceAllele = "-";}
	else {$VariantAllele = "-";}
	
	return ($ReferenceAllele, $VariantAllele, $Position);
}

sub DetermineVariantInVCF_Line{
				
	(my $Line,
	 my $Caller)
	= @_;
	
	my @LineValues 			= split ("\t", $Line);
	
	my $Chrom				= $LineValues[0];
	my $Pos					= $LineValues[1];
	my $ReferenceBase		= $LineValues[3];
	my $VariantBases		= $LineValues[4];
	my $Info 				= $LineValues[7];
	my $Depth				= "/";
	my $Variant				= "/";
	
	my @AFs					= ();
	my @Variants			= ();
	my @VariantBases		= split (m/\,/, $VariantBases);
	
	foreach my $VariantBase (@VariantBases){
	
		### Determine Variant ###

		if (length($ReferenceBase) == length ($VariantBase)){ 		# SUBSTITUTIONS
			
			$Variant = $Chrom . "_" . $Pos . "_" . $Pos . "_" .  "$ReferenceBase\/$VariantBase";
		}
		elsif (length($ReferenceBase) < length ($VariantBase)){ 	# INSERTION
			
			($ReferenceBase, $VariantBase, $Pos) = FormatGATKCoordinatesOfVariant($ReferenceBase, $VariantBase, $Pos);
			
			my $StartPos 	= $Pos;
			my $EndPos  	= $Pos - 1;
			
			$Variant = $Chrom . "_" . $StartPos . "_" . $EndPos . "_" .  "$ReferenceBase\/$VariantBase";
		}
		else { 														# DELETION
		
			($ReferenceBase, $VariantBase, $Pos) = FormatGATKCoordinatesOfVariant($ReferenceBase, $VariantBase, $Pos);
			
			my $StartPos = $Pos;
			my $EndPos   = $Pos + length($ReferenceBase) - 1;
			
			$Variant = $Chrom . "_" . $StartPos . "_" . $EndPos . "_" .  "$ReferenceBase\/$VariantBase";
		}
		
		push (@Variants, $Variant);
	}
	
	##########################
	
	if 		($Caller =~	m/GATK/){
		
		my 	$Format 	= $LineValues[8];
		my 	$Sample 	= $LineValues[9];
		
		# print $Caller . "\t" . $Format . "\t" . $Sample . "\n";
		
		if ($Format eq "GT:AD:DP:GQ:PL"){
		
			$Sample =~ m/^.*:(.*):(.*):.*:.*$/;
			
			if (!$1){
			
				print $Caller . "\t" . $Line . "\n";
			
			}
			
			# print $1 . "\t" . $2 . "\n";
			
			my 	@ADs 	= split (m/\,/, $1);
				$Depth	= $2;
			
			my $RefAD  	= shift (@ADs);
			
			foreach my $AltAD (@ADs){
			
				if ($AltAD+$RefAD==0){
					push(@AFs, "/");
				}
				else {
					push(@AFs, ($AltAD/($AltAD+$RefAD)));
				}
			}
		}
		else {
		
			die ("Format does not match: GT:AD:DP:GQ:PL\n");
		}							
	}
	elsif 	($Caller eq "LoFreq"){
		
		($Depth, $AFs, undef, undef) = split (/\;/, $Info);
		$Depth =~ s/DP=//g;
		
		@AFs = ($AFs);
	}
	else {
	
		die ("Unknown caller\n");
	}
	
	return \@Variants, \@AFs, $Depth;
}