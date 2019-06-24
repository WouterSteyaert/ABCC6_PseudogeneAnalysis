#################################################################################################################
#################################################################################################################
###						        			 Compare Biological Variation			    				  	  ###
#################################################################################################################
#################################################################################################################

use strict;
use warnings;
use Getopt::Long;

#######################################################################################
###						      	   Initialize & declare		    					###
#######################################################################################

my $ClustalO_FP				= "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/CLustalO";
my $OutputDir 				= "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/outputscripts/";
my $BV_DAT_FilePath			= $OutputDir . "BioVar.dat";
my $VariantDNA_Rank			= 1;
my %Exons					= ();
my %SDsAll					= ();
my %Coverage				= ();
my %Analysis				= ();
my %AllVariants				= ();
my %VariantDNA_Rank			= ();
my %MaxVariantAF			= ();
my %DiffPos					= ();
my @AnalysisDirs
							= 	(		
									"/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/fordBWsamples/masking/analysis/",
									"/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/fordBWsamples/nomasking/analysis/",
									"/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/exomeBWsamples/masking/analysis/",
									"/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/exomeBWsamples/nomasking/analysis/"
								); 

my @MaskingAnalysisDirs 
							= 	(	
									"/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/fordBWsamples/masking/analysis/",
									"/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/exomeBWsamples/masking/analysis/"
								);


#######################################################################################
###                		   	determine all sequence differences          			###
#######################################################################################

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
		$DiffPos{$ABCC6_Pos} = undef;
	}
	
	elsif (	$ABCC6_Base 	=~ m/[A|C|G|T]/ 	&& 
			$ABCC6P1_Base 	=~ m/[A|C|G|T]/ 	&&
			$ABCC6P2_Base 	=~ m/[A|C|G|T]/){		# ABCC6P2 is present.
		
		if ($ABCC6_Base ne $ABCC6P1_Base){
		
			
			if 		($ABCC6P1_Base eq $ABCC6P2_Base){
				$SequenceDiff 				= "chr16_$ABCC6_Pos\_$ABCC6_Pos\_$ABCC6_Base\/$ABCC6P1_Base";
				$SDsAll{$SequenceDiff} 		= "0\tP1~P2";
				
				$DiffPos{$ABCC6_Pos} = undef;
			}
			elsif 	($ABCC6_Base eq $ABCC6P2_Base) {
				
				$SequenceDiff 				= "chr16_$ABCC6_Pos\_$ABCC6_Pos\_$ABCC6_Base\/$ABCC6P1_Base";
				$SDsAll{$SequenceDiff} 		= "0\tP1";
				
				$DiffPos{$ABCC6_Pos} = undef;
			}
			else {
			
				$SequenceDiff 				= "chr16_$ABCC6_Pos\_$ABCC6_Pos\_$ABCC6_Base\/$ABCC6P1_Base";
				$SDsAll{$SequenceDiff} 		= "0\tP1";
				
				$SequenceDiff 				= "chr16_$ABCC6_Pos\_$ABCC6_Pos\_$ABCC6_Base\/$ABCC6P2_Base";
				$SDsAll{$SequenceDiff} 		= "0\tP2";
				
				$DiffPos{$ABCC6_Pos} = undef;
			}
		}
		elsif ($ABCC6_Base ne $ABCC6P2_Base){
			
			$SequenceDiff 				= "chr16_$ABCC6_Pos\_$ABCC6_Pos\_$ABCC6_Base\/$ABCC6P2_Base";
			$SDsAll{$SequenceDiff} 		= "0\tP2";
			
			$DiffPos{$ABCC6_Pos} = undef;
		}
	}	
}
close CLUST;

# We add 1 SD from Amp to be compatible, i.e. insertion

$SDsAll{"chr16_16315529_16315528_-/A"} = "0\tP1";

$DiffPos{16315529} = undef;
$DiffPos{16315528} = undef;

#######################################################################################
### 			    			  	Read In Coverage								###
#######################################################################################

print "\nRead in coverage\n"; 

foreach my $AnalysisDir (@AnalysisDirs){
	
	$AnalysisDir =~ m/\/kyukon\/scratch\/gent\/vo\/000\/gvo00082\/research\/wst\/abcc6\/(.*)\/(.*)\/analysis\//;
	
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
			
			$Analysis{$Enrichment}{$Analysis}{$Reads} = undef;
			
			my 	$DNA = $CovFile;
				$DNA =~ s/\..*//g;
			
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

#######################################################################################
### 			    			    Read In Variants								###
#######################################################################################

print "\nRead in vcf\n"; 

foreach my $AnalysisDir (@AnalysisDirs){
	
	$AnalysisDir =~ m/\/kyukon\/scratch\/gent\/vo\/000\/gvo00082\/research\/wst\/abcc6\/(.*)\/(.*)\/analysis\//;
	
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
					
					if ($AF	ne "/" && $AF >= 0.05){
					
						if ($DNA eq "D1613848" && $Variant eq "chr16_16305928_16305928_A/G"){
							
							print "Found\t$AF\n";
						}
					
						# if (!exists $SDsAll{$Variant}){
					
							$AllVariants{$DNA}{$Variant}{$Enrichment}{$Analysis}{$Reads}{$Caller} = "$AF\_$Depth";
							
							if (!exists $MaxVariantAF{$DNA}{$Variant}{$Caller}){
								
								$MaxVariantAF{$DNA}{$Variant}{$Caller} = $AF;
							}
							elsif ($MaxVariantAF{$DNA}{$Variant}{$Caller} < $AF){
							
								$MaxVariantAF{$DNA}{$Variant}{$Caller} = $AF;
							}
						# }
					}
				}	
			}
			close VCF;
		}
	}
	closedir(DIR);
}

#######################################################################################
### 			    			 Read In Exon Positions								###
#######################################################################################

my $RegionsFilePath = "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/bed/ABCC6_ABCC6P1_ABCC6P2_Exons.bed";

open REG, "$RegionsFilePath" or die ("Can't open $RegionsFilePath\n");
while (<REG>){
	
	my 	$Line = $_;
		$Line =~ s/\n//g;
	
	next if ($Line !~ m/ABCC6_/);
		
	(undef, my $Start, my $End, my $Exon) = split (m/\t/, $Line);
	
	$Exon =~ s/Region_ABCC6_//g;
	$Exon =~ s/ABCC6_Exon//g;
	$Exon =~ s/Ex//g;
	
	for (my $I = $Start; $I <= $End; $I++){
		
		$Exons{$I} = $Exon;
	}
}
close REG;

#######################################################################################
### 			    		   Biological Variation DAT file						###
#######################################################################################

# only take valid variants: those ones where there is sufficient coverage. 
# if variant not present: insert a 0

my $Caller = "GATK_Def";

open BV, ">$BV_DAT_FilePath" or die ("Can't open $BV_DAT_FilePath\n");
print BV "FullAnalysis\tEnrichment\tAnalysis\tReads\tCaller\tVariant\tExon\tDNA\tVariant-DNA\tCol\tAF\tDepth\n";

foreach my $DNA (keys %AllVariants){
	foreach my $Variant (keys %{$AllVariants{$DNA}}){
	
		next if (!exists $MaxVariantAF{$DNA}{$Variant}{$Caller} || 
				$MaxVariantAF{$DNA}{$Variant}{$Caller} < 0.25);
	
		(undef, my $Pos, undef, my $Change) = split (m/_/, $Variant);
		
		# next if (exists $DiffPos{$Pos});
	
		# next if (!exists $Exons{$Pos});
		
		# Check if variant has sufficient coverage in all analysis #
							
		my $SuffCov = 1;
		
		foreach my $Enr (keys %Coverage){
			foreach my $An (keys %{$Coverage{$Enr}}){
				foreach my $R (keys %{$Coverage{$Enr}{$An}}){
					
					if (!exists $Coverage{$Enr}{$An}{$R}{$Pos}{$DNA} || 
						$Coverage{$Enr}{$An}{$R}{$Pos}{$DNA} < 20){
						$SuffCov = 0;
					}
				}
			}
		}
		
		##############################################################
		
		if ($SuffCov){
		
			# next if (!exists $Exons{$Pos});
		
			foreach my $Enrichment (keys %Analysis){
				
				my $EnrichmentP = "";
				
				if 		($Enrichment =~ m/^f/){$EnrichmentP = "f";}
				elsif 	($Enrichment =~ m/^e/){$EnrichmentP = "e";}
			
				foreach my $Analysis (keys %{$Analysis{$Enrichment}}){
				
					my $AnalysisP = "";
				
					if 		($Analysis =~ m/^n/){$AnalysisP = "noma";}
					elsif 	($Analysis =~ m/^m/){$AnalysisP = "ma";}
				
					foreach my $Reads (keys %{$Analysis{$Enrichment}{$Analysis}}){
					
						my $ReadsP = "";
				
						if 		($Reads =~ m/^a/){$ReadsP = "a";}
						elsif 	($Reads =~ m/^u/){$ReadsP = "u";}
						
						if (!exists $VariantDNA_Rank{"$Variant\-$DNA"}){
							$VariantDNA_Rank{"$Variant\-$DNA"} = $VariantDNA_Rank;
							$VariantDNA_Rank++;
						}
						
						if (! exists $Exons{$Pos}){$Exons{$Pos} = "/";}
						
						if (exists $AllVariants{$DNA}{$Variant}{$Enrichment} 						&&
							exists $AllVariants{$DNA}{$Variant}{$Enrichment}{$Analysis} 			&&
							exists $AllVariants{$DNA}{$Variant}{$Enrichment}{$Analysis}{$Reads}{$Caller}){
							
							(my $AF, my $Depth) = split (m/_/, $AllVariants{$DNA}{$Variant}{$Enrichment}{$Analysis}{$Reads}{$Caller});							
								
							if (!exists $SDsAll{$Variant}){
								
								print BV "$EnrichmentP\-$AnalysisP\-$ReadsP\t$Enrichment\t$Analysis\t$Reads\t$Caller\t$Variant\t$Exons{$Pos}\t$DNA\t$Variant\-$DNA\t$VariantDNA_Rank{\"$Variant\-$DNA\"}\t$AF\t$Depth\tNoSD\n";
							}
							else {
								
								print BV "$EnrichmentP\-$AnalysisP\-$ReadsP\t$Enrichment\t$Analysis\t$Reads\t$Caller\t$Variant\t$Exons{$Pos}\t$DNA\t$Variant\-$DNA\t$VariantDNA_Rank{\"$Variant\-$DNA\"}\t$AF\t$Depth\tSD\n";
							}
						}
						else {
						
							if (!exists $SDsAll{$Variant}){
							
								print BV "$EnrichmentP\-$AnalysisP\-$ReadsP\t$Enrichment\t$Analysis\t$Reads\t$Caller\t$Variant\t$Exons{$Pos}\t$DNA\t$Variant\-$DNA\t$VariantDNA_Rank{\"$Variant\-$DNA\"}\t0\t$Coverage{$Enrichment}{$Analysis}{$Reads}{$Pos}{$DNA}\tNoSD\n";
								
							}
							else {
							
								print BV "$EnrichmentP\-$AnalysisP\-$ReadsP\t$Enrichment\t$Analysis\t$Reads\t$Caller\t$Variant\t$Exons{$Pos}\t$DNA\t$Variant\-$DNA\t$VariantDNA_Rank{\"$Variant\-$DNA\"}\t0\t$Coverage{$Enrichment}{$Analysis}{$Reads}{$Pos}{$DNA}\tSD\n";
							
							}
						}
					}
				}
			}
		}
	}
}

close BV;

#######################################################################################
### 			    		   		Analyze Variation								###
#######################################################################################

foreach my $DNA (keys %AllVariants){
	foreach my $Variant (keys %{$AllVariants{$DNA}}){
	
		next if ("$Variant\-$DNA" eq "chr16_16313545_16313545_C/T-D1500757");
		next if ("$Variant\-$DNA" eq "chr16_16313398_16313398_C/T-D1613541");
		next if ("$Variant\-$DNA" eq "chr16_16306145_16306145_G/T-D1500754");
	
		next if (!exists $MaxVariantAF{$DNA}{$Variant}{$Caller} || 
				$MaxVariantAF{$DNA}{$Variant}{$Caller} < 0.25);
	
		(undef, my $Pos, undef, my $Change) = split (m/_/, $Variant);
		
		if (! exists $Exons{$Pos}){$Exons{$Pos} = "/";}
									
		my $SuffCov = 1;
		
		foreach my $Enr (keys %Coverage){
			foreach my $An (keys %{$Coverage{$Enr}}){
				foreach my $R (keys %{$Coverage{$Enr}{$An}}){
					
					if (!exists $Coverage{$Enr}{$An}{$R}{$Pos}{$DNA} || 
						$Coverage{$Enr}{$An}{$R}{$Pos}{$DNA} < 20){
						$SuffCov = 0;
					}
				}
			}
		}
		
		##############################################################
		
		if ($SuffCov){
		
			if ( exists $AllVariants{$DNA}{$Variant}{"exomeBWsamples"}{"nomasking"}{"unique"}{"GATK_Def"} 	|| 
				 exists $AllVariants{$DNA}{$Variant}{"fordBWsamples"}{"masking"}{"unique"}{"GATK_Def"}){
				 
				my $AF_WES 		= 0;
				my $AF_GSEA		= 0;
				my $Depth_WES	= 0;
				my $Depth_GSEA	= 0;
				
				my $AF_WES_M	= 0;
				my $Depth_WES_M	= 0;
				
				if (exists $AllVariants{$DNA}{$Variant}{"exomeBWsamples"}{"nomasking"}{"unique"}{"GATK_Def"}){
				
					($AF_WES, 
					 $Depth_WES) = split (m/_/, $AllVariants{$DNA}{$Variant}{"exomeBWsamples"}{"nomasking"}{"unique"}{"GATK_Def"});
				}
				
				if (exists $AllVariants{$DNA}{$Variant}{"fordBWsamples"}{"masking"}{"unique"}{"GATK_Def"}){
				
					($AF_GSEA, 
					 $Depth_GSEA) = split (m/_/, $AllVariants{$DNA}{$Variant}{"fordBWsamples"}{"masking"}{"unique"}{"GATK_Def"});
				}
				
				if (exists $AllVariants{$DNA}{$Variant}{"exomeBWsamples"}{"masking"}{"unique"}{"GATK_Def"}){
				
					($AF_WES_M, 
					 $Depth_WES_M) = split (m/_/, $AllVariants{$DNA}{$Variant}{"exomeBWsamples"}{"masking"}{"unique"}{"GATK_Def"});
				}
				
				if ($AF_WES >= 0.25 || $AF_GSEA >=0.25){
				
					if (!exists $SDsAll{$Variant}){
					
						# no SDs #
						
						print "exomeBWsamples-unmasked-unique\twes\tunmasked\tunique\tGATK\t$Variant\t$Exons{$Pos}\t$DNA\t$Variant\-$DNA\t$VariantDNA_Rank{\"$Variant\-$DNA\"}\t$AF_WES\t$Coverage{exomeBWsamples}{nomasking}{unique}{$Pos}{$DNA}\t$Depth_WES\tNoSD\n";
						
						print "fordBWsamples-masked-unique\tgsea\tmasked\tunique\tGATK\t$Variant\t$Exons{$Pos}\t$DNA\t$Variant\-$DNA\t$VariantDNA_Rank{\"$Variant\-$DNA\"}\t$AF_GSEA\t$Coverage{fordBWsamples}{masking}{unique}{$Pos}{$DNA}\t$Depth_GSEA\tNoSD\n";
					}
					else {
					
						# SDs #
						
						if ($Variant eq "chr16_16302586_16302586_T/C"){
						
							print "exomeBWsamples-unmasked-unique\twes\tunmasked\tunique\tGATK\t$Variant\t$Exons{$Pos}\t$DNA\t$Variant\-$DNA\t$VariantDNA_Rank{\"$Variant\-$DNA\"}\t$AF_WES\t$Coverage{exomeBWsamples}{nomasking}{unique}{$Pos}{$DNA}\t$Depth_WES\tSD-TP\n";
						
							print "fordBWsamples-masked-unique\tgsea\tmasked\tunique\tGATK\t$Variant\t$Exons{$Pos}\t$DNA\t$Variant\-$DNA\t$VariantDNA_Rank{\"$Variant\-$DNA\"}\t$AF_GSEA\t$Coverage{fordBWsamples}{masking}{unique}{$Pos}{$DNA}\t$Depth_GSEA\tSD-TP\n";
							
						}
						elsif (exists $AllVariants{$DNA}{$Variant}{"exomeBWsamples"}{"masking"}{"unique"}{"GATK_Def"} && 
									$AF_WES_M >= 0.70) {
							
							# TP
							
							print "exomeBWsamples-unmasked-unique\twes\tunmasked\tunique\tGATK\t$Variant\t$Exons{$Pos}\t$DNA\t$Variant\-$DNA\t$VariantDNA_Rank{\"$Variant\-$DNA\"}\t$AF_WES\t$Coverage{exomeBWsamples}{nomasking}{unique}{$Pos}{$DNA}\t$Depth_WES\tSD-TP\n";
						
							print "fordBWsamples-masked-unique\tgsea\tmasked\tunique\tGATK\t$Variant\t$Exons{$Pos}\t$DNA\t$Variant\-$DNA\t$VariantDNA_Rank{\"$Variant\-$DNA\"}\t$AF_GSEA\t$Coverage{fordBWsamples}{masking}{unique}{$Pos}{$DNA}\t$Depth_GSEA\tSD-TP\n";
							
						}
						else {
						
							# FP 
						
							print "exomeBWsamples-unmasked-unique\twes\tunmasked\tunique\tGATK\t$Variant\t$Exons{$Pos}\t$DNA\t$Variant\-$DNA\t$VariantDNA_Rank{\"$Variant\-$DNA\"}\t$AF_WES\t$Coverage{exomeBWsamples}{nomasking}{unique}{$Pos}{$DNA}\t$Depth_WES\tSD-FP\n";
						
							print "fordBWsamples-masked-unique\tgsea\tmasked\tunique\tGATK\t$Variant\t$Exons{$Pos}\t$DNA\t$Variant\-$DNA\t$VariantDNA_Rank{\"$Variant\-$DNA\"}\t$AF_GSEA\t$Coverage{fordBWsamples}{masking}{unique}{$Pos}{$DNA}\t$Depth_GSEA\tSD-FP\n";

						}
					}
				}
			}
		}
	}
}


#######################################################################################
### 			    			       Subroutines								 	###
#######################################################################################

sub DetermineVariantInVCF_Line{
				
	(my $Line,
	 my $Caller)
	= @_;
	
	my @LineValues 			= split ("\t", $Line);
	
	my $Chrom				= $LineValues[0];
	my $Pos					= $LineValues[1];
	my $ReferenceBase		= $LineValues[3];
	my $VariantBases		= $LineValues[4];
	my $Filter				= $LineValues[6];
	my $Info 				= $LineValues[7];
	my $Depth				= "/";
	my $Variant				= "/";
	my $QD					= 3;
	
	my @AFs					= ();
	my @Variants			= ();
	my @VariantBases		= split (m/\,/, $VariantBases);
	
	# ($QD) 					= $Info =~ m/.*\;QD\=(.*?)\;.*/ if ($Info =~ m/.*\;QD\=.*\;.*/);
	
	 # if ($QD >= 3){
	
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
	# }
	
	##########################
	
	if 		($Caller =~	m/GATK/){
		
		my 	$Format 	= $LineValues[8];
		my 	$Sample 	= $LineValues[9];
		
		if ($Format eq "GT:AD:DP:GQ:PL"){
		
			$Sample =~ m/^.*:(.*):(.*):.*:.*$/;
			
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
		elsif ($Format eq "GT:GQ:PL"){
			next;
		}
		else {
		
			die ("Format does not match: $Line\n");
		}							
	}
	elsif 	($Caller eq "LoFreq"){
		
		($Depth, my $AFs, undef, undef) = split (/\;/, $Info);
		$Depth =~ s/DP=//g;
		$AFs =~ s/AF=//g;
		
		@AFs = ($AFs);
	}
	else {
	
		die ("Unknown caller\n");
	}
	
	return \@Variants, \@AFs, $Depth;
}

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