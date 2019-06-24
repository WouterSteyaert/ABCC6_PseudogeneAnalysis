#######################################################################################
###						      Calculate Allelic Fractions SDs			    		###
#######################################################################################

use strict;
use warnings;
use Getopt::Long;

###################################################################################
### 			    		   INITIALIZE & DECLARE 							###
###################################################################################

my $OutputDir 				= "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/outputscripts/";
my $GSEA_SpecFilePath		= $OutputDir . "fordBWsamples.nomasking.unique.exoncoverage.norm";
my $WES_AF_FilePath			= $OutputDir . "exomeBWsamples.masking.all.LoFreq.SDsAmp";
my %CompareDatFormat		= ();
my %SDs						= ();

###################################################################################
### 			    		   	   Loop over files 								###
###################################################################################

my $GSEA_SpecMatrix = ReadInFileWithHeadersAsAMatrix($GSEA_SpecFilePath, "\t");
my $WES_AF_Matrix	= ReadInFileWithHeadersAsAMatrix($WES_AF_FilePath, "\t");

foreach my $LineNr (sort {$a <=> $b} keys %$GSEA_SpecMatrix){

	if ($LineNr >= 9 && $LineNr <= 19){
		
		my $ExonRank 	= $GSEA_SpecMatrix->{$LineNr}{"Gene~Exon"};
		my $Gene 		= "";
		
		if ($LineNr < 17){
			$Gene = "ABCC6P1";
		}
		else {
			$Gene = "ABCC6P2";
		}
		
		foreach my $Column (keys %{$GSEA_SpecMatrix->{$LineNr}}){
		
			if ($Column =~ m/^D.*/){
			
				$CompareDatFormat{$Gene}{$ExonRank}{"$Column\-1"}{"GSEA"} = $GSEA_SpecMatrix->{$LineNr}{$Column};
			}
		}	
	}
}

foreach my $LineNr (sort {$a <=> $b} keys %$WES_AF_Matrix){
		
	my $Variant_Exon_Pi 	= $WES_AF_Matrix->{$LineNr}{"Variant~Exon~Pi"};
	
	my @Variant_Exon_Pi_Values = split (m/\~/, $Variant_Exon_Pi);
	
	my $Variant = shift (@Variant_Exon_Pi_Values);
	my $Exon 	= shift (@Variant_Exon_Pi_Values);
	
	if ($Exon == 3 || $Exon == 4){
		$Exon = "3~4";
	}
	
	next if (scalar @Variant_Exon_Pi_Values == 2);
	
	my $GeneAbr = shift (@Variant_Exon_Pi_Values);
	my $Gene 	= "ABCC6" . $GeneAbr;
	my $DupNr	= 0;
	
	if (scalar @Variant_Exon_Pi_Values != 2){
	
		$SDs{$Gene}{$Exon}{$Variant} = undef;
	
		$DupNr = scalar (keys %{$SDs{$Gene}{$Exon}});
	
	}
	
	foreach my $Column (keys %{$WES_AF_Matrix->{$LineNr}}){
	
		if ($Column =~ m/^D.*/){
		
			my $Value = $WES_AF_Matrix->{$LineNr}{$Column};
			
			# if ($Value =~ m/\d+/){
				# $Value -= 0.5;
			# }
					
			if (scalar @Variant_Exon_Pi_Values == 2){
				
				# $CompareDatFormat{"ABCC6P1"}{$Exon}{$Column}{"WES"}{$Variant} = $Value;
				# $CompareDatFormat{"ABCC6P2"}{$Exon}{$Column}{"WES"}{$Variant} = $Value;
			}
			else {
			
				$CompareDatFormat{$Gene}{$Exon}{"$Column\-$DupNr"}{"WES"} = $Value;
				$CompareDatFormat{$Gene}{$Exon}{"$Column\-$DupNr"}{"GSEA"} = $CompareDatFormat{$Gene}{$Exon}{"$Column\-1"}{"GSEA"};
			}
		}
	}	
}


###################################################################################
### 			    		   	   	  Write Out									###
###################################################################################

my @Exons = ("ABCC6P1_1", "ABCC6P1_2","ABCC6P1_3~4","ABCC6P1_5","ABCC6P1_6","ABCC6P1_7","ABCC6P1_8","ABCC6P1_9","ABCC6P2_1","ABCC6P2_2","ABCC6P2_3~4");

foreach my $Exon (@Exons){

	(my $Gene, my $ExonRank) = split (m/_/, $Exon);

	my $ComparFilePath = $OutputDir . "ComparisonGSEA_SpecWithWES_AF\.$Exon";
	
	open COMP, ">$ComparFilePath" or die ("Can't open $ComparFilePath\n");
	print COMP "DNA\tType\tRatio\n";
	
	foreach my $DNA (keys %{$CompareDatFormat{$Gene}{$ExonRank}}){
		foreach my $Enrichment (keys %{$CompareDatFormat{$Gene}{$ExonRank}{$DNA}}){
			
			
			print COMP "$DNA\t$Enrichment\t$CompareDatFormat{$Gene}{$ExonRank}{$DNA}{$Enrichment}\n";
			
		}
	}
	
	close COMP;
}





###################################################################################
### 			    		   	   	  Subroutines 								###
###################################################################################

#######################################################
# 		  	 ReadInFileWithHeadersAsAMatrix 	   	  #
#######################################################

sub ReadInFileWithHeadersAsAMatrix {
	
	my $FilePath 				= $_[0];
	my $Delimiter				= $_[1];
	my $LineNr					= 0;
	my %VariantFileAsAMatrix 	= ();
	my %Headers					= ();
	
	open FILE, "$FilePath" or die ("Can't open $FilePath\n");
	while (<FILE>){
		
		$_ =~ s/\n//;
		$_ =~ s/\r//;
		
		if ($_ ne ""){
		
			if ($LineNr == 0){
				
				my 	@Headers 	= split($Delimiter, $_);
				my	$ColumnNr	= 0;
								
				foreach my $Header (@Headers){
					
					$Header =~ s/ //g;
					
					$Headers{$ColumnNr} = StripQuotes(StripBracket($Header));
					$ColumnNr++;
				}
			}
			else {
				
				my @LineValues 	= split ($Delimiter, $_);
				my $ColumnNr	= 0;
				
				foreach my $LineValue (@LineValues){
					
					$LineValue = StripQuotes(StripBracket($LineValue));
					
					if ($LineValue ne ""){
						$VariantFileAsAMatrix{$LineNr}{$Headers{$ColumnNr}} = $LineValue;						
					}
					else {
						$VariantFileAsAMatrix{$LineNr}{$Headers{$ColumnNr}} = "/";
					}
					$ColumnNr++;	
				}
			}
			$LineNr++;
		}
	}
	close FILE;

	return \%VariantFileAsAMatrix;
}

#######################################################
# 	  			   	  StripBracket 	   		  	  	  #
#######################################################

sub StripBracket{
	
	my $string = shift;
	$string =~ s/[\[\]\(\)]//g;
	return $string;
}

#######################################################
# 	  			   	  StripQuotes 	   		  	  	  #
#######################################################

sub StripQuotes{
	
	my $string = shift;
	$string =~ s/\"//g;
	return $string;
}