#######################################################################################
###						      Calculate Allelic Fractions SDs			    		###
#######################################################################################

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

#######################################################################################
###						      	   Initialize & declare		    					###
#######################################################################################

my $BCBIO_MAIN_DIR 					= "/user/data/gent/gvo000/gvo00082/vsc41234/bcbio/exome/";
my $BedFileGenesHg38				= "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/bed/ABCC6_ABCC6P1_ABCC6P2_Genes_HG38.bed";
my $BedFileGenes					= "/user/scratch/gent/gvo000/gvo00082/research/wst/abcc6/bed/ABCC6_ABCC6P1_ABCC6P2_Genes.bed";
my $POP_MAIN_DIR					= "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/population/";
my $POP_BAM_DIR						= $POP_MAIN_DIR . "Hg38_Bams/";
my $TEMP_DIR						= "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/temp/";
my $ReferenceSequence				= "/user/scratch/gent/gvo000/gvo00082/genomes/GRCh37/homo_sapiens_GRCh37.fa";
my $CreateVCF						= "no";
my $AnalyzeVCF						= "no";
my $CallerSelect					= "LoFreq";
my %Samples							= ();

#######################################################################################
###				 				  	 Get Options 			    		     		###
#######################################################################################

GetOptions ("CreateVCF=s"			=>	\$CreateVCF,
			"AnalyzeVCF=s"			=> 	\$AnalyzeVCF);

#######################################################################################
###						      		  Fetch BAM hg38		    					###
#######################################################################################

# system("cd $BCBIO_MAIN_DIR");
# print "find . -name *bam | while read FilePath; do samtools view -b -L $BedFileGenesHg38 \$FilePath > $POP_MAIN_DIR\${FilePath##*/} ;done\n";
# system("find . -name *bam | while read FilePath; do samtools view -b -L $BedFileGenesHg38 \$FilePath > $POP_MAIN_DIR\${FilePath##*/} ;done");

if ($CreateVCF eq "yes"){
	
	opendir(POP_DIR, $POP_BAM_DIR) or die ("Can't open $POP_BAM_DIR\n");
	while(my $BamFile = readdir(POP_DIR)){
		
		if ($BamFile =~ m/\.bam$/){
		
			(my	$Sample) 		= $BamFile =~ m/^(.*)\.bam$/g;
			 my $JobScript 		= $POP_BAM_DIR . $Sample . ".sh";
			 my $BamFilePath	= $POP_BAM_DIR . $BamFile;
			 
			next if (exists $Samples{$Sample});
			
			$Samples{$Sample} = undef;
			
			open SH, ">$JobScript" or die ("Can't open $JobScript\n");
			
			print SH "#PBS -N $Sample\n";
			print SH "#PBS -l nodes=1:ppn=1\n";
			print SH "#PBS -l walltime=01:00:00\n";
			print SH "#PBS -l mem=6GB\n";
			print SH "#PBS -d $POP_BAM_DIR\n";	

			print SH "module load Trimmomatic/0.32-Java-1.7.0_40\n";
			print SH "module load SAMtools/1.4-intel-2016b\n";	
			
			### Bam To Fastq ###
			
			print SH "samtools collate -uO $BamFilePath $TEMP_DIR$Sample | samtools fastq -1 $Sample\.abcc6.1.fastq -2 $Sample\.abcc6.2.fastq -0 $Sample\.abcc6.0.fastq -t -\n";
		
			print SH "gzip $Sample\.abcc6.1.fastq\n";
			print SH "gzip $Sample\.abcc6.2.fastq\n";
			
			### Trim ###
			
			print SH "java -jar \${EBROOTTRIMMOMATIC}\/trimmomatic-0.32.jar PE $Sample\.abcc6.1.fastq.gz $Sample\.abcc6.2.fastq.gz $Sample\.abcc6.q30.1.fastq.gz  $Sample\.abcc6.single.1.fastq.gz $Sample\.abcc6.q30.2.fastq.gz $Sample\.abcc6.single.2.fastq.gz LEADING:30 TRAILING:30\n";
			
			print SH "rm $Sample\.abcc6.single.1.fastq.gz\n";
			print SH "rm $Sample\.abcc6.single.2.fastq.gz\n";
			print SH "rm $Sample\.abcc6.0.fastq\n";	
			
			print SH "rm $Sample\.abcc6.1.fastq.gz\n";
			print SH "rm $Sample\.abcc6.2.fastq.gz\n";
			
			### Load New Modules ###
			
			print SH "module purge\n";
			print SH "module load BBMap/36.62-intel-2016b-Java-1.8.0_112\n";
			
			### Repair ###
			
			print SH "repair.sh -in1=$Sample\.abcc6.q30.1.fastq.gz in2=$Sample\.abcc6.q30.2.fastq.gz out1=$Sample\.abcc6.q30.1.rep.fastq.gz out2=$Sample\.abcc6.q30.2.rep.fastq.gz outsingle=$Sample\.abcc6.q30.0.rep.fastq.gz\n";
			
			print SH "rm $Sample\.abcc6.q30.1.fastq.gz\n";
			print SH "rm $Sample\.abcc6.q30.2.fastq.gz\n";
			
			print SH "module purge\n";
			print SH "module load BWA/0.7.15-intel-2016b\n";
			print SH "module load picard/2.1.1-Java-1.8.0_74\n";
			print SH "module load GATK/3.5-Java-1.8.0_74\n";
			print SH "module load LoFreq/2.1.2-intel-2016b-Python-2.7.12\n";
				
			### MASKING ###
			
			print SH "cd $POP_MAIN_DIR\/masking\/\n";
			
				# ### BWA ###
				
				print SH "bwa mem -t 1 -M /user/scratch/gent/gvo000/gvo00082/genomes/GRCh37_Masked_ABCC6P1P2/homo_sapiens_GRCh37_MaskedForABCC6P1P2.fa $POP_BAM_DIR$Sample\.abcc6.q30.1.rep.fastq.gz $POP_BAM_DIR$Sample\.abcc6.q30.2.rep.fastq.gz > $Sample\.abcc6.q30.sam 2>> $Sample\.abcc6.q30.log\n";
			
				print SH "samtools view -Sb $Sample\.abcc6.q30.sam > $Sample\.abcc6.q30.bam 2>> $Sample\.abcc6.q30.log\n";
				
				print SH "rm $Sample\.abcc6.q30.sam\n";
				
				print SH "java -jar \${EBROOTPICARD}\/picard.jar SortSam I=$Sample\.abcc6.q30.bam O=$Sample\.abcc6.q30.sorted.bam SO=coordinate 2>> $Sample\.abcc6.q30.log\n";
				
				print SH "rm $Sample\.abcc6.q30.bam\n";
				
				# ### REMDUP ###
				
				print SH "java -jar \${EBROOTPICARD}\/picard.jar MarkDuplicates I=$Sample\.abcc6.q30.sorted.bam O=$Sample\.abcc6.q30.sorted.remdup.bam M=$Sample\.abcc6.q30.sorted.remdup.met REMOVE_DUPLICATES=true 2>> $Sample\.abcc6.q30.log\n";
				
				print SH "rm $Sample\.abcc6.q30.sorted.bam\n";
				
				# ### READGROUPS ###
				
				print SH "java -jar \${EBROOTPICARD}\/picard.jar AddOrReplaceReadGroups I=$Sample\.abcc6.q30.sorted.remdup.bam O=$Sample\.abcc6.q30.sorted.remdup.rg.bam RGID=$Sample RGLB=lib1 RGPL=illumina RGPU=$Sample RGSM=$Sample\n";
				
				print SH "rm $Sample\.abcc6.q30.sorted.remdup.bam\n";
				
				# ### BQSR ###
				
				print SH "lofreq  indelqual --dindel -f $ReferenceSequence -o $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.bam $Sample\.abcc6.q30.sorted.remdup.rg.bam\n";
				
				print SH "rm $Sample\.abcc6.q30.sorted.remdup.rg.bam\n";
				
				# ### INDEX ###
				
				print SH "samtools index $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.bam\n";
				
				# ### UNIQUE ###
				
				print SH "samtools view -bh -q 10 $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.bam > $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.bam\n";
				
				print SH "rm $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.bam\n";
				
				# ### INDEX ###
				
				print SH "samtools index $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.bam\n";
				
				# ### CALL ###
				
				print SH "lofreq call --call-indels -f $ReferenceSequence -o $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.vcf $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.bam\n";
				
				print SH "java -jar \$EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ReferenceSequence -L \$VSC_SCRATCH_VO\/research/wst/abcc6/bed/ABCC6_ABCC6P1_ABCC6P2_Genes.bed -I $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.bam -o $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.gatk.vcf\n";
				
				### DEPTH ###
				
				print SH "samtools depth -b $BedFileGenes $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.bam > $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.cov\n";
							
			### NO MASKING ###
			
			print SH "cd $POP_MAIN_DIR\/nomasking\/\n";
			
				# ### BWA ###
				
				print SH "bwa mem -t 1 -M /user/scratch/gent/gvo000/gvo00082/genomes/GRCh37/homo_sapiens_GRCh37.fa $POP_BAM_DIR$Sample\.abcc6.q30.1.rep.fastq.gz $POP_BAM_DIR$Sample\.abcc6.q30.2.rep.fastq.gz > $Sample\.abcc6.q30.sam 2>> $Sample\.abcc6.q30.log\n";
				
				print SH "rm $POP_BAM_DIR$Sample\.abcc6.q30.1.rep.fastq.gz\n";
				print SH "rm $POP_BAM_DIR$Sample\.abcc6.q30.2.rep.fastq.gz\n";
			
				print SH "samtools view -Sb $Sample\.abcc6.q30.sam > $Sample\.abcc6.q30.bam 2>> $Sample\.abcc6.q30.log\n";
				
				print SH "rm $Sample\.abcc6.q30.sam\n";
				
				print SH "java -jar \${EBROOTPICARD}\/picard.jar SortSam I=$Sample\.abcc6.q30.bam O=$Sample\.abcc6.q30.sorted.bam SO=coordinate 2>> $Sample\.abcc6.q30.log\n";
				
				print SH "rm $Sample\.abcc6.q30.bam\n";
				
				# ### REMDUP ###
				
				print SH "java -jar \${EBROOTPICARD}\/picard.jar MarkDuplicates I=$Sample\.abcc6.q30.sorted.bam O=$Sample\.abcc6.q30.sorted.remdup.bam M=$Sample\.abcc6.q30.sorted.remdup.met REMOVE_DUPLICATES=true 2>> $Sample\.abcc6.q30.log\n";
				
				print SH "rm $Sample\.abcc6.q30.sorted.bam\n";
				
				# ### READGROUPS ###
				
				print SH "java -jar \${EBROOTPICARD}\/picard.jar AddOrReplaceReadGroups I=$Sample\.abcc6.q30.sorted.remdup.bam O=$Sample\.abcc6.q30.sorted.remdup.rg.bam RGID=$Sample RGLB=lib1 RGPL=illumina RGPU=$Sample RGSM=$Sample\n";
				
				print SH "rm $Sample\.abcc6.q30.sorted.remdup.bam\n";
				
				# ### BQSR ###
				
				print SH "lofreq  indelqual --dindel -f $ReferenceSequence -o $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.bam $Sample\.abcc6.q30.sorted.remdup.rg.bam\n";
				
				print SH "rm $Sample\.abcc6.q30.sorted.remdup.rg.bam\n";
				
				# ### INDEX ###
				
				print SH "samtools index $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.bam\n";
				
				# ### UNIQUE ###
				
				print SH "samtools view -bh -q 10 $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.bam > $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.bam\n";
				
				print SH "rm $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.bam\n";
				
				# ### INDEX ###
				
				print SH "samtools index $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.bam\n";
				
				# ### CALL ###
				
				print SH "lofreq call --call-indels -f $ReferenceSequence -o $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.vcf $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.bam\n";
				
				print SH "java -jar \$EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ReferenceSequence -L \$VSC_SCRATCH_VO\/research/wst/abcc6/bed/ABCC6_ABCC6P1_ABCC6P2_Genes.bed -I $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.bam -o $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.gatk.vcf\n";
				
				### DEPTH ###
				
				print SH "samtools depth -b $BedFileGenes $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.bam > $Sample\.abcc6.q30.sorted.remdup.rg.bqsr.unique.cov\n";
			
			close SH;
			
			system("qsub $JobScript");
		}
	}
	closedir(POP_DIR);
}

elsif ($AnalyzeVCF eq "yes"){

	###################################################################################
	###                		  		 initialize and declare         				###
	###################################################################################
	
	my $ClustalO_FP			= "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/CLustalO";
	my $OutputDir 			= "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/outputscripts/";
	
	
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
	
	my %SDsAll				= ();
	
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
	###                		 	 		read in coverage          					###
	###################################################################################
	
	print "Read in coverage\n";
	
	my $MaskingDir 		= $POP_MAIN_DIR . "masking/";
	my $NoMaskingDir 	= $POP_MAIN_DIR . "nomasking/";
	my $Analysis		= "/";
	
	my @AnalysisDirs 	= ($MaskingDir, $NoMaskingDir);
	
	foreach my $AnalysisDir (@AnalysisDirs){
	
		if 		($AnalysisDir =~ m/nomasking/){
		
			$Analysis = "nomasking";
		}
		elsif 	($AnalysisDir =~ m/masking/){
		
			$Analysis = "masking";
		}
		
		my %Coverage			= ();
		my %SD_Variants			= ();
		my %AllVariants			= ();
		my %SD_Positions		= ();
		my %AllPositions		= ();
		my %VarsToExclude		= ();
		my %Samples				= ();
	
		opendir(AN_DIR, $AnalysisDir) or die ("Can't open $AnalysisDir\n");
		while(my $CovFile = readdir(AN_DIR)){
		
			if ($CovFile =~ m/^(.*)\.abcc6\.q30\.sorted\.remdup\.rg\.bqsr\.unique\.cov$/){
			
				$CovFile =~ m/^(.*)\.abcc6\.q30\.sorted\.remdup\.rg\.bqsr\.unique\.cov$/;
				
				my $Sample = $1;
				
				$Samples{$Sample} = undef;
				
				my $CovFilePath = $AnalysisDir . $CovFile;
			
				open COV, "$CovFilePath" or die ("Can't open $CovFile\n");
				while (<COV>){
				
					 my $Line = $_;
						$Line =~ s/\n//g;
						
					(undef, my $Position, my $Cov) = split (m/\t/, $Line);
					
					$Coverage{$Sample}{$Position} = $Cov;		
				}
				close COV;
			}
		}
		closedir(AN_DIR);
		
		###################################################################################
		###                		 	 		 loop over vcf          					###
		###################################################################################
		
		print "loop over vcf\n";
		
		opendir(AN_DIR, $AnalysisDir) or die ("Can't open $AnalysisDir\n");
		while(my $VcfFile = readdir(AN_DIR)){
		
			if ($VcfFile =~ m/^(.*)\.abcc6\.q30\.sorted\.remdup\.rg\.bqsr\.unique\.vcf$/){
			
				$VcfFile =~ m/^(.*)\.abcc6\.q30\.sorted\.remdup\.rg\.bqsr\.unique\.vcf$/;
				
				my $Sample = $1;
				
				my $VcfFilePath = $AnalysisDir . $VcfFile;
			
				open VCF, "$VcfFilePath" or die ("Can't open $VcfFilePath\n");
				while (<VCF>){
				
					next if ($_ =~ m/^#/);
				
					 my $Line = $_;
						$Line =~ s/\n//g;
											
					(my $VariantsRef, 
					 my $AFsRef, 
					 my $Depth) = DetermineVariantInVCF_Line($Line, $CallerSelect);	

					
					for (my $I = 0; $I < scalar @$VariantsRef; $I++){
					
						my $Variant = $VariantsRef->[$I];
						my $AF		= $AFsRef->[$I];
						
						if ($AF	ne "/"){

							if ($SDsInAmp{$Variant} || $SDsAll{$Variant}){
						
								$SD_Variants{"exome"}{$Analysis}{"unique"}{$CallerSelect}{$Variant}{$Sample} = "$AF\_$Depth";
								
								if ($AF >= 0.05){
									
									(undef, my $Pos, undef, undef) = split (m/\_/, $Variant);
									
									$SD_Positions{$Pos} = $Variant;
								}
							}
							else {
								
								$AllVariants{"exome"}{$Analysis}{"unique"}{$CallerSelect}{$Variant}{$Sample} = "$AF\_$Depth";
								
								if ($AF >= 0.05){
								
									(undef, my $Pos, undef, undef) = split (m/\_/, $Variant);
									
									$AllPositions{$Pos} = $Variant;
								}
							}
						}
					}
				}
				close VCF;
			}
		}
		closedir(AN_DIR);
		
		###################################################################################
		###                	 check for what variants there are > 2 alleles          	###
		###################################################################################
		
		print "check for what variants there are > 2 alleles\n";
		
		foreach my $SD_Pos (sort keys %SD_Positions){
		
			if ($AllPositions{$SD_Pos}){
				
				print $AllPositions{$SD_Pos} . "\n";
				
				foreach my $Sample (keys %{$AllVariants{"exome"}{$Analysis}{"unique"}{$CallerSelect}{$AllPositions{$SD_Pos}}}){
					
					# print "\t" . $Sample . "\t" . $AllVariants{"exome"}{$Analysis}{"unique"}{$CallerSelect}{$AllPositions{$SD_Pos}}{$Sample} . "\n";
				
				}
				
				# $VarsToExclude{$SD_Positions{$SD_Pos}} = undef;
			}
		}
		
		###################################################################################
		###                	 		    Write Out Variant Matrix         				###
		###################################################################################
		
		print "Write Out Variant Matrix \n";
		
		my $VariantMatrixFP = $OutputDir . "$Analysis\.VariantMatrix.txt";
		
		open MAT, ">$VariantMatrixFP" or die ("Can't open $VariantMatrixFP\n");
		
		foreach my $VarPos (sort { $a <=> $b } keys %SD_Positions){
		
			my $Variant = $SD_Positions{$VarPos};
			
			next if (exists $VarsToExclude{$Variant});
			
			# print MAT $Variant . "\t";
			
			foreach my $Sample (sort keys %Samples){			
				
				if (exists $SD_Variants{"exome"}{$Analysis}{"unique"}{$CallerSelect}{$Variant} && 
					exists $SD_Variants{"exome"}{$Analysis}{"unique"}{$CallerSelect}{$Variant}{$Sample}){
				
					(my $AF, my $Depth) = split (m/\_/, $SD_Variants{"exome"}{$Analysis}{"unique"}{$CallerSelect}{$Variant}{$Sample});
											
					if ($Depth >= 60){
						
						print MAT $Sample . "\t" . $Variant . "\t" . $AF . "\n";
					}
					else {
						
						print MAT $Sample . "\t" . $Variant . "\t" . "NA" . "\n"; 
					}
				}
				elsif (exists $Coverage{$Sample}{$VarPos} && $Coverage{$Sample}{$VarPos} >= 60){

					print MAT $Sample . "\t" . $Variant . "\t" . "0" . "\n"; 
				}
				else {
					
					print MAT $Sample . "\t" . $Variant . "\t" . "NA" . "\n"; 
				}
			}
			
			# print MAT "\n";
		}
		
		close MAT;
		
		
		###################################################################################
		###                	 		    Write Out Genotype Counts         				###
		###################################################################################
		
		print "Write Out Genotype Counts \n";
		
		my $GenotypeFile1_4 = $OutputDir . "$Analysis\.GenotypesExon1-4.txt";
		my $GenotypeFile5_9 = $OutputDir . "$Analysis\.GenotypesExon5-9.txt";
		
		open _1_4, ">$GenotypeFile1_4" or die ("Can't open $GenotypeFile1_4\n");
		open _5_9, ">$GenotypeFile5_9" or die ("Can't open $GenotypeFile5_9\n");
		print _1_4 "Variant\tHomAltFreq\tHomRefFreq\tOneHetCount\tNrOfSamples\tLowCov\n";
		print _5_9 "Variant\tHomAltFreq\tHomRefFreq\tNrOfSamples\tLowCov\n";
			
		foreach my $VarPos (sort { $a <=> $b } keys %SD_Positions){
		
			my $Variant = $SD_Positions{$VarPos};
			
			next if (exists $VarsToExclude{$Variant});
			
			my $HomRefCount = 0;
			my $HomAltCount = 0;
			my $OneHetCount = 0;
			my $NrOfSamples	= 0;
			my $LowCov		= 0;
					
			if ($VarPos <= 16318328 && $VarPos >= 16312404){
			
				###################################################################################
				###                	 		determine frequencies exon 1 - 4          			###
				###################################################################################
				
				foreach my $Sample (keys %Samples){			
				
					if (exists $SD_Variants{"exome"}{$Analysis}{"unique"}{$CallerSelect}{$Variant} && 
						exists $SD_Variants{"exome"}{$Analysis}{"unique"}{$CallerSelect}{$Variant}{$Sample}){
					
						(my $AF, my $Depth) = split (m/\_/, $SD_Variants{"exome"}{$Analysis}{"unique"}{$CallerSelect}{$Variant}{$Sample});
												
						if ($Depth >= 60){
							
							if 		($AF >= 1100/12){ # (500/6 + 600/6)/2 # hom alt in ABCC6
													
								$HomAltCount++;
								$NrOfSamples++;
							}
							elsif ($AF >= 900/12){ # most likely het in ABCC6
								
								$OneHetCount++;
								$NrOfSamples++;
							}
							elsif 	($AF <= 100/12){
							
								$HomRefCount++;
								$NrOfSamples++;
							}
							else {
								
								$NrOfSamples++;
							}
						}
						else {
							
							$LowCov++;
						}
					}
					elsif (exists $Coverage{$Sample}{$VarPos} && $Coverage{$Sample}{$VarPos} >= 60){

						$HomRefCount++;
						$NrOfSamples++;
					}
					elsif (exists $Coverage{$Sample}{$VarPos}){
						
						$LowCov++;
					}
				}
				
				if ($NrOfSamples && $NrOfSamples >=20){
					
					my $HomAltFreq = $HomAltCount/$NrOfSamples;
					my $HomRefFreq = $HomRefCount/$NrOfSamples;
					
					print _1_4 "$Variant\t$HomAltFreq\t$HomRefFreq\t$OneHetCount\t$NrOfSamples\t$LowCov\n";
				}
			}	
			
			else {
			
				###################################################################################
				###                	 		determine frequencies exon 5 - 9          			###
				###################################################################################
				
				foreach my $Sample (keys %Samples){			
				
					if (exists $SD_Variants{"exome"}{$Analysis}{"unique"}{$CallerSelect}{$Variant} && 
						exists $SD_Variants{"exome"}{$Analysis}{"unique"}{$CallerSelect}{$Variant}{$Sample}){
					
						(my $AF, my $Depth) = split (m/\_/, $SD_Variants{"exome"}{$Analysis}{"unique"}{$CallerSelect}{$Variant}{$Sample});
						
						if ($Depth >= 20){
							
							if 		($AF >= ((3/4+4/4)/2)){
													
								$HomAltCount++;
								$NrOfSamples++;
							}
							elsif 	($AF <= 1/8){
							
								$HomRefCount++;
								$NrOfSamples++;
							}
							else {
								
								$NrOfSamples++;
							}
						}
						else {
							
							$LowCov++;
						}
					}
					elsif (exists $Coverage{$Sample}{$VarPos} && $Coverage{$Sample}{$VarPos} >= 20){

						$HomRefCount++;
						$NrOfSamples++;
					}
					elsif (exists $Coverage{$Sample}{$VarPos}){
						
						$LowCov++;
					}
				}
				
				if ($NrOfSamples && $NrOfSamples >=20){
					
					my $HomAltFreq = $HomAltCount/$NrOfSamples;
					my $HomRefFreq = $HomRefCount/$NrOfSamples;
					
					print _5_9 "$Variant\t$HomAltFreq\t$HomRefFreq\t$NrOfSamples\t$LowCov\n";
				}	
			}
		}
		
		close _1_4;
		close _5_9;
	}
}

###################################################################################
###                		 	 		   subroutines          					###
###################################################################################

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