#######################################################################################
###     		     	Merge, Reduce, Filter BAM and Call Variants   	 			###
#######################################################################################
#######################################################################################
###     		     			  	  Load Libraries  	  		   					###
#######################################################################################

use strict;
use warnings;

#######################################################################################
### 	  	 				   	  INITIALIZE AND DECLARE    	 					###
#######################################################################################

my @AnalysisDirs 				= ( 	
									"/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/fordValidationSeq2/masking/analysis/",
									"/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/fordValidationSeq2/nomasking/analysis/"
								  );

my $ReferenceSequence			= 	"/kyukon/scratch/gent/vo/000/gvo00082/genomes/GRCh37/homo_sapiens_GRCh37.fa";
my $BedFileRegions				=	"/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/bed/ABCC6_ABCC6P1_ABCC6P2_Regions.bed";
my $BedFileGenes				= 	"/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/bed/ABCC6_ABCC6P1_ABCC6P2_Genes.bed";
my $BedFileExons				= 	"/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/bed/ABCC6_ABCC6P1_ABCC6P2_Exons.bed";
my $NonPseudoRegionsFilePath 	= 	"/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/bed/NonPseudo.bed";
my $TempDir						= 	"/kyukon/scratch/gent/vo/000/gvo00082/research/wst/temp/";

my %BamFiles					= ();
my %Regions						= ();

#######################################################################################
### 	  	 				   	     Read in BED file    	 						###
#######################################################################################

open BED, "$BedFileRegions" or die ("Can't open $BedFileRegions\n");
while (<BED>){
	
	 my	$Line = $_;
		$Line =~ s/\n//g;
		
	(my $Chromosome, my $Start, my $End, my $Region) = split ("\t", $Line);
		
		$Region = lc($Region);
		$Region =~ s/region_//g;
		$Region =~ s/_/\./g;
		$Region =~ s/ex/r/g;
	 
	 $Regions{$Region} = "$Chromosome\t$Start\t$End\t$Region";
}
close BED;

#######################################################################################
### 	  	 				   	 		  LOOP    	 								###
#######################################################################################

foreach my $AnalysisDir (@AnalysisDirs){

	print $AnalysisDir . "\n";
	
	opendir (ANADIR, $AnalysisDir) or die ("Can't open $AnalysisDir!\n");
	while (my $BamFile = readdir(ANADIR)){
	
		my $Sample	  		= "/";
	
		if 		($AnalysisDir =~ m/exome/){
		
			next if ($BamFile !~  m/^(.*)_S\d+_L\d+\.q30\.sorted\.remdup\.rg\.bam$/);
		
			$BamFile =~ m/^(.*)_S\d+_L\d+\.q30\.sorted\.remdup\.rg\.bam$/;
			
			$Sample	 = $1;
			
		}
		elsif 	($AnalysisDir =~ m/ford/){
		
			next if ($BamFile !~  m/^(.*)\.q30\.sorted\.remdup\.rg\.bam$/);
		
			$BamFile =~ m/^(.*)\.q30\.sorted\.remdup\.rg\.bam$/;
			
			$Sample	 = $1;
			
		}
		else {
		
			die ("No Ford & Exome\n");
		}
		
		my $BamFilePath   	= $AnalysisDir . $BamFile;
		
		print $Sample . "\n";

		if (defined $BamFiles{$Sample}){
			$BamFiles{$Sample} .= " $BamFilePath";
		}
		else {
			$BamFiles{$Sample} = "$BamFilePath";
		}
	}
	closedir(ANADIR);

	foreach my $Sample (keys %BamFiles){

        my 	$JobScript           				= $AnalysisDir . $Sample . ".sh";
		
		print $JobScript . "\n";
		
   		open JOB, ">$JobScript" or die ("Can't open $JobScript\n");

		print JOB "#PBS -N $Sample\n";
        print JOB "#PBS -l nodes=1:ppn=1\n";
        print JOB "#PBS -l walltime=09:00:00\n";
		print JOB "#PBS -l mem=8gb\n";
		print JOB "#PBS -d $AnalysisDir\n";		
		
		print JOB "module load BamUtil/1.0.13-intel-2016b\n";
        print JOB "module load LoFreq/2.1.2-intel-2016b-Python-2.7.12\n";
        print JOB "module load picard/2.1.1-Java-1.8.0_74\n";
        print JOB "module load GATK/3.5-Java-1.8.0_74\n";
		
		if 	  ($AnalysisDir =~ m/exome/){
		
			print JOB "samtools merge -f $Sample\.bam $BamFiles{$Sample}\n";
		}
		elsif ($AnalysisDir =~ m/ford/) {
			
			print JOB "mv $BamFiles{$Sample} $Sample\.bam\n";
		}
		else {
			die ($!);
		}
       	
		my @BamsPerLane = split(m/ /, $BamFiles{$Sample});
		foreach my $BamPerLane (@BamsPerLane){
		
			# print JOB "rm $BamPerLane\n";
		}
		
		print JOB "samtools index $Sample\.bam\n";
		
		print JOB "java -jar \${EBROOTPICARD}\/picard.jar SortSam I=$Sample\.bam O=$Sample\.sorted.bam SO=coordinate TMP_DIR=$TempDir 2>> $Sample\.log\n";
		
		# print JOB "rm $Sample\.bam\n";
		
		print JOB "samtools view -b -L $BedFileGenes $Sample\.sorted.bam > $Sample\.sorted.abcc6.bam\n";
		print JOB "samtools view -b -L $NonPseudoRegionsFilePath $Sample\.sorted.bam > $Sample\.sorted.nonpseudo.bam\n";
		
		print JOB "samtools index $Sample\.sorted.abcc6.bam\n";
		print JOB "samtools index $Sample\.sorted.nonpseudo.bam\n";
		
		print JOB "bam clipOverlap --in $Sample\.sorted.abcc6.bam --out $Sample\.sorted.abcc6.clipped.bam 2> /dev/null\n";
		print JOB "bam clipOverlap --in $Sample\.sorted.nonpseudo.bam --out $Sample\.sorted.nonpseudo.clipped.bam 2> /dev/null\n";
		
		print JOB "samtools index $Sample\.sorted.abcc6.clipped.bam\n";
		print JOB "samtools index $Sample\.sorted.nonpseudo.clipped.bam\n";
		
		print JOB "samtools idxstats $Sample\.sorted.abcc6.clipped.bam > $Sample\.sorted.nonpseudo.clipped.stat\n";
		print JOB "samtools idxstats $Sample\.sorted.nonpseudo.clipped.bam > $Sample\.sorted.nonpseudo.clipped.stat\n";
		
		print JOB "samtools depth -a -b $BedFileGenes $Sample\.sorted.abcc6.clipped.bam > $Sample\.sorted.abcc6.clipped.cov\n";
		
		print JOB "samtools view -bh -q 10 $Sample\.sorted.abcc6.clipped.bam > $Sample\.sorted.abcc6.clipped.unique.bam\n";
		print JOB "samtools view -bh -q 10 $Sample\.sorted.nonpseudo.clipped.bam > $Sample\.sorted.nonpseudo.clipped.unique.bam\n";
		
		print JOB "samtools index $Sample\.sorted.abcc6.clipped.unique.bam\n";
		print JOB "samtools index $Sample\.sorted.nonpseudo.clipped.unique.bam\n";
		
		print JOB "samtools depth -a -b $BedFileGenes $Sample\.sorted.abcc6.clipped.unique.bam > $Sample\.sorted.abcc6.clipped.unique.cov\n";
		print JOB "samtools depth -a -b $NonPseudoRegionsFilePath $Sample\.sorted.nonpseudo.clipped.unique.bam > $Sample\.sorted.nonpseudo.clipped.unique.cov\n";
		
		# only BQSR for ABCC6 because this is only needed for variant calling #
		
		print JOB "lofreq  indelqual --dindel -f $ReferenceSequence -o $Sample\.sorted.abcc6.bqsr.bam $Sample\.sorted.abcc6.bam\n";
		
		print JOB "samtools index $Sample\.sorted.abcc6.bqsr.bam\n";
		
		print JOB "samtools view -bh -q 10 $Sample\.sorted.abcc6.bqsr.bam > $Sample\.sorted.abcc6.bqsr.unique.bam\n";
		
		print JOB "samtools index $Sample\.sorted.abcc6.bqsr.unique.bam\n";
		
		# LoFreq #
		
		print JOB "lofreq call-parallel --call-indels --pp-threads 4 -f $ReferenceSequence -o $Sample\.sorted.abcc6.bqsr.unique.vcf $Sample\.sorted.abcc6.bqsr.unique.bam\n";
		
		print JOB "lofreq call-parallel --call-indels  --pp-threads 4 -f $ReferenceSequence -o $Sample\.sorted.abcc6.bqsr.vcf $Sample\.sorted.abcc6.bqsr.bam\n";	
		
		# GATK #
		
		print JOB "java -jar \$EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ReferenceSequence -L \$VSC_SCRATCH_VO\/research/wst/abcc6/bed/ABCC6_ABCC6P1_ABCC6P2_Genes.bed -I $Sample\.sorted.abcc6.bqsr.unique.bam -o $Sample\.sorted.abcc6.bqsr.unique.gatk35.default.vcf\n";

		print JOB "java -jar \$EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ReferenceSequence -L \$VSC_SCRATCH_VO\/research/wst/abcc6/bed/ABCC6_ABCC6P1_ABCC6P2_Genes.bed -mmq 0 -I $Sample\.sorted.abcc6.bqsr.unique.bam -o $Sample\.sorted.abcc6.bqsr.unique.gatk35.mmq0.vcf\n";
		
		print JOB "java -jar \$EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ReferenceSequence -L \$VSC_SCRATCH_VO\/research/wst/abcc6/bed/ABCC6_ABCC6P1_ABCC6P2_Genes.bed -I $Sample\.sorted.abcc6.bqsr.bam -o $Sample\.sorted.abcc6.bqsr.gatk35.default.vcf\n";

		print JOB "java -jar \$EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ReferenceSequence -L \$VSC_SCRATCH_VO\/research/wst/abcc6/bed/ABCC6_ABCC6P1_ABCC6P2_Genes.bed -mmq 0 -I $Sample\.sorted.abcc6.bqsr.bam -o $Sample\.sorted.abcc6.bqsr.gatk35.mmq0.vcf\n";
		
		foreach my $Region (keys %Regions){
			
			my 	$RegionFilePath 		= "/user/scratch/gent/gvo000/gvo00082/research/wst/abcc6/bed/$Region\.bed";
			
			unless (-f $RegionFilePath){
				
				open REG, ">$RegionFilePath" or die $!;
				print REG "$Regions{$Region}\n";
				close REG;
			}
						
			print JOB "samtools view -b -L $RegionFilePath $Sample\.sorted.abcc6.clipped.bam > $Sample\.sorted.abcc6.clipped.$Region\.bam\n";
			
			print JOB "samtools index  $Sample\.sorted.abcc6.clipped.$Region\.bam\n";
			
			print JOB "samtools idxstats  $Sample\.sorted.abcc6.clipped.$Region\.bam > $Sample\.sorted.abcc6.clipped.$Region\.stat\n";
			
			print JOB "samtools view -bh -q 10 $Sample\.sorted.abcc6.clipped.$Region\.bam > $Sample\.sorted.abcc6.clipped.$Region\.unique.bam\n";
			
			print JOB "samtools index $Sample\.sorted.abcc6.clipped.$Region\.unique.bam\n";	
			
			print JOB "samtools idxstats $Sample\.sorted.abcc6.clipped.$Region\.unique.bam >   $Sample\.sorted.abcc6.clipped.$Region\.unique.stat\n";
		}
		
		close JOB;

		system ("qsub $JobScript");
	}
	
	%BamFiles = ();
}
 

