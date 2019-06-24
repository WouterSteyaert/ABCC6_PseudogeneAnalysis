#######################################################################################
###     		     			    	   MAP 										###
#######################################################################################
#######################################################################################
###     		     			  	  Load Libraries 								###
#######################################################################################

use strict; 
use warnings; 
use Getopt::Long;
		
#######################################################################################
### 	  	 				   	  INITIALIZE AND DECLARE 							###
#######################################################################################

my $FastqDir 			= "/user/scratch/gent/gvo000/gvo00082/research/wst/abcc6/exomeBWsamples/fastqs/";
my $TrimmedFastqDir     = "/user/scratch/gent/gvo000/gvo00082/research/wst/abcc6/exomeBWsamples/trimmedfastqs/";
my $AnalysisDir_1       = "/user/scratch/gent/gvo000/gvo00082/research/wst/abcc6/exomeBWsamples/masking/analysis/";
my $AnalysisDir_2       = "/user/scratch/gent/gvo000/gvo00082/research/wst/abcc6/exomeBWsamples/nomasking/analysis/";

#######################################################################################
### 	  	 				    	  	LOOP 										###
#######################################################################################

opendir (FASTQDIR, $FastqDir) or die ("Can't open $FastqDir!\n"); 
while (my $FirstFastQGZFile = readdir(FASTQDIR)){
	
	if ($FirstFastQGZFile =~ m/.*_R1_.*\.fastq\.gz$/){

		# next if ($FirstFastQGZFile =~ m/^Undeter/);
		# next if ($FirstFastQGZFile =~ m/^D1010417/);
		# next if ($FirstFastQGZFile =~ m/^D1509853/);
		# next if ($FirstFastQGZFile =~ m/^D1009445/);
		
		 my $FirstFilePath					= $FastqDir . $FirstFastQGZFile;
		 my $SecondFastQFile 				= $FirstFastQGZFile;
		    $SecondFastQFile 				=~ s/_R1/_R2/g;
		 my $SecondFilePath 				= $FastqDir . $SecondFastQFile;
		(my $FilePathBasis_1) 				= $FirstFastQGZFile =~ m/^(.*)_R1_001\.fastq\.gz$/;
		 my $FilePathBasis_2				= $FilePathBasis_1;
		    $FilePathBasis_1				= $AnalysisDir_1 . $FilePathBasis_1;
		    $FilePathBasis_2            	= $AnalysisDir_2 . $FilePathBasis_2;
		 		 		
		 my $FirstTrimmedPairedFilePath		=  $TrimmedFastqDir . $FirstFastQGZFile;
		    $FirstTrimmedPairedFilePath		=~ s/\.fastq\.gz$/\.q30\.paired\.fastq\.gz/g;
		 my $FirstTrimmedUnpairedFilePath 	=  $FirstTrimmedPairedFilePath;
		    $FirstTrimmedUnpairedFilePath   =~ s/\.q30\.paired\.fastq\.gz$/\.q30\.unpaired\.fastq\.gz/g;

		 my $SecondTrimmedPairedFilePath	=  $TrimmedFastqDir . $SecondFastQFile;
		    $SecondTrimmedPairedFilePath	=~ s/\.fastq\.gz$/\.q30\.paired\.fastq\.gz/g;
		 my $SecondTrimmedUnPairedFilePath	=  $SecondTrimmedPairedFilePath;
		    $SecondTrimmedUnPairedFilePath  =~ s/\.q30\.paired\.fastq\.gz$/\.q30\.unpaired\.fastq\.gz/g;

		 my $ShellScript                	= $FirstTrimmedPairedFilePath;
		    $ShellScript					=~ s/\.q30\.paired\.fastq\.gz/\.sh/g;
			
		$FirstFilePath =~ m/^(.*)\_S\d+\_L(\d+)\_.*gz$/;
		
		my 	$Sample 	= $1;
		my 	$RG 		= $2;
			$RG			=~ s/^0+//g;
	
		open JOBSCRIPT, ">$ShellScript" or die ("Can't open ShellScript\n");
		
		print JOBSCRIPT "#PBS -N $FirstFastQGZFile\n";
		print JOBSCRIPT	"#PBS -l nodes=1:ppn=8\n";
		print JOBSCRIPT	"#PBS -l walltime=09:00:00\n";
		print JOBSCRIPT	"#PBS -d $TrimmedFastqDir\n";
		
		print JOBSCRIPT "module load Trimmomatic/0.32-Java-1.7.0_40\n";
		
		print JOBSCRIPT "java -jar \${EBROOTTRIMMOMATIC}\/trimmomatic-0.32.jar PE $FirstFilePath $SecondFilePath $FirstTrimmedPairedFilePath  $FirstTrimmedUnpairedFilePath $SecondTrimmedPairedFilePath $SecondTrimmedUnPairedFilePath LEADING:30 TRAILING:30\n";
		
		# print JOBSCRIPT "rm $FirstFilePath\n";
		# print JOBSCRIPT "rm $SecondFilePath\n";
		
		print JOBSCRIPT "module unload Java/1.7.0_40\n";		
		print JOBSCRIPT	"module load BWA/0.7.13-intel-2016a\n";
		print JOBSCRIPT	"module load SAMtools/1.3-intel-2016a\n";
		print JOBSCRIPT	"module load Java/1.8.0_74\n";
		print JOBSCRIPT	"module load picard/2.1.1-Java-1.8.0_74\n";
		
		# masking #
		
		print JOBSCRIPT "bwa mem -t 8 -M /user/scratch/gent/gvo000/gvo00082/genomes/GRCh37_Masked_ABCC6P1P2/homo_sapiens_GRCh37_MaskedForABCC6P1P2.fa $FirstTrimmedPairedFilePath $SecondTrimmedPairedFilePath > $FilePathBasis_1\.q30.sam 2>> $FilePathBasis_1\.q30.log\n";
		
		print JOBSCRIPT "samtools view -Sb $FilePathBasis_1\.q30.sam > $FilePathBasis_1\.q30.bam 2>> $FilePathBasis_1\.q30.log\n";
		
		print JOBSCRIPT "rm $FilePathBasis_1\.q30.sam\n";
		
		print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar SortSam I=$FilePathBasis_1\.q30.bam O=$FilePathBasis_1\.q30.sorted.bam SO=coordinate 2>> $FilePathBasis_1\.q30.log\n";
		
		print JOBSCRIPT "rm $FilePathBasis_1\.q30.bam\n";
		
		if ($FastqDir =~ m/research\/wst\/abcc6\/ford/){
			
			 print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar MarkDuplicates I=$FilePathBasis_1\.q30.sorted.bam O=$FilePathBasis_1\.q30.sorted.remdup.bam M=$FilePathBasis_1\.q30.sorted.remdup.met REMOVE_DUPLICATES=true 2>> $FilePathBasis_1\.q30.log\n";
		}
		else {

			 print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar MarkDuplicates I=$FilePathBasis_1\.q30.sorted.bam O=$FilePathBasis_1\.q30.sorted.remdup.bam M=$FilePathBasis_1\.q30.sorted.remdup.met REMOVE_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 2>> $FilePathBasis_1\.q30.log\n";
		}
		
		print JOBSCRIPT "rm $FilePathBasis_1\.q30.sorted.bam\n";
		
		print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar AddOrReplaceReadGroups I=$FilePathBasis_1\.q30.sorted.remdup.bam O=$FilePathBasis_1\.q30.sorted.remdup.rg.bam RGID=$RG RGLB=lib1 RGPL=illumina RGPU=$Sample RGSM=$Sample\n";
		
		# nomasking #
		
		print JOBSCRIPT "bwa mem -t 8 -M /user/scratch/gent/gvo000/gvo00082/genomes/GRCh37/homo_sapiens_GRCh37.fa $FirstTrimmedPairedFilePath $SecondTrimmedPairedFilePath > $FilePathBasis_2\.q30.sam 2>> $FilePathBasis_2\.q30.log\n";
		print JOBSCRIPT "samtools view -Sb $FilePathBasis_2\.q30.sam > $FilePathBasis_2\.q30.bam 2>> $FilePathBasis_2\.q30.log\n";
		print JOBSCRIPT "rm $FilePathBasis_2\.q30.sam\n";
		print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar SortSam I=$FilePathBasis_2\.q30.bam O=$FilePathBasis_2\.q30.sorted.bam SO=coordinate 2>> $FilePathBasis_2\.q30.log\n";
		
		print JOBSCRIPT "rm $FilePathBasis_2\.q30.bam\n";

		if ($FastqDir =~ m/research\/wst\/abcc6\/ford/){

				 print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar MarkDuplicates I=$FilePathBasis_2\.q30.sorted.bam O=$FilePathBasis_2\.q30.sorted.remdup.bam M=$FilePathBasis_2\.q30.sorted.remdup.met REMOVE_DUPLICATES=true 2>> $FilePathBasis_2\.q30.log\n";
		}
		else {

				 print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar MarkDuplicates I=$FilePathBasis_2\.q30.sorted.bam O=$FilePathBasis_2\.q30.sorted.remdup.bam M=$FilePathBasis_2\.q30.sorted.remdup.met REMOVE_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 2>> $FilePathBasis_2\.q30.log\n";
		}
		
		print JOBSCRIPT "rm $FilePathBasis_2\.q30.sorted.bam\n";
		
		print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar AddOrReplaceReadGroups I=$FilePathBasis_2\.q30.sorted.remdup.bam O=$FilePathBasis_2\.q30.sorted.remdup.rg.bam RGID=$RG RGLB=lib1 RGPL=illumina RGPU=$Sample RGSM=$Sample\n";
				
		close JOBSCRIPT;
		
		system ("qsub $ShellScript");
	}
}
closedir (FASTQDIR);
