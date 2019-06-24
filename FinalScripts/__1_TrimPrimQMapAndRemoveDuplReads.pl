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

my $FastqDir   	 		= "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/fordValidationSeq2/fastqs/";
my $TrimmedFastqDir		= "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/fordValidationSeq2/trimmedfastqs/";
my $AnalysisDir_1		= "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/fordValidationSeq2/masking/analysis/";
my $AnalysisDir_2       = "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/fordValidationSeq2/nomasking/analysis/";

my $FivePrimeTrim		= "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/fordprimers/5Prime.fasta";
my $ThreePrimeTrim		= "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/fordprimers/3Prime.fasta";

system ("rm  $TrimmedFastqDir\/*");
system ("rm  $AnalysisDir_1\/*");
system ("rm  $AnalysisDir_2\/*");

#######################################################################################
### 	  	 				    	  		LOOP 									###
#######################################################################################

opendir (FASTQDIR, $FastqDir) or die ("Can't open $FastqDir!\n"); 
while (my $FirstFastQGZFile = readdir(FASTQDIR)){
	
	if ($FirstFastQGZFile =~ m/.*_R1_001\.fastq\.gz/){
	# if ($FirstFastQGZFile =~ m/.*_1\.fastq\.gz/){
		
		 my $FirstFilePath						= $FastqDir . $FirstFastQGZFile;
		 my $SecondFastQFile 					= $FirstFastQGZFile;
		    # $SecondFastQFile 					=~ s/_1/_2/g;
		    $SecondFastQFile 					=~ s/_R1_/_R2_/g;
		 my $SecondFilePath 					= $FastqDir . $SecondFastQFile;
		# (my $FileBasis) 						= $FirstFastQGZFile =~ m/(.*)_1\.fastq\.gz$/;
		(my $FileBasis) 						= $FirstFastQGZFile =~ m/(.+?)_.*\.fastq\.gz$/;
		 my $FilePathBasis_1					= $AnalysisDir_1 . $FileBasis;
		 my $FilePathBasis_2					= $AnalysisDir_2 . $FileBasis;
		 
		 my $ShellScript                		= $TrimmedFastqDir . $FileBasis . ".sh";
		 
		print $FileBasis . "\n";
		
		# $FileBasis =~ m/(.+?)_.*\.fastq\.gz$/;
		
		my 	$Sample 	= $FileBasis;
		my 	$RG 		= 1;
			# $RG			=~ s/^0+//g;
	
		open JOBSCRIPT, ">$ShellScript" or die ("Can't open ShellScript\n");
		
		print JOBSCRIPT "#PBS -N $FileBasis\n";
		print JOBSCRIPT	"#PBS -l nodes=1:ppn=8\n";
		print JOBSCRIPT	"#PBS -l walltime=09:00:00\n";
		print JOBSCRIPT	"#PBS -d $TrimmedFastqDir\n";
		print JOBSCRIPT "cd $TrimmedFastqDir\n";
		
		print JOBSCRIPT "module load Trimmomatic/0.32-Java-1.7.0_40\n";
		print JOBSCRIPT "module load cutadapt/1.8.1-intel-2015a-Python-2.7.9\n";
		
		print JOBSCRIPT "cutadapt -g file:$FivePrimeTrim -a file:$ThreePrimeTrim -o $FileBasis\_1\.prim.fastq.gz -p $FileBasis\_2\.prim.fastq.gz $FirstFilePath $SecondFilePath --minimum-length 10 -e 0.0001\n";
		
		print JOBSCRIPT "java -jar \${EBROOTTRIMMOMATIC}\/trimmomatic-0.32.jar PE $FileBasis\_1\.prim.fastq.gz $FileBasis\_2\.prim.fastq.gz $FileBasis\_1\.prim.q30.paired.fastq.gz $FileBasis\_1\.prim.q30.unpaired.fastq.gz $FileBasis\_2\.prim.q30.paired.fastq.gz $FileBasis\_2\.prim.q30.unpaired.fastq.gz LEADING:30 TRAILING:30\n";
		
		# print JOBSCRIPT "rm $FirstFilePath\n";
		# print JOBSCRIPT "rm $SecondFilePath\n";
		
		print JOBSCRIPT "module purge\n";
		print JOBSCRIPT	"module load BWA/0.7.13-intel-2016a\n";
		print JOBSCRIPT	"module load SAMtools/1.3-intel-2016a\n";
		print JOBSCRIPT	"module load Java/1.8.0_74\n";
		print JOBSCRIPT	"module load picard/2.1.1-Java-1.8.0_74\n";
		
		# masking #
		
		print JOBSCRIPT "bwa mem -t 8 -M /user/scratch/gent/gvo000/gvo00082/genomes/GRCh37_Masked_ABCC6P1P2/homo_sapiens_GRCh37_MaskedForABCC6P1P2.fa $FileBasis\_1\.prim.q30.paired.fastq.gz $FileBasis\_2\.prim.q30.paired.fastq.gz > $FilePathBasis_1\.q30.sam 2>> $FilePathBasis_1\.q30.log\n";
		
		print JOBSCRIPT "samtools view -Sb $FilePathBasis_1\.q30.sam > $FilePathBasis_1\.q30.bam 2>> $FilePathBasis_1\.q30.log\n";
		
		print JOBSCRIPT "rm $FilePathBasis_1\.q30.sam\n";
		
		print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar SortSam I=$FilePathBasis_1\.q30.bam O=$FilePathBasis_1\.q30.sorted.bam SO=coordinate 2>> $FilePathBasis_1\.q30.log\n";
		
		if ($FastqDir =~ m/research\/wst\/abcc6\/ford/){
			
			 print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar MarkDuplicates I=$FilePathBasis_1\.q30.sorted.bam O=$FilePathBasis_1\.q30.sorted.remdup.bam M=$FilePathBasis_1\.q30.sorted.remdup.met REMOVE_DUPLICATES=true 2>> $FilePathBasis_1\.q30.log\n";
		}
		else {

			 print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar MarkDuplicates I=$FilePathBasis_1\.q30.sorted.bam O=$FilePathBasis_1\.q30.sorted.remdup.bam M=$FilePathBasis_1\.q30.sorted.remdup.met REMOVE_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 2>> $FilePathBasis_1\.q30.log\n";
		}
		
		print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar AddOrReplaceReadGroups I=$FilePathBasis_1\.q30.sorted.remdup.bam O=$FilePathBasis_1\.q30.sorted.remdup.rg.bam RGID=$RG RGLB=lib1 RGPL=illumina RGPU=$Sample RGSM=$Sample\n";
		
		# nomasking #
		
		print JOBSCRIPT "bwa mem -t 8 -M /user/scratch/gent/gvo000/gvo00082/genomes/GRCh37/homo_sapiens_GRCh37.fa $FileBasis\_1\.prim.q30.paired.fastq.gz $FileBasis\_2\.prim.q30.paired.fastq.gz > $FilePathBasis_2\.q30.sam 2>> $FilePathBasis_2\.q30.log\n";
		
		print JOBSCRIPT "samtools view -Sb $FilePathBasis_2\.q30.sam > $FilePathBasis_2\.q30.bam 2>> $FilePathBasis_2\.q30.log\n";
		
		print JOBSCRIPT "rm $FilePathBasis_2\.q30.sam\n";
		
		print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar SortSam I=$FilePathBasis_2\.q30.bam O=$FilePathBasis_2\.q30.sorted.bam SO=coordinate 2>> $FilePathBasis_2\.q30.log\n";
		
		if ($FastqDir =~ m/research\/wst\/abcc6\/ford/){
			
			 print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar MarkDuplicates I=$FilePathBasis_2\.q30.sorted.bam O=$FilePathBasis_2\.q30.sorted.remdup.bam M=$FilePathBasis_2\.q30.sorted.remdup.met REMOVE_DUPLICATES=true 2>> $FilePathBasis_2\.q30.log\n";
		}
		else {

			 print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar MarkDuplicates I=$FilePathBasis_2\.q30.sorted.bam O=$FilePathBasis_2\.q30.sorted.remdup.bam M=$FilePathBasis_2\.q30.sorted.remdup.met REMOVE_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 2>> $FilePathBasis_2\.q30.log\n";
		}
		
		print JOBSCRIPT "java -jar \${EBROOTPICARD}\/picard.jar AddOrReplaceReadGroups I=$FilePathBasis_2\.q30.sorted.remdup.bam O=$FilePathBasis_2\.q30.sorted.remdup.rg.bam RGID=$RG RGLB=lib1 RGPL=illumina RGPU=$Sample RGSM=$Sample\n";
		
		close JOBSCRIPT;
		
		# system ("qsub $ShellScript");
	}
}
closedir (FASTQDIR);
