#######################################################################################
###						      	Calculate ABCC6 BAMs HG38			    			###
#######################################################################################

use strict;
use warnings;
use Getopt::Long;

#######################################################################################
###						      	   Initialize & declare		    					###
#######################################################################################

my $BCBIO_MAIN_DIR 		= "/user/data/gent/gvo000/gvo00082/vsc41234/bcbio/exome/";
my $BedFileGenesHg38	= "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/bed/ABCC6_ABCC6P1_ABCC6P2_Genes_HG38.bed";
my $Hg38_Bams_Dir		= "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/population/Hg38_Bams/";

#######################################################################################
###						      	   Initialize & declare		    					###
#######################################################################################

opendir(BCBIO_MAIN, $BCBIO_MAIN_DIR) or die ("Can't open $BCBIO_MAIN_DIR\n");
while (my $RunDir = readdir(BCBIO_MAIN)){

	my $BCBIO_RUN_DIR = $BCBIO_MAIN_DIR . $RunDir . "/";
	
	print $BCBIO_RUN_DIR . "\n";
	
	next if (! -d $BCBIO_RUN_DIR);
	
	opendir(BCBIO_RUN, $BCBIO_RUN_DIR) or die ("Can't open $BCBIO_RUN_DIR\n");
	while (my $PatientDir = readdir(BCBIO_RUN)){
	
		my $PatientDirFP = $BCBIO_RUN_DIR . $PatientDir . "/final/" . $PatientDir . "/";
	
		if ($PatientDir =~ m/^D/ && -d $PatientDirFP){
		
			print "\t" . $PatientDir . "\n";
		
			my $BamFilePath 	= $PatientDirFP . "$PatientDir\-ready.bam";
			
			my $Command 		= "samtools view -b -L $BedFileGenesHg38 $BamFilePath > $Hg38_Bams_Dir$PatientDir\.abcc6.bam";
			my $JobScript		= "$Hg38_Bams_Dir$PatientDir\.abcc6.sh";
			
			open JOB, ">$JobScript" or die ("Can't open Proband_$JobScript\n");
				
			print JOB "#PBS -l nodes=1:ppn=1\n";
			print JOB "#PBS -l walltime=01:00:00\n";
			print JOB "#PBS -d $Hg38_Bams_Dir\n";
			
			print JOB "ml load SAMtools/1.8-intel-2018a\n";
			
			print JOB "samtools view -b -L $BedFileGenesHg38 $BamFilePath > $Hg38_Bams_Dir$PatientDir\.abcc6.bam\n";

			close JOB;
			
			system ("qsub $JobScript");

		}
	}
	closedir(BCBIO_RUN);
}
closedir(BCBIO_MAIN);