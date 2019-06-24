# ml load SAMtools/1.8-intel-2018a
# ml load picard/2.18.5-Java-1.8.0_162
use strict;

my $PoolsToDNA_FilePath = "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/SequencingShanaFordPools.txt";
my %PoolsToDNAs			= ();

open P, "$PoolsToDNA_FilePath" or die ("Can't open $PoolsToDNA_FilePath\n");
while (<P>){

	my 	$Line = $_;
		$Line =~ s/\n//g;
	(my $DNA, my $Pool) = split (m/\t/, $Line);
	
	$PoolsToDNAs{$Pool} = $DNA;
}
close P;

my @Dirs = ("/user/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/nomasking/analysis/",
			"/user/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/masking/analysis/");

foreach my $Dir (@Dirs){

	my %DNA_ToBam = ();

	opendir(DIR, $Dir) or die $!;
	while (my $BamFile = readdir(DIR)){

		next if ($BamFile !~ m/bam$/);

		my $OldBamFilePath 	= $Dir . $BamFile;
		
		$BamFile =~ m/(.*)\.q30\.sorted\.remdup\.rg\.bam$/;
		
		if (exists $PoolsToDNAs{$1}){
		
			if ($DNA_ToBam{$PoolsToDNAs{$1}}){
				$DNA_ToBam{$PoolsToDNAs{$1}} .= " $OldBamFilePath";
			}
			else {
			
				$DNA_ToBam{$PoolsToDNAs{$1}} = "$OldBamFilePath";
			}	
		}
		else {
		
			system ("rm $OldBamFilePath");
		}
	}
	
	foreach my $DNA (keys %DNA_ToBam){
	
		print "$Dir$DNA\.q30.sorted.remdup.rg.bam $DNA_ToBam{$DNA}\n\n";

		system ("samtools merge -f $Dir$DNA\.q30.sorted.remdup.bam $DNA_ToBam{$DNA}");

		system ("java -jar \$EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=$Dir$DNA\.q30.sorted.remdup.bam O=$Dir$DNA\.q30.sorted.remdup.rg.bam RGID=$DNA RGLB=1 RGPL=ILLUMINA RGPU=$DNA RGSM=$DNA");
	}
}


