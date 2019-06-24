use strict;

my $PoolsToDNA_FilePath = "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/SequencingShanaFordPools.txt";
my %PoolsToDNAs			= ();

open P, "$PoolsToDNA_FilePath" or die ("Can't open $PoolsToDNA_FilePath\n");
while (<P>){

	my 	$Line = $_;
		$Line =~ s/\n//g;
	(my $DNA, my $Pool) = split (m/\t/, $Line);
	
	$PoolsToDNAs{$Pool} = $DNA;
	
	system ("cp /kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/fastqsall/$Pool\* /kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/fastqs/");
}
close P;

opendir(DIR, "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/fastqs/") or die $!;
while (my $SeqFile = readdir(DIR)){

	next if ($SeqFile !~ m/gz$/);

	my @SeqFileValues = split(m/_/, $SeqFile);
	
	my $DNA 	= $PoolsToDNAs{$SeqFileValues[0]};
	my $SeqDir	= "/";
	
	if 		($SeqFile =~ m/_R1_/){$SeqDir = 1;}
	elsif 	($SeqFile =~ m/_R2_/){$SeqDir = 2;}
	else	 					 {die ($!);}
		
	my $OldFilePath 	= "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/fastqs/" . $SeqFile;
	my $NewFilePath		= "/kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/fordBWsamples/fastqs/" . $DNA . "_$SeqDir\.fastq.gz";
	
	if (-f $NewFilePath){
		
		system ("cat $OldFilePath >> $NewFilePath");
		system ("rm $OldFilePath");
	}
	else {
		
		system ("mv $OldFilePath $NewFilePath");
	}
}
closedir (DIR);


