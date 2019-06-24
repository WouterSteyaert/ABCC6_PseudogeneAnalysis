use strict;

my $MainDir = "/kyukon/scratch/gent/vo/000/gvo00082/research/wst/abcc6/AllAlignments/";

my @Enr 	= ("GSEA", "WES");
my @Anas	= ("masking", "nomasking");
my %DNAs	= ();

foreach my $Enr (@Enr){
	foreach my $Ana (@Anas){
	
		my $Dir = $MainDir . $Enr . "/" . $Ana . "/";
		
		opendir(DIR, $Dir) or die ("Can't open $Dir\n");
		while (my $File = readdir(DIR)){
			
			next if ($File =~ m/^\./);
			
			my 	$DNA = $File;
				$DNA =~ s/\..*//g;
				
			$DNAs{$DNA} = undef;
		}
		
		closedir(DIR);
	}
}

my $Count = 1;

foreach my $DNA (sort keys %DNAs){

	$DNAs{$DNA} = "DNA$Count";
	
	$Count++;
}

foreach my $Enr (@Enr){
	foreach my $Ana (@Anas){
	
		my $Dir = $MainDir . $Enr . "/" . $Ana . "/";
		
		opendir(DIR, $Dir) or die ("Can't open $Dir\n");
		while (my $File = readdir(DIR)){
			
			next if ($File =~ m/^\./);
			
			my 	$DNA = $File;
				$DNA =~ s/\..*//g;
			
			my 	$DNAQM 	= quotemeta($DNA);
			my 	$NRQM	= quotemeta($DNAs{$DNA});
				
			my 	$NewFileName = $File;
				$NewFileName =~ s/$DNAQM/$NRQM/g;
			
			print $File . "\t" . $NewFileName . "\n";
			
			my $FilePath 	= $Dir . $File;
			my $NewFilePath = $Dir . $NewFileName;
			
			rename $FilePath, $NewFilePath;
				
			
		}
		
		closedir(DIR);
	}
}