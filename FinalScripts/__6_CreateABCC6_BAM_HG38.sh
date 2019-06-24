#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
ml load SAMtools/1.8-intel-2018a
cd /user/data/gent/gvo000/gvo00082/vsc41234/bcbio/exome/NHHW160062b/
find . -name *bam | while read FilePath; do samtools view -b -L /kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/bed/ABCC6_ABCC6P1_ABCC6P2_Genes_HG38.bed $FilePath > /kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/population/Hg38_Bams/${FilePath##*/} ;done