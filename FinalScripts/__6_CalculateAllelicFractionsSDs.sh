#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -l mem=16GB

perl /kyukon/scratch/gent/gvo000/gvo00082/research/wst/abcc6/scripts/__6_CalculateAllelicFractionsSDs.pl --AnalyzeVCF=yes