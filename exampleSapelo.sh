#PBS -S /bin/bash
#PBS -q batch
#PBS -N RATT_pipeline
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=4gb
#PBS -l walltime=00:15:00


module load mummer/3.23
module load ratt/10052017
module load python/2.7.8
module load samtools/1.3.1
module load emboss/6.5.7

cd /lustre1/myUsrName/workingDir/


sub_fasta="reference.fasta"
sub_gff="reference.gff"
query_fasta="query.fasta"

RATTwithGFF.py $sub_gff $sub_fasta $query_fasta run1 Assembly


























