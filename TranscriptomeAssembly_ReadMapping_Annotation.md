## Transcriptome Assembly

This workbook details transcriptome assembly and mapping.

First, assembly of the _Ranitomeya imitator_ transcriptome: 

```
#!/bin/bash
#SBATCH -N 1
#SBATCH -t 140:00:00
# echo commands to stdout
set -x

# first, concatenate reads by imitator treatment (egg fed, lab fed, not fed)
for i in reads/KW_Ri_fi*_1.fq.gz; do zcat $i; done > imitator_eggs_reads_R1.fastq &
for f in reads/KW_Ri_fi*_2.fq.gz; do zcat $f; done > imitator_eggs_reads_R2.fastq &


for g in reads/KW_Ri_lf*_1.fq.gz; do zcat $g; done > imitator_lab_reads_R1.fastq &
for h in reads/KW_Ri_lf*_2.fq.gz; do zcat $h; done > imitator_lab_reads_R2.fastq &


for l in reads/KW_Ri_nf*_1.fq.gz; do zcat $l; done > imitator_notfed_reads_R1.fastq &
for m in reads/KW_Ri_nf*_2.fq.gz; do zcat $m; done > imitator_notfed_reads_R2.fastq


# subsample each treatment to equal coverage
seqtk sample -s100 imitator_eggs_reads_R1.fastq 13400000 > subsamp.imitator.eggs.R1.fastq &
seqtk sample -s100 imitator_eggs_reads_R2.fastq 13400000 > subsamp.imitator.eggs.R2.fastq &

seqtk sample -s100 imitator_lab_reads_R1.fastq 13400000 > subsamp.imitator.lab.R1.fastq &
seqtk sample -s100 imitator_lab_reads_R2.fastq 13400000 > subsamp.imitator.lab.R2.fastq &

seqtk sample -s100 imitator_notfed_reads_R1.fastq 13400000 > subsamp.imitator.notfed.R1.fastq &
seqtk sample -s100 imitator_notfed_reads_R2.fastq 13400000 > subsamp.imitator.notfed.R2.fastq



# concatenate forward and reverse reads, each should be a little over 40 million reads

for i in subsamp.imitator*R1.fastq; do cat $i; done > subsamp.imitator.R1.fastq &
for f in subsamp.imitator*R2.fastq; do cat $f; done > subsamp.imitator.R2.fastq



# delete intermediate files
rm imitator_eggs* imitator_lab* imitator_not* subsamp.imitator.eggs* subsamp.imitator.lab* subsamp.imitator.notfed*



source activate orp_v2

/mnt/lustre/macmaneslab/ams1236/Oyster_River_Protocol/oyster.mk main \
MEM=120 \
CPU=24 \
READ1=subsamp.imitator.R1.fastq \
READ2=subsamp.imitator.R2.fastq \
RUNOUT=imi_sub
```


Second, assembly of the _Ranitomeya variabilis_ transcriptome: 

```
#!/bin/bash
#SBATCH -N 1
#SBATCH -t 140:00:00
# echo commands to stdout
set -x

# variabilis tadpole gut transcriptome assembly and downstream bash work

# start in kayladata/variabilis folder

# first, concatenate reads by variabilis treatment (egg fed, lab fed, not fed)
for i in reads/KW_Rv_fi*_1.fq.gz; do zcat $i; done > variabilis_field_reads_R1.fastq &
for f in reads/KW_Rv_fi*_2.fq.gz; do zcat $f; done > variabilis_field_reads_R2.fastq &


for g in reads/KW_Rv_lf*_1.fq.gz; do zcat $g; done > variabilis_lab_reads_R1.fastq &
for h in reads/KW_Rv_lf*_2.fq.gz; do zcat $h; done > variabilis_lab_reads_R2.fastq 




# subsample each treatment to equal coverage
seqtk sample -s100 variabilis_field_reads_R1.fastq 20000000 > subsamp.variabilis.field.R1.fastq &
seqtk sample -s100 variabilis_field_reads_R2.fastq 20000000 > subsamp.variabilis.field.R2.fastq &

seqtk sample -s100 variabilis_lab_reads_R1.fastq 20000000 > subsamp.variabilis.lab.R1.fastq &
seqtk sample -s100 variabilis_lab_reads_R2.fastq 20000000 > subsamp.variabilis.lab.R2.fastq 




# concatenate forward and reverse reads, each should be a little over 40 million reads

for i in subsamp.variabilis*R1.fastq; do cat $i; done > subsamp.variabilis.R1.fastq &
for f in subsamp.variabilis*R2.fastq; do cat $f; done > subsamp.variabilis.R2.fastq



# delete intermediate files
rm variabilis_field* variabilis_lab* subsamp.variabilis.field* subsamp.variabilis.lab* 

source activate orp_v2


# Reads managed. Run the Oyster River Protocol to build an variabilis tadpole gut transcriptome
/mnt/lustre/macmaneslab/ams1236/Oyster_River_Protocol/oyster.mk main \
MEM=120 \
CPU=24 \
READ1=subsamp.variabilis.R1.fastq \
READ2=subsamp.variabilis.R2.fastq \
RUNOUT=var_sub
```

Third, merge the species-specific transcriptome assemblies using `OrthoFuser`:

```
#!/bin/bash

### Orthofuse the two transcriptomes together
# combine readsets
cat imitator/subsamp.imitator.R1.fastq variabilis/subsamp.variabilis.R1.fastq > subsamp.imivar.R1.fastq
cat imitator/subsamp.imitator.R2.fastq variabilis/subsamp.variabilis.R2.fastq > subsamp.imivar.R2.fastq

# copy assemblies over
mkdir assemblies2/
cp imitator/assemblies/*.ORP.fasta
cp variabilis/assemblies/*.ORP.fasta

# rename the assemblies so they work
cd assemblies2/
awk '/^>/{print ">Ranitomeya_imitator_" ++i; next}{print}'  imi_sub.ORP.fasta > tmp
mv tmp imi_sub.ORP.fasta
awk '/^>/{print ">Ranitomeya_variabilis_" ++i; next}{print}'  var_sub.ORP.fasta > tmp
mv tmp var_sub.ORP.fasta 


# run orthofuser:
/home/summerslab/programs/Oyster_River_Protocol/orthofuser.mk all \
READ1=subsamp.imivar.R1.fastq \
READ2=subsamp.imivar.R2.fastq \
CPU=16 \
MEM=550 \
RUNOUT=Rimivartrans \
FASTADIR=/home/summerslab/kayla2.0/assemblies2
```

## RNA mapping

Pseudo-quantify RNAseq readsets to the assembly using `kallisto`:

```
#!/bin/bash

### Pseudo-quantification
# Make directories for each sample:
cd imitator/reads/
imitator=$(ls *1.fq.gz | sed "s/1.fq.gz//g")
cd ../..
#mkdir kallisto_quants
cd kallisto_quants
for i in $imitator; do mkdir $i; done
cd ..

kallisto index -i Rimivartrans.idx Rimivartrans.ORP.fasta

echo $imitator

for i in $imitator
do
kallisto quant -i Rimivartrans.idx -o kallisto_quants/${i} -t 8 -b 100 \
imitator/reads/${i}1.fq.gz imitator/reads/${i}2.fq.gz
done

cd variabilis/reads/
variabilis=$(ls *1.fq.gz | sed "s/1.fq.gz//g")
cd ../..
cd kallisto_quants
for i in $variabilis; do mkdir $i; done
cd ..

echo $variabilis

for i in $variabilis
do
kallisto quant -i Rimivartrans.idx -o kallisto_quants/${i} -t 8 -b 100 \
variabilis/reads/${i}1.fq.gz variabilis/reads/${i}2.fq.gz
done
```

## Transcriptome Annotation

Annotate the assembly with `Diamond`:

```
#### Transcriptome annotation
diamond makedb -in allpep.fa -d allpep.dmnd

diamond blastx -d /home/summerslab/peptidedb/allpep.dmnd -q Rimivartrans.ORP.fasta -o tad_annotation.txt --threads 16

# sort by top hits
sort tad_annotation.txt -k 1,1 -k11,11g | sort -u -k 1,1 --merge > tad_annotation_tophit.txt
```
