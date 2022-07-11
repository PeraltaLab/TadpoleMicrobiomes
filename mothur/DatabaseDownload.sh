#!/bin/bash

wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v132.tgz
tar -xzvf silva.nr_v132.tgz

wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138_1.tgz
tar -xzvf silva.nr_v138_1.tgz

wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.gold.bacteria.zip
unzip -j silva.gold.bacteria.zip

wget https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset18_062020.rdp.tgz
tar -xzvf  trainset18_062020.rdp.tgz

source /center1/MICROBES/Mothur/mothur-v.1.47.0/setEnvironment.sh

# Test pcr seqs with small files: 
# head -200 BNZ16S.trim.contigs.good.unique.fasta > align.test.fasta
# align.seqs(fasta=align.test.fasta, reference=silva.nr_v132.align)
# summary.seqs(fasta=align.test.align)
mothur "#pcr.seqs(fasta=silva.nr_v132.align, start=12000, end=30000, keepdots=F, processors=1)"
mv silva.nr_v132.pcr.align silva_v132.v4v5.fasta

# align.seqs(fasta=align.test.fasta, reference=silva.nr_v138_1.align)
# summary.seqs(fasta=align.test.align)
mothur "#pcr.seqs(fasta=silva.nr_v138_1.align, start=12000, end=30000, keepdots=F, processors=1)"
mv silva.nr_v138_1.pcr.align silva_v138_1.v4v5.fasta
