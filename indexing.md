### Indexing methods

> From genome directory

```
for i in bwa bowtie2 bwa segemeh1
do
mkdir $i && ln -s ../ce10.fa ./$i
done

```

> Then index by tool config in each director created:

```
# SEGEMEHL
module add  segemehl/0.2.0 
segemehl.x -d ce10.fa -x ce10.idx

# BOWTIE2
module add  bowtie2/2.3.4.1 
bowtie2-build ce10.fa ce10

# BWA
module add bwa/0.7.17
bwa index ce10.fa

```

> Let's mapping

```
segemehl.x -d ce10.fa -i ../segemehl/ce10.idx -q ../../clipped/SRR359063_1.trimmed.fastq.gz -p ../../clipped/SRR359063_2.trimmed.fastq.gz -t 4 -S > SRR359063.segemeh1.sam
```

> Next bowtie2

```
bowtie2 -x ../bowtie2/ce10 -1 ../../clipped/SRR359063_1.trimmed.fastq.gz -2 ../../clipped/SRR359063_2.trimmed.fastq.gz > SRR359063.bt2.sam
```

> Finally, bwa

```
for i in ../../clipped/*trimmed.fastq.gz;  
do
echo "bwa aln ../bwa/ce10.fa $i > ${i%.trimmed.fastq.gz}.sai"; done | sh
 
# second step (sam paired end samples mapping)
bwa sampe ce10.fa SRR359063_1.sai SRR359063_2.sai ../../clipped/SRR359063_1.trimmed.fastq.gz ../../clipped/SRR359063_2.trimmed.fastq.gz > SRR359063.bwa.sam

```

## Open SAM/BAM md

