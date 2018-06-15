>  Selene Fernandez ([Syllabus](https://liz-fernandez.github.io/Escuela_Bioinfo_2018_Non_Model/))

### Learning aims:

* What are model and non-model organims
* why and how study model organisms
* High-thorughput sequencing techniques
* ...

**Model organisms**

Are non-human species tht are used in the lab (easy to culture) to help scientists understand biological processes.

**Non model organisms**

Life is enormously diverse!!. Many biological phenomena of interest can be very specific characterustucs of certain organims that can not be easly maintained in the laboratory. The revolution in genomic techniques allows us to acquire genetic information of all types of organisms, even if they can only be collected in the field and cannot be grown in a laboratory settings.

**High-thorughput sequencing techniques**

- Pyrosequencing (Roche 454, long sequencing 500 bp, )
- _Sequencing by synthesis_ (Ilumina - Process start by joining adapters to the DNA or RNA fragments that we want to sequence)
- Sequencing by ligation (pacBio, is expensive, 860-1199 bp, )
- Ion semiconductor 
- Nanopore (wowowow)

### Source of error in the field:

There are two main source of error:

- Human: Mixing of samples (in lab or when files were received) , errors in the protocol.
- Technical error: Errors inherent ot the platform (e.g mononucleotide sequence iin pyrosequencing) - all platforms have some level of error than must been taken in account:
  - Sample prep error
    - User error (e.g mis-labeling )

## Practice: Assembly genome

- Learn how to assemble a genome of a non-model organism *de novo*.
- Learn the descriptive genome statistics to assess an assembly.

In this session, we will be working with the non model organism *Leucothrix mucor*, which is a marine bacteria with a epiphyte lifestyle and has an estimated genome size of ~5 Mb. The sequencing was performed using Illumina’s HiSeq 2000 with a paired-end configuration of 2x150.

> Starting in lavish:

```
ssh ricardo@kneipe.lavis.unam.mx
ricardo@kneipe.lavis.unam.mx's password: ***
```



> Copy directory with dataset into your home:

```
 cp -r /home/sfernandez/DATA/GENOMICS .
 cd GENOMICS
```



```bash
ln -s /home/sfernandez/TOOLS/Trimmomatic-0.38/trimmomatic-0.38.jar .

java -jar trimmomatic-0.38.jar PE Reads_1.fq.gz Reads_2.fq.gz  Reads_1.trimm.P.fq.gz Reads_1.trimm.UP.fq.gz Reads_2.trimm.P.fq.gz Reads_2.trimm.UP.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:24 MINLEN:60
```

As result:

Input Read Pairs: 2000000 Both Surviving: 1998803 (99.94%) Forward Only Surviving: 1191 (0.06%) Reverse Only Surviving: 0 (0.00%) Dropped: 6 (0.00%) 

**Note**: Both, unpaired (UP) and Paired (P) files contains the filtered reads than trimmomatic trims, but **UP contain reads than only one strand pass the filter.**

> lets cat paired files

```bash
cat *.UP.* > Reads.UP.fq.gz
cat *.P.* > Reads.P.fq.g
```

> Then, use SPAdes

```bash
module add spades/3.12.0
```

Code the spades

```bash
spades.py --careful --pe1-1 Reads_1.trimm.P.fq.gz --pe1-2 Reads_2.trimm.P.fq.gz --pe1-s Reads.UP.fq.gz -t 1 -m 8 -k 21,33,55 -o  /media/public/ricardo/
```

The parameters we have used do the following:

- `--careful` : This option tells SPAdes to try to eliminate mismatches, insertion and deletions that might be due to sequencing errors. It is important to mention that this option should not be used with big genomes.
- `--pe#-1` / `--pe#-2`: With this parameter we will pass onto SPAdes each of our files containing our sequencing data, being # the number of sample. Unfortunately, SPAdes only accepts up to 9 samples, however, this can be worked around by merging all the files into one.
- `--pe#-s`: With this parameter we can make use of our unpaired reads to further improve our assembly. For example, paired reads can become unpaired during the error correction procedure.
- `-t`: This parameters establishes the number of threads SPAdes is allowed to use.
- `-m`: This option will set the parameter of how much memory will be used in the assembly (always in Gb).
- `-k`: This might be one of the most important parameters of the software, where we can specify multiple k-mer sizes.
- `-o`: Here we will find the directory where we want SPAdes to deposit the output.

> Good!

Let's figure out the number of scaffolds and contains:

```
 cat scaffolds.fasta | grep '>' -c 
```

> stats
>
> ```
> quast.py SPAdes/scaffolds.fasta
> quast.py SPAdes/contigs.fasta
> cat quast_results/results_2018_06_13_12_35_59/report.txt 
> ```
>
> | Assembly       | Contigs | Scaffolds |
> | -------------- | ------- | --------- |
> | Largest contig | 26558   | 26558     |
> | Total length   | 5141357 | 5144321   |
> | GC (%)         | 47.75   | 47.75     |
> | N50            | 4347    | 310       |
> | L50            | 317     | 743       |
>
> 

### Genome _ab initio_ annotation

https://liz-fernandez.github.io/Escuela_Bioinfo_2018_Non_Model/16-ab_initio_Annotation.html



# De novo transcriptome Assembly

- Learn how the *de novo* RNA-Seq assembly works.
- Learn how to use Trinity to assemble a transcriptome *de novo*.

It's to take RNA-seq libraries and assembled them into complete transcriptss **without the guide of a reference genome**

> Let's download data

```bash
path=liz-fernandez.github.io/Escuela_Bioinfo_2018_Non_Model/datasets

wget https://$path/Sp_ds.left.fq.gz 
wget https://$path/Sp_ds.right.fq.gz
wget https://$path/Sp_hs.right.fq.gz 
wget https://$path/Sp_hs.left.fq.gz   
wget https://$path/Sp_plat.left.fq.gz
wget https://$path/Sp_plat.right.fq.gz
wget https://$path/Sp_log.left.fq.gz  
wget https://$path/Sp_log.right.fq.gz

```

> Running trinity

```bash
$ module load bowtie2
$ module load jellyfish

/home/sfernandez/TOOLS/trinityrnaseq-Trinity-v2.5.0/Trinity --seqType fq --SS_lib_type RF \
--left Sp_log.left.fq.gz,Sp_hs.left.fq.gz \
--right Sp_log.right.fq.gz,Sp_hs.right.fq.gz \
--CPU 2 --max_memory 1G --no_normalize_reads
```

This command will take approximately 15 minutes to assemble a transcriptome.

The options (flags) that we have used are the following:

- –seqType fq - We indicate that we are using fastq type files
- –SS_lib_type RF - We indicate that the library was built using paired-end reads in the RF orientation (reverse-forward)
- –left - Left side reads (o R1)
- –right - Right side reads (o R2)
- –CPU 2 - Use 2 CPUs
- –max_memory 1G - Indicate that the maximum RAM to use is 1GB
- –no_normalize_reads - Do not carry out digital normalization of reads

These are the most essential options to carry out the analysis. Therefore it is very important to know the orientation of the library which depends on the experimental protocol that was used to generate it (ex. dUTP).



## Align good reads vs reference 

We will use the fastq files that we used in the previous practice, as well as the reference genome of our organism, *Saccharomyces pombe*:

- Align sequencing data to a reference (genomes or transcriptomes)
- Understand how to interpret mass sequencing alignment data
- First approach to the SAM and BAM coordinate formats

```bash
wget https://liz-fernandez.github.io/PBI_transcriptomics/datasets/genome/Sp_genome.fa ./reference
```

> Within the WorkDirectory use tophat (first index the reference):

```bash
bowtie2-build ./reference/Sp_genome.fa ./reference/Sp_genome 

tophat2 -I 300 -i 20 reference/Sp_genome \
Sp_log.left.fq.gz,Sp_hs.left.fq.gz,Sp_ds.left.fq.gz,Sp_plat.left.fq.gz \
Sp_log.right.fq.gz,Sp_hs.right.fq.gz,Sp_ds.right.fq.gz,Sp_plat.right.fq.gz -o tophat_reference
```

**Note**: Option `-i` means min intron size whereas `-I` max Intron size. (It values seem well for _S.pompe_ **but values depend upon your model**)

### Exploring SAM results 

```bash
module load samtools
samtools view tophat_out/accepted_hits.bam | head
```

## Align good reads vs transcriptome assembly

We have aligned the readings to the genome but we also want to align them directly to the transcriptome. We are going to review the TopHat2 manual and use the options that allows us to map readings directly to transcriptomes.

> indexing assembly first

```bash
 bowtie2-build ./trinity_out_dir/Trinity.fasta ./trinity_out_dir/Trinity
```

> Run aligment with tophat

```bash
tophat2 -I 300 -i 20 ./trinity_out_dir/Trinity Sp_ds.left.fq.gz,Sp_hs.left.fq.gz,Sp_log.left.fq.gz,Sp_plat.left.fq.gz Sp_ds.right.fq.gz,Sp_hs.right.fq.gz,Sp_log.right.fq.gz,Sp_plat.right.fq.gz -o tophat_assembly -o tophat_assembly
```

> Good!

Let's compare both alignments, the transcriptomic and reference one.

Sort first

```
samtools sort tophat_assembly/accepted_hits.bam > assembly.bam
# 2
samtools sort tophat_reference/accepted_hits.bam > reference.bam
```

Then index:

```bash
samtools index assembly.bam
samtools index reference.bam
```

> Stats!!

```bash
module load bamtools/2.5.1 
bamtools stats -in reference.bam
bamtools stats -in assembly.bam
```

Results in!

|                   | Assembly | reference |
| ----------------- | -------- | --------- |
| Total reads       | 577.499  | 646.837   |
| Mapped reads      | 577.499  | 646.837   |
| Forward strand    | 292.524  | 320.824   |
| Reverse strand    | 284.975  | 326.013   |
| Failed QC         | 0        | 0         |
| Duplicates        | 0        | 0         |
| Paired-end reads  | 577.499  | 646.837   |
| Proper-pairs'     | 69.128   | 80.432    |
| Both pairs mapped | 529.994  | 619.436   |
| Read 1            | 287.111  | 319.135   |
| Read 2            | 290.388  | 327.702   |
| Singletons        | 47.505   | 27.401    |

# Guide transcriptome assembly

- Learn how the guided assembly works.
- Learn to use Cufflinks to assemble de novo data.

```
module load cufflinks/2.2.1
```

Index one more time:

```
mkdir CUFFLINKS && cd CUFFLINKS
# 2
ln -s ../Sp_genome.fa .
bowtie2-build Sp_genome.fa Sp_genome
```

We run separately each library in the guide transcriptome assembly

```
tophat2 -I 1000 -i 20 --library-type fr-firststrand -o tophat.Sp_ds.dir Sp_genome ../Sp_ds.left.fq.gz ../Sp_ds.right.fq.gz
```

Then index and run cufflinks

```
samtools index tophat.Sp_ds.dir/accepted_hits.bam
# 2

/home/sfernandez/TOOLS/cufflinks-2.2.1.Linux_x86_64/cufflinks --overlap-radius 1 \
             --library-type fr-firststrand \
             -o cufflinks.Sp_ds.dir tophat.Sp_ds.dir/accepted_hits.bam
```



> Try a loop for the rest set of libs :)

```bash
for i in hs log plat
do
tophat2 -I 1000 -i 20 --library-type fr-firststrand -o tophat.$i Sp_genome ../Sp_$i.left.fq.gz ../Sp_$i.right.fq.gz
done
```

2nd loop step

```bash
for i in hs log plat
do
mv tophat.$i/accepted_hits.bam tophat.$i/${i}_hits.bam;
samtools index tophat.$i/${i}_hits.bam;
done
```

then,

```bash
for i in hs log plat
do
/home/sfernandez/TOOLS/cufflinks-2.2.1.Linux_x86_64/cufflinks --overlap-radius 1 \
             --library-type fr-firststrand \
             -o cufflinks.${i}.dir tophat.${i}/${i}_hits.bam
             
done
```

Finally rename the cufflink output file (gtf)

```bash
for i in hs log plat ds
do
mv cufflinks.${i}.dir/transcripts.gtf cufflinks.${i}.dir/Sp_${i}.transcripts.gtf
done
```

And then, make a final input matrix

```bash
ls cufflinks.*/*transcripts.gtf > assemblies.txt
# 2
cuffmerge -s Sp_genome.fa assemblies.txt
```

> good!

or

```bash
for i in hs log plat ds
do
echo "tophat2 -I 1000 -i 20 --library-type fr-firststrand -o tophat.Sp.${i}.dir Sp_genome ../Sp.${i}.left.fq.gz ../Sp.${i}.right.fq.gz; mv tophat.Sp.${i}.dir/accepted_hits.bam tophat.Sp.${i}.dir/Sp.${i}.bam; samtools index tophat.Sp.${i}.dir/Sp.${i}.bam; /home/sfernandez/TOOLS/cufflinks-2.2.1.Linux_x86_64/cufflinks --overlap-radius 1 \
             --library-type fr-firststrand \
             -o cufflinks.Sp.${i}.dir tophat.Sp.${i}.dir/Sp.${i}.bam ; mv cufflinks.Sp.${i}.dir/transcripts.gtf cufflinks.Sp.${i}.dir/Sp.${i}.transcripts.gtf"
done | sh

ls cufflinks.Sp*/*transcripts.gtf > assemblies.txt
# 2
cuffmerge -s Sp_genome.fa assemblies.txt
```

### Visualization

in IGV:

Open > File > Load FILE > Select within CUFFLINKS directory to load

1. Referece than use to align reads (genome/transcriptome assembly)
2. Merged *gtf file
3. Merge independent each *gtf file from condition (ex.cufflinks.ds.dir/Sp_ds.transcripts.gtf)
4. hits bam files from tophat (ex. Tophat.hs/hs_hits.bam) 

Use your reference