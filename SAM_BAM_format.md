# SAM file format

> The apprach with SAM is visualize aligments results by ordering:
>
> 1. convert sam to a compressed - binary (bam) file
> 2. Sort the bam format
> 3. Indexing
> 4. Input in visualizers like [picard](http://broadinstitute.github.io/picard/) IGB or [IGV](https://software.broadinstitute.org/software/igv/) 

The first step orders the results by their coordinates and the second one creates indexes to speed up the display using a browser.

```bash
module add samtools/1.8 
samtools view -bS SRR359063.segemeh1.sam > SRR359063.segemeh1.bam

# display file type:
file SRR359063.segemeh1.bam
```

> display the head from the binary file and sort also

```bash
samtools view -H SRR359063.segemeh1.bam
samtools sort SRR359063.segemeh1.bam -o SRR359063.segemeh1.sorted
```

> Then, lets index file in order to input file in a **Visualizer**

```bash
samtools index SRR359063.segemeh1.sorted 
```

> lets use column 2 (bitwise flag) in order to filter (-f) bam file

```
samtools view -f 0x10 SRR359063.segemeh1.sorted  | wc -l
```

La opcion -f reconoce reconoce los resultados de las lecturas reverse, mientras que -F las forward

```bash
samtools view -F 0x10 SRR359063.segemeh1.sorted  | wc -l
```

### bamtools stats (perfect end)

```
module load bamtools/2.5.1 
bamtools stats -in SRR359063.segemeh1.sorted 
```



### cigar string (important when analyze SNP's)

this values report the aligment information . ie. which nucleotide of the read is aligned and where insertions and deletions occur.

> Read report [here](https://genome.sph.umich.edu/wiki/SAM#What_is_a_CIGAR.3F) 

### samtools tview

```
samtools tview SRR359063.segemeh1.sorted ../../../genome/chr/ce10.fa
```

## The SAM/BAM format

The SAM format is a plain text format that allows us to save sequencing data in ASCII format delimited by tabs.

It is made up of two core sections:

- The headers
- The alignment

The section of the **header** starts with the character `@` followed by one of the codes of two letters that denote the characteristics of the alignments in this file. Each line is delimited by tabs and, in addition to the lines that begin with `@CO`, each data field has the format`TAG:VALUE`, where `TAG` is a string Two characters that define the format and content of `VALUE`.

The header is not indispensable but contains information about the version of the file as well as if it is ordered or not. Therefore it is advisable to include it.

```
@HD The header line; VN: Format version; SO: Sorting order of alignments

@

@PG

```

The aligment section contains the follow information:

1. QNAME .... https://liz-fernandez.github.io/Escuela_Bioinfo_2018_Non_Model/04-mapping.html

### The alignment part

In the sam file each column display details from the alignment:

1. Query template name
2. bitwise flag:  determine good aligments (details [here](http://broadinstitute.github.io/picard/explain-flags.html))
3. reference sequence name (ie. chromosome)
4. .
5. .
6. .
7. Ref. name of the mate/ next segment
8. position of the mate/ next segment
9. Observed template length
10. segment sequence (query sequence)
11. ASCIII of the phred-scale 33 score

