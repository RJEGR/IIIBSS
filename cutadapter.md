```
cutadapt -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG \
		-A ACACTCTTTCCCTACACGACGCTCTTCCGATCT \
		-q 20 -O 10 -m 25 \
		-o SRR359063_1.trimmed.fastq.gz -p SRR359063_2.trimmed.fastq.gz \
		../raw/SRR359063_1.fastq.gz ../raw/SRR359063_2.fastq.gz
```

