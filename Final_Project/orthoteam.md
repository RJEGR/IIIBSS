# Inference on the evolutionary history of embryo cerebral related genes in vertebrates using a novel method.

## General objetive

- Infer the evolutionary history of genes associated with the development of the brain in the embryonic stage of three vertebrates through an orthologous gene function analysis.
- Inferir la histora evolutiva de genes asociados al desarrollo del cerebro en la etapa embrionaria de tres vertebrados mediante un análisis de funcoón de genes ortólogos.

## Particular objectives

- From [ensemble.org](http://www.ensembl.org/info/data/ftp/index.html) download the transcripts coding sequences resulting from Ensembl gene predictions (CDS) for the spices *Mus musculus*, *Gallus gallus* and *Homo sapiens*.
- Build an orthologous-genes graph through the use of the tool *proteinortho*.
- Identify in the graph those genes related with the development of the brain in the embryonic stage.
- Identificar en el grafo a los genes que están relacionados con el desarrollo del cerebro en la etapa embrionaria.
- Determinar eventos evolutivos (especiasion y duplicación) de los genes 
- ... para hacer la reconciliación con el árbol de especies de los eventos evolutivos.
- Comparar resultados con los árboles de eventos generados por otros métodos.



- 







## Methods

a. Download the data set

a. Annotate with gene ontology terms the cds genes from human, mouse and Gallus.

b. subset the genes by the DEGs list

### a. Downloading datasets

First download the data from the ftp (using the peptide instead)

```bash
wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/pep/*pep.all.fa.gz
wget  ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/pep/*pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/gallus_gallus/pep/*pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/danio_rerio/pep/*pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/petromyzon_marinus/pep/*pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/pep/*pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/xenopus_tropicalis/pep/*pep.all.fa.gz
```

Lets label the sets of data

```bash
for i in *.fa; 
	do 
awk '/^>/{print $1, $4, $7; next}{print}' < $i > $(echo $i | cut -d"." -f1).rm; 
done
# second 
for i in *.rm; 
    do
cat $i | sed 's/gene://g' | sed 's/gene_symbol://g' > ${i%.rm}.tmp
done

# rename headers (in process)
# Use double quotes to make the shell expand variables while preserving whitespace:

sed -i "s/$var1/ZZ/g" "$file"

for i in *tmp;
	do
	S = $(echo $i | cut -d"." -f1)
sed "s/^\(>.*\)$/\1 $S/" ${S}.fasta) > "$S".fst
done

# lets rename headers


sed 's/^\(>.*\)$/\1 Danio_rerio/' Danio_rerio.tmp > Danio_rerio.fst
sed 's/^\(>.*\)$/\1 Xenopus_tropicalis/' Xenopus_tropicalis.tmp > Xenopus_tropicalis.fst
sed 's/^\(>.*\)$/\1 Drosophila_melanogaster/' Drosophila_melanogaster.tmp > Drosophila_melanogaster.fst
sed 's/^\(>.*\)$/\1 Petromyzon_marinus/' Petromyzon_marinus.tmp > Petromyzon_marinus.fst
sed 's/^\(>.*\)$/\1 Gallus_gallus/' Gallus_gallus.tmp > Gallus_gallus.fst
sed 's/^\(>.*\)$/\1 Homo_Sapients/' Homo_sapiens.tmp > Homo_sapiens.fst
sed 's/^\(>.*\)$/\1 Mus_musculus/' Mus_musculus.tmp > Mus_musculus.fst

rm .*tmp
cat *fst > cds.fasta



```

And also filter sequence by length

```python
from Bio import SeqIO
from Bio.Alphabet import generic_dna

filename = "./cds.fasta"

L= [[record.id,len(record)] for record in SeqIO.parse(filename, "fasta", generic)]
o = open("cds.fasta.lengths", "w") # create out file
# write the file by pipe the replace and upper function through your data content
for i in range(len(L)):
    o.write(str(L[i][0])+ " " +str(L[i][1])+'\n')

```

Then:

```bash
awk '{if($2 < 200) print $1 "\t" $2}' cds.fasta.lengths > min.lengths
awk '{if($2 >= 200) print $1}' cds.fasta.lengths  > max.lengths
```

And also prepare this input in order to plot a density plot colored by species

```bash
grep "^>" cds.fasta |  awk '{print $2}' > species
paste cds.fasta.lengths species > cds.fasta.lengths.species && rm species
```

And open in R

```R

```



And run the follow script from the bash `python get_seqs.py max.lengths cds.fasta input.fasta`

```python
#!/usr/bin/python
"""
# Execute by: python seq_extractor.py readsList fastafile outputfile
"""
from Bio import SeqIO
import sys

readsList = open(sys.argv[1], 'rU')
fastafile = sys.argv[2]
outputfile = open(sys.argv[3], 'w')

wanted = set()
with readsList as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.add(line)

fasta_sequences = SeqIO.parse(open(fastafile),'fasta')

with outputfile as i:
    for seq in fasta_sequences:
        if seq.id in wanted:
            SeqIO.write([seq], i, "fasta")

readsList.close()
outputfile.close()
```

And, figure out about those sequence less than 200 pb in size:

```bash
 awk '{print $1}'  min.lengths > minList
 python get_seqs.py minList cds.fasta minLengths.fasta
```

Lets figure out the proportions of this minLengths

```bash
grep "^>" minLengths.fasta | awk '{print $2}' | sort | uniq -c
```

| Seq length < 200 | Specie                  |
| ---------------- | ----------------------- |
| 1156             | Danio_rerio             |
| 468              | Drosophila_melanogaster |
| 366              | Gallus_gallus           |
| 10290            | Homo_Sapients           |
| 4650             | Mus_musculus            |
| 377              | Petromyzon_marinus      |
| 53               | Xenopus_tropicalis      |

Finally:

```bash
 awk '{print $1}'  max.lengths > maxList
 python get_seqs.py maxList cds.fasta input.fasta
 #
 grep "^>" input.fasta | awk '{print $2}' | sort | uniq -c
```

Now lets use the `input.fasta` file for further analysis.

Also, use the lengths files in order to filter sequence per species

```bash
awk '{if($2 >= 200) print $1}' cds_lenghts.txt
```

Then, separate by organism

```bash
for i in Danio_rerio Drosophila_melanogaster Gallus_gallus Homo_Sacdpients Mus_musculus Petromyzon_marinus Xenopus_tropicalis;
do
grep $i input.fasta | sed 's/>//g' | awk '{print $1}'  > ${i}.list;
done
```

and get, seqs

```bash

for i in Danio_rerio Drosophila_melanogaster Gallus_gallus Homo_Sacdpients Mus_musculus Petromyzon_marinus Xenopus_tropicalis;
do
python ./get_seqs.py ${i}.list cds.fasta ${i}.fasta
done
```

prove we have the same number of sequence than table above

```bash
grep -c "^>" *.fasta
```



## Annotate the sequence

lets use some tools in order to annotate 

```bash
 module load hmmer/3.1.2 trinotate/3.1.1 ncbi_blast+/2.7.1
```

Finally run the annotation: 

>  turn on the `qlogin` option before everything

1. First , download dataset

```bash
nohup ./Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate &
# in case of: install `cpan DBI.pm`
```

2. Also prepare both  databases `uniprot_sprot.pep` and `Pfam-A.hmm.gz`

```bash
makeblastdb -in uniprot_sprot.pep -dbtype prot
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

2. Time to run the annotation step:

```bash
nohup blastp -query Mus_musculus_DEGS.fasta -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > swissprot.blastp.outfmt6 &
```

3. and Pfam method by hide markov model:

```bash
nohup hmmscan --cpu 4 --domtblout TrinotatePFAM.out Pfam-A.hmm Mus_musculus_DEGS.fasta > pfam.log &
```



And also use this tools

```bash
wget https://github.com/Trinotate/Trinotate/archive/Trinotate-v3.1.1.tar.gz # <---
wget https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-v5.3.0.tar.gz

```

## Subseting up-regulated genes expressed in Mice

in R lets get the list of up genes

```R
x <- read.csv("mouse_degs.csv")
table(x$mouse_deg)
#Down   Up
# 1083 1174
x[ x$mouse_deg == "Up",]
up <- x[ x$mouse_deg == "Up",]

#
write.table(x[,c(2, 5:7)], 
                file = "mouse_degs_up.txt", sep = "\t",
                 row.names = FALSE)
```

## Running proteinorth

You can easly run by: `./proteinortho5.pl -p=blastp+ -project="yourprojectname" -singles -clean -graphe` (define whrer is of nucleotide `blastn`or peptide `blastp` mode) and aslo perfor in a specific step analysis:

- 1 → generate indices
- 2 → perform pairwise analyses
- 3 → perform clustering
- 0 → perform all steps

In case you run it in a server mode:

> module load gcc-5.4.0 # in  omica
> module avail ncbi_blast+/2.7.1 # in lavis ::: use screen to for task manager :) screen mycode &

```bash
./proteinortho5.pl -p=blastp+ -project="orthopep" -singles -verbose -clean -step=1 -graph *.fst
```

Then , perform in a parallel mode the pairwise alignment:

```bash

for i in $(seq 4)
do
echo "./proteinortho5.pl -p=blastp+ -project="orthopep" -singles -verbose -clean -step=2 -job=${i}/4 -graph *.fst"
done 
```

And final step: Clustering,

```bash
#:::::: finaly clustering
./proteinortho5.pl -p=blastn+ -project="orthoteam" -singles -verbose -clean -step=3 -graph *.fst
```

follow steps for:

1. Decomposition
2. reconcilation (https://github.com/Bloodfield/reconciliation)