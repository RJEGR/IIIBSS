#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
#
# @uthor:      rgomez@cicese.edu.mx
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#
# This file contains code that helps to wragling fasta files
#
#------------------------------
from Bio import SeqIO
from Bio.Alphabet import generic_dna

filename = "./cds.fasta"

L= [[record.id,len(record)] for record in SeqIO.parse(filename, "fasta", generic_dna)]
o = open("cds_lengths.txt", "w") # create out file
# write the file by pipe the replace and upper function through your data content
for i in range(len(L)):
    o.write(str(L[i][0])+ " " +str(L[i][1])+'\n')


# second approach


filename = "./cds.fasta"
i = open(filename, "r")

def seqLen(filename):
    f1 = open(filename, "r" )
    o = open("cds_lenghts.txt", "w")
    L= [[record.id,len(record)] for record in SeqIO.parse(filename, "fasta", generic_dna)]
    for i in range(len(L)):
        o.write(str(L[i][0])+ " " +str(L[i][1])+'\n')
    o.close
    return print("Analysis doing ...")


# ok, continue by filtering the data

tools = /cm/shared/apps/
/cm/shared/apps/ncbi_blast+/2.7.1/bin
/cm/shared/apps/hmmer/3.1.2/src #<-
/cm/shared/apps/hmmer/3.1.2
/cm/shared/apps/trinotate/3.1.1/Trinotate-Trinotate-v3.1.1/:
/cm/shared/apps/


nohup ./Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate &


autoTrinotate.pl --Trinotate_sqlite Trinotate.sqlite \
			--transcripts input.fasta \
			--conf ./conf.txt



/media/public/orthoTeam/fastas/cds.fasta.lengths.species

for i in Danio_rerio Drosophila_melanogaster Gallus_gallus Homo_Sacdpients Mus_musculus Petromyzon_marinus Xenopus_tropicalis;
do
grep $i input.fasta | sed 's/>//g' | awk '{print $1}'> ${i}.list
done

# Finally

for i in Danio_rerio Drosophila_melanogaster Gallus_gallus Homo_Sacdpients Mus_musculus Petromyzon_marinus Xenopus_tropicalis;
do
python ./seq_len/get_seqs.py ${i}.list.tmp cds.fasta ${i}.fasta
done


grep "^>" cds.fasta |  awk '{print $2}' > species
paste cds.fasta.lengths species > cds.fasta.lengths.species && rm species

x <- read.csv("input.file", sep=" ", header=FALSE)
colnames(x) <- c("length", "Specie")
library(ggpubr)

min = x[x$length < 200, ]

den <- ggdensity(x, x = "length", y = "..density..",
          #add = "mean", 
          rug = TRUE,
          #color = " Specie", 
          fill = "Specie",
          #palette = c("#FC4E07", "#00AFBB", "#E7B800")
          )

den + facet_grid(.~ Specie)

#

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

library(tidyverse)

peptide <- read.csv("peptide_db.tmp2", sep=" ", header=F, stringsAsFactors=F)
DEGs <- read.csv("DEGsList", header=F, stringsAsFactors=F)
colnames(peptide) <- c("peptide", "genes")
colnames(DEGs) <- "genes"

savefile <- merge(DEGs, peptide, suffix = c("genes"), all=FALSE)
nomatch <- anti_join(DEGs, peptide, by ="genes")
write.table(savefile, file="DEGS2peptide.csv", 
                sep = "\t", col.names = FALSE,
                row.names = FALSE, quote = FALSE)

write.table(nomatch, file="DEGS2nomatch.csv", 
                sep = "\t", col.names = FALSE,
                row.names = FALSE, quote = FALSE)

# in bash
awk '{print $2}' DEGS2peptide.csv > DEGS2peptide.csv.tmp
python3.6 get.seqs.py DEGS2peptide.csv.tmp Mus_musculus.GRCm38.PEP.nucl.all.fa Mus_musculus_DEGS.fst


One of the easy ways is not to store the output in a variable, but directly iterate over it with a while/read loop.

Something like:

grep xyz abc.txt | while read -r line ; do
    echo "Processing $line"
    # your code goes here
done
There are variations on this scheme depending on exactly what you're after.

If you need to change variables inside the loop (and have that change be visible outside of it), you can use process substitution as stated in fedorqui's answer:


# other stuffs

wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/cds/*.gz
wget  ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/cds/*.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/gallus_gallus/cds/*.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/danio_rerio/cds/*.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/petromyzon_marinus/cds/*.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/cds/*.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/xenopus_tropicalis/cds/*.gz

gunzip *.gz



for i in *fa; 
	do 
awk '/^>/{print $1, $4; next}{print}' < $i > $(echo $i | cut -d"." -f1).tmp; 
done

sed 's/^\(>.*\)$/\1 Danio_rerio/' Danio_rerio.tmp > Danio_rerio.fst
sed 's/^\(>.*\)$/\1 Xenopus_tropicalis/' Xenopus_tropicalis.tmp > Xenopus_tropicalis.fst
sed 's/^\(>.*\)$/\1 Drosophila_melanogaster/' Drosophila_melanogaster.tmp > Drosophila_melanogaster.fst
sed 's/^\(>.*\)$/\1 Petromyzon_marinus/' Petromyzon_marinus.tmp > Petromyzon_marinus.fst
sed 's/^\(>.*\)$/\1 Gallus_gallus/' Gallus_gallus.tmp > Gallus_gallus.fst
sed 's/^\(>.*\)$/\1 Homo_Sapients/' Homo_sapiens.tmp > Homo_sapiens.fst
sed 's/^\(>.*\)$/\1 Mus_musculus/' Mus_musculus.tmp > Mus_musculus.fst

module load gcc-5.4.0 # omica
module avail ncbi_blast+/2.7.1 # in lavis ::: use screen to for task manager :) screen mycode &
# then run in bash
# start here indexitn
./proteinortho5.pl -p=blastp+ -project="orthopep" -singles -verbose -clean -step=1 -graph *.fst
# then the pair-wise blast
./proteinortho5.pl -p=blastp+ -project="orthopep" -singles -verbose -clean -step=2 -job=1/4 -graph *.fst 
#**Step 2**
#Running blast analysis: 100% (5/5, 21 in total)
#[OUTPUT] -> written to orthoteam.blast-graph_1_4
./proteinortho5.pl -p=blastn+ -project="orthopep" -singles -verbose -clean -step=2 -job=2/4 -threads=168 -graph *.fst 
#Running blast analysis: 100% (5/5, 21 in total)
#[OUTPUT] -> written to orthoteam.blast-graph_1_4
./proteinortho5.pl -p=blastn+ -project="orthopep" -singles -verbose -clean -step=2 -job=3/4 -threads=168 -graph *.fst 
#Running blast analysis: 100% (5/5, 21 in total)
#[OUTPUT] -> written to orthoteam.blast-graph_3_4
./proteinortho5.pl -p=blastn+ -project="orthopep" -singles -verbose -clean -step=2 -job=4/4 -threads=168 -graph *.fst
#Running blast analysis: 100% (6/6, 21 in total)
#[OUTPUT] -> written to orthoteam.blast-graph_4_4

#:::::: finaly clustering
./proteinortho5.pl -p=blastn+ -project="orthoteam" -singles -verbose -clean -step=3 -graph *.fst
 82198 orthoteam.blast-graph_1_4
  16768 orthoteam.blast-graph_2_4
   5905 orthoteam.blast-graph_3_4
   4621 orthoteam.blast-graph_4_4

scp rgomez@omica:/home/rgomez/bin/proteinortho_v5.16b/orthoteam/run/nucleotide/orthoteam.proteinortho* .
# path /home/rgomez/bin/proteinortho_v5.16b/orthoteam/run

for i in $(seq 4)
do
echo "./proteinortho5.pl -p=blastp+ -project="orthopep" -singles -verbose -clean -step=2 -job=${i}/4 -graph *.fst"
done 


grep  -Fxv -f Mus_musculus.db2 DEGS2peptide.csv  | while read -r line ; do
    echo "Processing $line"
    # your code goes here
done

grep -Fxv -f first-file.txt second-file.txt

x <- read.csv("Mus_musculus.db3", sep=" ", header=F, stringsAsFactors=F)


for gene in $(cut -f1 DEGS2peptide.csv)
do   
    grep $gene Mus_musculus.db2 | column -t >> output.2
done


out <- read.table("mus_musculos.output.csv", header=F, sep = " ", stringsAsFactors=FALSE)
colnames(out) <- c("transcript", "isoform", "gene_name", "gene")

# use the go annotation in order to merge datasets

#

