
# ensembl: ENSMUSP00000121394.1 
# source: pep 
# chromosome:GRCm38:2
# startPosition:181811197
# endPosition:181820099:1 
# gene:ENSMUSG00000010505.16 # <--- this id match to the DEGs list
# transcript:ENSMUST00000135744.2 
# gene_biotype:protein_coding 
# transcript_biotype:protein_coding 
# gene_symbol:Myt1 
# description:Myelin transcription factor 1  [Source:UniProtKB/Swiss-Prot;Acc:Q8CFC2]

# in bash
cat Mus_musculus.GRCm38.PEP.all.headers | tr ":" " "  > Mus_musculus.GRCm38.PEP.all.headers2 
# in r
y <- read.delim("Mus_musculus.GRCm38.PEP.all.headers2 ", sep=" ", header=FALSE, row.names = NULL)
savefile <- y[,c(1,10,12,18, 5:7)]
colnames(savefile) <- c("ensembl", "gene", "transcript", "gene_symbol", "chr", "start", "end")

write.table(savefile, file="peptide_db.csv", 
                sep = "\t", 
                col.names = FALSE, 
                row.names = FALSE, quote = FALSE)


# python get_seqs.py DEGsList  output.fasta


