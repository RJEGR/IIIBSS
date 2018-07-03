#library("GO.db")
#pfam <- read.delim("TrinotatePFAM.out")
blastp <- read.csv("swissprot.blastp.outfmt6.tmp", sep="\t", header=F)
colnames(blastp)[1:2] <- c("ENSEMBLPROT", "SYMBOL")

# bitr: Biological Id TranslatoR
# clusterProfiler provides bitr and bitr_kegg for converting ID types. Both bitr and bitr_kegg support many species including model and many non-model organisms.
# ref: https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#bitr-biological-id-translator
library("org.Mm.eg.db")
keytypes(org.Mm.eg.db)

head(keys(org.Mm.eg.db, keytype="ENSEMBLPROT"))

y <- stringr::str_replace(blastp[,1], ".[0-9]$", "")
# y <- read.csv("DEGs.list")
colnames(y) <- NULL

go = clusterProfiler::bitr(y, fromType="ENSEMBLPROT", toType=c("GO", "ONTOLOGY", "GENENAME", "ENSEMBLTRANS" ), OrgDb="org.Mm.eg.db")
str(go)

# prepare bridge file 
db <- out <- read.table("genes.pep.csv", header=F, sep = "\t", stringsAsFactors=FALSE)
colnames(db) <- c("ENSEMBLGENE", "ENSEMBLPROT")
z <- merge(go, db, suffix = c("ENSEMBLPROT"), all=FALSE)
go <- z
table(go['ONTOLOGY'])
MF <- go[go['ONTOLOGY']=="MF",c(6,7,1,5, 4)]
CC <- go[go['ONTOLOGY']=="CC",c(6,7,1,5, 4)]
BP <- go[go['ONTOLOGY']=="BP",c(6,7,1,5, 4)]


MF <- MF[!duplicated(MF$ENSEMBLPROT),]
CC <- CC[!duplicated(CC$ENSEMBLPROT),]
BP <- BP[!duplicated(BP$ENSEMBLPROT),]

savefile <- rbind(CC,BP, MF)

write.table(savefile, file=paste0("Annot.ontology.db.csv"), 
                sep = "\t", col.names = TRUE,
                row.names = TRUE, quote = FALSE)

write.table(MF, file=paste0("Annot.MF.db.csv"), 
                sep = "\t", col.names = TRUE,
                row.names = FALSE, quote = FALSE)

write.table(CC, file=paste0("Annot.CC.db.csv"), 
                sep = "\t", col.names = TRUE,
                row.names = FALSE, quote = FALSE)

write.table(BP, file=paste0("Annot.BP.db.csv"), 
                sep = "\t", col.names = TRUE,
                row.names = FALSE, quote = FALSE)
# /Users/cigom/Documents/orthoteam
newh <- read.csv("cds.4.mus_musculus.headers", sep=" ", header=F, stringsAsFactors=FALSE)
colnames(newh)[4:5] <- c("ENSEMBLTRANS", "ENSEMBLGENE")
newh2 <- merge(go, newh, suffix = c("ENSEMBLTRANS"), all=FALSE)
colnames(newh2)[8:10] <- c("ORTHOTRANS", "ORTHOGENE", "GENENAME2")
CC2 <- newh2[newh2['ONTOLOGY']=="CC",]
CC2 <- newh2[!duplicated(newh2$ORTHOTRANS),c(8,9,10,7)]


write.table(CC2, file=paste0("Annot.CC.mm.csv"), 
                sep = "\t", col.names = TRUE,
                row.names = FALSE, quote = FALSE)

library(DT)
 z <- data.frame(x2)
 #z$pfam <- paste0('<a href="http://pfam.xfam.org/family/', z$pfam, '">', z$pfam,  '</a>')
datatable(db , escape=1, options = list( pageLength = 25 ), filter="top")
ENSMUST00000131374.7	ENSMUSG00000008683.16	Rps15a
ENSMUST00000119705.7	ENSMUSG00000037426.17	Depdc5 GATOR complex protein DEPDC5 
increased circulating interferon-gamma level	MGI	4 (Unknown)	Igtp_tm1Gvw
increased circulating interleukin-12b level	MGI	4 (Unknown)	Igtp_tm1Gvw
increased susceptibility to parasitic infection
ENSMUST00000106588.7	ENSMUSG00000008683.16	Rps15a 40S ribosomal protein S15a



      datatable(data, options = list(), class = "display", callback = JS("return table;"),
         rownames, colnames, container, caption = NULL, filter = c("none",
             "bottom", "top"), escape = TRUE, style = "default", width = NULL,
         height = NULL, elementId = NULL, fillContainer = getOption("DT.fillContainer",
             NULL), autoHideNavigation = getOption("DT.autoHideNavigation",
             NULL), selection = c("multiple", "single", "none"), extensions = list(),
         plugins = NULL, editable = FALSE)