dir <- c("/Users/cigom/Documents/orthoteam/")
library("treeio")
library("ggtree")

species <- read.csv(paste0(dir,'species.id'), header=F, stringsAsFactors=F, sep=" ")
annot <- read.csv(paste0(dir, 'Annot.ontology.db.csv'), header = TRUE, row.names = 1, sep='\t')
colnames(species) <- c("ENSEMBLTRANS", "ENSEMBLGENE", "GENBANK", "SPECIES")


dplyr::as.tbl(annot) %>% 
  filter(paste(ENSEMBLTRANS) %in% paste(species$ENSEMBLTRANS, species$SPECIES))



for gene in $(awk '{print $1}' species.id)
do
        echo $(grep $gene Annot.ontology.db.csv) >> test
done


# Load tree :::::
netwick <- args[0]
tree <- treeio::read.newick(netwick)


tr <- read.tree(text = "((a,(b,c)),d);")
genus <- c("Gorilla", "Pan", "Homo", "Pongo")
species <- c("gorilla", "spp.", "sapiens", "pygmaeus")
geo <- c("Africa", "Africa", "World", "Asia")
d <- data.frame(label = tr$tip.label, genus = genus,
                species = species, geo = geo)
d

ggtree(tr) %<+% d + xlim(NA, 5) +
    geom_tiplab(aes(label=paste0('italic(', genus, ')~bolditalic(', species, ')~', geo)), parse=T)



tree <- read.tree("~/Downloads/phylotree.txt")
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)

ggtree(tree)+
    geom_point(aes(shape=isTip, 
    color=isTip), size=3) + geom_tiplab2(color='blue')

ggtree(tree, layout="circular") + geom_tiplab2(color='blue')
ggtree(tree) +
    geom_text(aes(label=node), hjust=-.3)


geom_tiplab(mapping = NULL, hjust = 0, align = FALSE,q

           seq          cat
1         E71T     Database
2        JS100     Outbreak


ggtree(tree) %<+% savefile + 
  geom_tiplab(aes(fill = ONTOLOGY),
              color = "black", # color for label font
              geom = "label",  # labels not text
              #label.padding = unit(0.15, "lines"), # amount of padding around the labels
              label.size = 0) + # size of label border


  theme(legend.position = c(0.5,0.2), 
        legend.title = element_blank(), # no title
        legend.key = element_blank()) # no keys



filaname <- 
netwick <- treeio::read.newick()

