
for i in *.fa; 
	do 
awk '/^>/{print $1, $4, $7; next}{print}' < $i > $(echo $i | cut -d"." -f1).rm; 
done
#
for i in *.rm; 
    do
cat $i | sed 's/gene://g' | sed 's/gene_symbol://g' > ${i%.rm}.tmp
done

cat *.fa > cds.1.fasta
awk '/^>/{print $1, $4, $7; next}{print}' cds.1.fasta > cds.2.fasta
grep "^>" cds.2.fasta | sed 's/^>//g' > cds.2.fasta.header
cat cds.2.fasta.header | sed 's/gene://g' | sed 's/gene_symbol://g' > cds.3.fasta.header

awk '{print $1}'  cds.3.fasta.header| sed 's/.[0-9]$//' > transcript
awk '{print $2}' cds.3.fasta.header | sed 's/.[0-9]$//' > genes

cds.4.mus_musculus.headers
# in r

#### ::::
sed 's/^\(>.*\)$/\1 Danio_rerio/' Danio_rerio.tmp > Danio_rerio.fst
sed 's/^\(>.*\)$/\1 Xenopus_tropicalis/' Xenopus_tropicalis.tmp > Xenopus_tropicalis.fst
sed 's/^\(>.*\)$/\1 Drosophila_melanogaster/' Drosophila_melanogaster.tmp > Drosophila_melanogaster.fst
sed 's/^\(>.*\)$/\1 Petromyzon_marinus/' Petromyzon_marinus.tmp > Petromyzon_marinus.fst
sed 's/^\(>.*\)$/\1 Gallus_gallus/' Gallus_gallus.tmp > Gallus_gallus.fst
sed 's/^\(>.*\)$/\1 Homo_Sapients/' Homo_sapiens.tmp > Homo_sapiens.fst
sed 's/^\(>.*\)$/\1 Mus_musculus/' Mus_musculus.tmp > Mus_musculus.fst
# rename headers (in process)
for i in *tmp;
	do
S = echo $i | cut -d"." -f1
sed "s/^\(>.*\)$/\1 $S/" ${S}.fasta) > "$S".fst
done

Use double quotes to make the shell expand variables while preserving whitespace:

sed -i "s/$var1/ZZ/g" "$file"