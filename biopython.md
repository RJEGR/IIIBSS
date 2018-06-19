## Disclaim

Where to start a project?

- ensemble.org
- Ncbi.nlm.nih.gov
- Gene ontology.org

other databases:

- Epigenetics
- regulatory data
- wikipedia offers a good source list of databases [here](https://en.wikipedia.org/wiki/List_of_biological_databases) 

And Tools:

- galaxy.com
- Improve skills in the bash [linux](https://www.learnenough.com/command-line-tutorial) 

# Biopython

...

How to install it

> in mac first `sudo easy_install pip` before it

```python
pip install biopython
pip install --upgrade biopython
pip uninstall bipython
```

Then

```python
from Bio.Seq import Seq
```

Let's make an biological example

```python
mySeq = Seq("ACGGTAGATACACGGTAGATAC")
print(mySeq)
print(mySeq.alphabet)
#Alphabet()
```

How do we get the complement of this sequence:

```python
print(mySeq.complement())
# TGCCATCTATGTGCCATCTATG
```

```python
print(mySeq.reverse_complement())
GTATCTACCGTGTATCTACCGT
#
print(mySeq.reverse_complement().reverse_complement())
GTATCTACCGTGTATCTACCGT
```

### Sequence string

> Concept: an arragement of characters

Let's manage this characters

```python
for nucleotide in mySeq:
    print (nucleotide)
# and also:    
print(mySeq[5])
print("Length of your sequence is: ", len(mySeq))
```

Also can count 

> WARNING: this is not an overlaping method!!

```python
Seq("ACGTTGCA").count("AC")
# 1
```

So, lets find at least one time the statement from above:

```python
Seq("AC") in Seq("ACGTTGCA")
# true
```



Great!, let's do some exercise:

In the next `seq` object count the patterns:

Seq = ACGGTAGCTAGCAACTAGATAGACA

```python
seq = Seq("ACGGTAGCTAGCAACTAGATAGACA")
seq.count("GA") # 2
seq.count("AC") # 3
seq.count("TA") # 4
seq.count("TAG") # 4
seq.count("GATA") # 1
```

## count the GC %

in the manual way

```python
mGCpercentage = 100 * (seq.count("G") + seq.count("C")) / len(seq)

print(mGCpercentage)
# 44.0
```



```python
from Bio.SeqUtils import GC
GCpercentage = GC(seq)
print(GCpercentage)
# 44.0
```

## Slicing a sequence

Let's get a slice of the sequence. Create sliced object using the next Sequence `GATCGATGGGCCTATATAGGATCGAAAATCGC` . Split in the half the sequence and save them in 2 separate objects (variables).

```python
seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
z = len(seq)

sliced1 = seq[:int(z / 2)]
sliced2 = seq[int(z / 2): ]

print("Both", sliced1 , "and", sliced2, "are halfs from: ", sliced1 + sliced2
```

### Uper and lower method (to serch motifs:

```python
if "ATATAGG" in seq:
    print("Found")
else:
    print("Not Found")
    
 # lets secure all string are from same case (lower/uper)

if "ATATAGG" in seq.upper():
    print("Found")
else:
    print("Not Found")
```

### Sequence stride

```python
seq[0::3] # fist nucleotide of each codon 
seq[1::3] # second nucleotide of each codon 
seq[3::3] # last nucleotide of each codon

```

## Getting sequence from a fasta file

```python
from Bio import SeqIO
openedFastaFile = SeqIO.parse("./ls_orchid.fa", "fasta")
```

Now, that we have the fasta object, lets get metrix from the sequences:

```python
for record in openedFastaFile:
    ## ID of the read
    print( record.id)
    
    ## Full Sequence
    print( record.seq)
    
    ## Short sequene print
    print( repr( record.seq))
    
    #Length of the sequence
    print( len( record))
```

let's continue with some skills

```python
from Bio.seq import Seq
from Bio.Alphabet import IUPAC

template_dna = Seq("GCCGCTGAAACT", IUPAC.unambiguous_dna)
coding_dna = template_dna.reverse_complement()
messenger_rna = coding_dna.transcribe()
```

reverse transcription 

```python
protein = messenger_rna.translate()
```

> to choose open reading frame tables go to genetic code [here](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

Example:

```python
messenger_rna = essenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG", IUPAC.unambiguous_rna)

messenger_rna.translate()
# Seq('MAIVMGR*KGAR*', HasStopCodon(IUPACProtein(), '*'))

messenger_rna.translate(table= "Vertebrate Mitochondrial")
# Seq('MAIVMGRWKGAR*', HasStopCodon(IUPACProtein(), '*'))

messenger_rna.translate(table = "Bacterial")
# Seq('MAIVMGR*KGAR*', HasStopCodon(IUPACProtein(), '*'))
```

