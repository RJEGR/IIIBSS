--------------------------------
Introduction to programming in R
--------------------------------
  
*by Katja Nowick, University of Leipzig*
  
Interactive R
-------------
  
# The easiest way to use R is in its interactive mode. You start R in its interactive mode by 
# simply typing "R" at the command line::
  
R

# You will see how R starts and welcomes you. Each line begins with a ">". In this mode we 
# can "talk" to R. We can ask it some questions, press "Enter" and R will answer us. 
# For instance, we can use R as a calculator:
  
3+5
13-8
2*(7+2)

# Or let's do something a little bit more complicated:

sqrt(16)
log(4)
sin(90)

# Ok. That's easy, isn't it?

# In the last three lines you used *functions*. Most of the time when you are programming 
# in R you will use functions. You always recognize functions on their parentheses. 
# Functions are operations that always work the same way. In the parentheses you write 
# the arguments with which you want to run the function. 
# In our example above, sqrt() is a function. You want to know the square root of 16, 
# so you stick 16 into your function. Likewise, log() and sin() are functions. 
# Here we calculated the logarithm of 4 and the sine of 90.
# We are going to use many more functions today.

Variables
~~~~~~~~~

# A variable is a symbolic name to which a value may be associated. The associated value 
# may be changed. You should always pick a meaningful name for your variables. Never use 
# white space or special characters in a variable name. Never start a variable name with 
# a number. 
# Let's define the variable a:
  
a=2

# If you forgot what a is, simply type "a" again, and R will tell you:
  
a

# Note, there are two ways to assign a value to a variable. We just did it by using 
# the "=" sign. You will also often see that people use "<-" instead.
  
a<-3
a

# Now you see that a equals 3 instead of 2, because we assigned the value 3 to a. 
# So, what we just did, was overwriting the old value of a and assigning it a new value. 
# Some people prefer to use "<-" to clearly indicate that they want to ASSIGN a value 
# to a variable and not TEST if the variable is equal to a certain value. To check 
# if a variable is equal to a certain value, we use "=="
  
a==3
a==2

# R answers us with "TRUE" or "FALSE", respectively.
# "TRUE" and "FALSE" are special values that indicate that a test was passed (TRUE) 
# or not (FALSE). 

# There are different types of variables. Our variable a is a "numeric" variable. 
# You can test of which type a variable is by asking:
  
class(a)

# Another type of variable is character. Let's define some character variables:

name1="Anna"
name2="Felix"

# Note how we always have to use quotes when assigning a character to a variable, 
# otherwise we get an error message.

# We can use our variables to do some calculations:

a*a
a^2
b=6
b/a
a<b
a==b

# A list of all operators can be found here: http://cran.r-project.org/doc/manuals/R-lang.html#Operators

# What happens if we type:

a=b


# Check what the value of a is now:

a

# What happens if we try to do calculations with our character variables:

name1*name2

# Oops, we can of course not do calculations with characters. But we can ask for instance, 
# if the two character variables are the same:

name1==name2

# Ok. So far, we talked about two types of variables, numeric variables and character 
# variables. Let's move on to some more complicated variable types: vectors and matrices.

Vectors
~~~~~~~
  
# A vector is a variable that contains more than one element. 
# These can be numbers or characters.
# Let's say we want to create a vector that contains the numbers from 1 to 9:

Numbers=1:9
Numbers

# Alternative ways to achieve the same are to use the functions "seq" or combine "c":

Numbers=seq(1,9)
Numbers
Numbers=c(1,2,3,4,5,6,7,8,9)
Numbers

# In this example, the first way was of course easier. But "c" is more flexible, and 
# would allow us to assign the numbers in any order we want:

Numbers=c(3,6,2,10,2647849,1,999)
Numbers

# Or to create vectors containing characters:

Names=c("Anna","Felix")
Names

# To get a certain element of your vector, you have to tell R which element you want. 
# We use the square brackets to tell R the index of the element we want. 
# Important: R always starts counting with 1:

Numbers[1]
Numbers[2]
Numbers[5]
Numbers[7]
Names[2]

# We can also ask R for more than one element:

Numbers[2:5]
Numbers[6:3]
Numbers[c(1,2,6)]

# If you only want to get the first or last elements of your vector, use:

head(Numbers)
tail(Numbers)

# This can be useful, if you have very long vectors and want to get an idea what your 
# vector contains. In this case you don't want R to print the whole vector on your screen.

# You may have wondered, if the different kinds of brackets mean something. We have 
# already used "()" and "[]". The square brackets are always for *indexing*, for 
# instance if we want to get certain elements from a vector. We will also use the square 
# brackets in the next part to get certain elements from a matrix. 
# The round parentheses are for *functions*. sqrt(), head(), tail(), c(), seq() are functions. 
# R has many, many, many functions. You can even write your own functions. Within the 
# parentheses your write the arguments, for instance, on which variable you want the 
# function to be performed, how exactly do you want the function to be performed etc. 
# You will see this, when we do more functions.
# There is a third type of brackets: {}. We are using these brackets to define *code blocks*, 
# for instance in loops or if-else statements. We will do this later.

# We can also change the values of the elements in our vectors. Let's go back to our 
# vector Numbers and do some changes:

Numbers
Numbers[2]=24
Numbers
Numbers[5]=231
Numbers

# We can add elements to our vector:

Numbers=c(Numbers, 602, 78,5)
Numbers
Numbers=c(-1,-392,52, Numbers)
Numbers

# And we can delete elements from our vector. To do this we need to tell R the index of 
# the element we want to delete. For instance, to delete the 3rd element of our vector:

Numbers=Numbers[-3]

# We can also do the following:

Numbers=Numbers[2:11]

# Here we overwrite our vector Numbers, saying that we only want the elements from 2 to 12 
# in the new vector also called Numbers. The outcome is similar to deleted elements.

# Some functions that can come in handy if you have a long vector and quickly want to 
# check some properties of the vector:

sum(Numbers)
min(Numbers)
max(Numbers)
mean(Numbers)
sort(Numbers)

# You can check whether the values in your vector meet some condition. Let's say, you want 
# to know if there is an element in your vector that has the value 3:
  
Numbers
Numbers==3

# The answer from R is another vector with TRUEs and FALSE. Each element in the result 
# vector corresponds to one element in the Numbers vector. Mostly the elements in the 
# result vector are FALSE, but the 3rd element is TRUE. This means that the 3rd 
# element in Numbers is equal to 3. Let's try some more tests:

Numbers==10
Numbers==6
Numbers<6
Numbers>=24

# Is there any element in your vector that is less than 6, less than -400, greater than 1000?

any(Numbers<6)
any(Numbers<(-400))
any(Numbers>1000)

# Are all elements less than 10?

all(Numbers<10)

# Which elements are less than 10?

which(Numbers<10)

# The function which() gives you the indices of the elements meeting the criteria.
# What is the value of the elements that are less than 10? Here we need to use indexing 
# to not get the number of the element but the value of the element:

Numbers[which(Numbers<10)]

# An alternative way to get the elements that are less than 10 is this:

Numbers[Numbers<10]

# What happens here internally is that R first creates a vector of TRUEs and FALSE by 
# testing each element of Numbers if it is smaller than 10:

Numbers<10

# This vector of TRUEs and FALSEs is then used for indexing and thus we can retrieve 
# the elements less than 10:

Numbers[Numbers<10]

# By now we already created a lot of variables. Do you still remember all of them? 
# If not, ask R to list them for you:

ls()

# If you want to remove some variables, you can use the function rm(). Often you don't 
# need to remove variables, because you can just assign a new values to them. But some 
# variables might be big and need a lot of memory, e.g. big tables. If you need more 
# free memory, remove some large variables:
  
  rm(a)


Exercise I
~~~~~~~~~~
  
# 1. Create a vector that contains the length of various alternative transcripts for a gene. 
# The lengths are in bp: 526,723,1064,821,697,743,1149,489.
# What is the length of the smallest transcript?
# Are there any transcripts longer than 1000 bp?
# Are all transcripts at least 500 bp long?
# You discovered a new splice variant, that creates a transcript of 762 bp length. 
# Add it to your vector.
# Sort all the transcripts by size.
# You sequenced the orthologous gene in the sister species. Most transcripts have the same 
# size. Create a new vector that is a duplicate of the first vector but for the gene of 
# the sister species.
# In the sister species you discovered that the transcript with 1149 bp contains a stop codon.
# Remove it from the vector.

# 2. You got a list of genes that are positively selected in your species of interest. 
# Create a vector that contains these genes: BAL, CAM, LSD, FUZ, ZERP, DING, NOP, YIN.
# Sort the genes alphabetically and save the sorted gene list in a new vector.
# Have a look at the first and last couple of elements of the sorted gene list vector 
# to see if they are indeed sorted.
# Which gene is the 3rd one in your sorted gene list?
# What is the index of DING in the unsorted gene list?
# Check if the gene "LSD" is in your gene list.

# 3. Your collaborator gives you two lists of genes and wants to know if they are among 
# your positively selected genes (from exercise 2): 
#  a) NOP, PHO, SEN; 
#  b) dns, hro, lamp, bal, cat, krr, zerp, yin, tir, jog, nop. 
# Try out how "%in%" works.
# Can you change the small letters to capital letters? Check the help pages on how the 
# function toupper() works.


Matrices and tables
~~~~~~~~~~~~~~~~~~~
  
# A matrix is a rectangular array of values, similar to a table. 
# You can easily create one from a vector. Let's first make a vector using some functions 
# we learned before:

NumberVector=rep(seq(1,9),2)
NumberVector
NumberMatrix=matrix(NumberVector, nrow=6)
NumberMatrix

# You see that you have to specify the number of rows you want. Alternatively, you could 
# have specified the number of columns:

NumberMatrix=matrix(NumberVector, ncol=6)
NumberMatrix

# You also see that in both cases R fills the matrix column-wise. 
# If you want to change it, you need to set the argument byrow to TRUE (it is FALSE by 
# default):

NumberMatrix=matrix(NumberVector, ncol=6, byrow=TRUE)

# In case you want to check, if a variable you are working with is a vector or a matrix, 
# you can ask:

is.vector(NumberVector)
is.vector(NumberMatrix)
is.matrix(NumberMatrix)

# Most of the time you are going to deal with much bigger matrices. Often they will be 
# tables that you got from somewhere. How do you get them into R? Here we need to learn a 
# new function: read.table(). 
# Have a look at the help page for read.table(). There are a lot of arguments you can set. 
# R is sometimes good in guessing the format of the table you are trying to read. At 
# other times it might mess up your table. Here we are going to read in a table that has 
# column names (header), is tab-separated, has no quotation signs, and we want R to keep 
# strings as strings (and not as factors):

Table1=read.table("Matrix1.txt", header=TRUE, sep="\t", quote="", stringsAsFactor=FALSE) 

# Always double-check that the table and the format of the table looks like what you expected.
# You can print the table to your screen:

Table1

# But if it's a big table you don't want to do this. It's wise to check the size (dimensions) 
# of the table first, before you print it:
  
dim(Table1)

# This gives you the number of 1. rows and 2. columns of your table.
# You can also ask for the number of rows and columns separately:
  
nrow(Table1)
ncol(Table1)

# Just like we have seen for vectors before, you can print the first and last part 
# (first/last 6 rows) of the table:
  
head(Table1)
tail(Table1)

# To get certain elements of your table, we need to use indexing. Remember the 
# square brackets for indexing of vectors? For matrices and tables it is exactly the same, 
# except that we have to give R two indices. R expects the row number to come first 
# and the column number to be specified second:
  
Table1[1,3]
Table1[3,4]
Table1[1,2]

# You can also get complete columns from your table. To do so, just leave the first index 
# blank, e.g.:
  
Table1[,1] 

# This will give you the complete first column.

# To get a complete row, leave the second index blank, e.g.:
  
Table1[1,]

# This will give you the complete first row.

# Of course, you can also get other kinds of subsets from your table, e.g.:
  
Table1[1:3,2:3]
Table1[4:10,]
Table1[6,2:4]
Table1[,1:2]

# Our table has column names (we specified this when reading the table by saying header=TRUE). 
# You might find it easier to work with the column names instead of the index for the column. 
# You need to tell R that you refer to the column name by using the $-symbol:
  
Table1$gene
Table1$brain1[1]
Table1$brain2[2:8]

# A table can also have row names. Just specify this when reading the table. If your 
# row names are in column 1, you would use:
  
Table1=read.table("Matrix1.txt", header=TRUE, sep="\t", quote="", row.names=1, 
                  stringsAsFactor=FALSE)

# What are the dimensions of the table now?

# Examples on how to use the row names:
  
rownames(Table1)
Table1["Pax6",2]
Table1["Pax6",]
Table1[c("Pax6", "Otx2"),]

# Let's get the old table back:

Table1=read.table("Matrix1.txt", header=TRUE, sep="\t", quote="", stringsAsFactor=FALSE)

# Since every column and every row of a matrix/table is a vector, you can use all the 
# functions we have learned for vectors, also for rows and columns of your table. 
# Be careful to only refer to one column or row for doing this, e.g.:

max(Table1[,3])
mean(Table1[,4])
sort(Table1[,1]) 
any(Table1[,2]=="Otx2") 

# And you could ask for all Expression1 values larger or smaller than a certain value. 
# For instance, to get all Expression1 values (column 3) larger than 10, you would ask for:

Table1[,3]>10

# This gives you a vector of TRUEs and FALSE, depending on if the Expression1 value is 
# larger than 10 or not. Often you would like to see the actual values and not the TRUEs 
# and FALSEs. You can achieve this by using:

Table1[Table1[,3]>10,3]

# Note, that we use a condition to create a Boolean vector to select the rows we want: 
# in this example we only want the rows in which the Expression1 value (column 3) is 
# larger than 10. So we use Table1[,3]>10 to select the rows. And we want to see the 
# Expression1 values, so we select column 3.

# If we want to do the same, but see the complete row, we would use:

Table1[Table1[,3]>10,]

# Let's suppose you don't want the Expression1 values larger than 10, but actually the 
# names of the genes that have Expression1 values larger than 10. Just specify which 
# column(s) you want:

Table1[Table1[,3]>10,1]
Table1[Table1[,3]>10,1:2]

# Of course, we can also change the values of elements in the table. Say which element 
# you want to change and what the new value should be. Let's change some expression values:
  
head(Table1)
Table1[5,3]=1
head(Table1)
Table1[3,3:4]=c(10, 5)
head(Table1)

# You can also add columns or rows to your table. To do so, you need to create a vector, 
# which you can then add as a column using cbind() or row using rbind():
  
newRow=c("NEW1", "NEW2", 3, 5)
Table1=rbind(Table1, newRow)


Exercise II
~~~~~~~~~~~
  
# Read in the table Matrix2.txt.
# # How many rows and columns does the table have?
# What are the column names?
# Which genes (column 1) have expression values larger than 30 (column 3)
# Is the gene Med31 in the table?
# How often is Med31 in the table?
# Which Expression 17 values belong to Med31?
# Make a second table that contains the first 6 rows of the first table.
# Make a vector containing the values 1, 2, 3, 4, 5, 6.
# Add this vector as another column to the newly created table.


Statistics
----------
  
# R offers a wealth of options for statistical data analysis. In fact, R is called 
# a "statistical programming language". Many statistical tests are implemented in 
# the standard R packages. And you will find even more tests in *R packages* that can be 
# downloaded from the Comprehensive R Archive network (CRAN): 
# http://cran.r-project.org/web/packages/. What R packages are, we will discuss next. 
# Here, as a teaser on statistics, just some very basic data descriptions and statistical 
# tests:


Descriptive statistics
~~~~~~~~~~~~~~~~~~~~~~
  
# Many functions in R are named very intuitively. We have seen how to use functions 
# earlier today. Remember the () brackets. 
# Let's take the vector with transcript lengths from your Exercise 1 (should look 
# like 526  723 1064  821  697  743 1149  489  762) and calculate some simple summary 
# statistics:

mean(Transcripts)
median(Transcripts)
var(Transcripts)
sd(Transcripts)
summary(Transcripts)
quantile(Transcripts)

# In the next part on differential gene expression we are going to use some selected 
# R packages and you can explore some sophisticated statistical tests they are implementing.



Graphics
--------

# R is also very powerful in terms of graphical data presentation. Talking about graphics 
# can easily fill a complete course day. Here are just a few examples to give you an idea.

# Let's plot Expression1 versus Expression2:
  
x=Table1[,3]
y=Table1[,4]
plot(x,y)

# To add axis labels, use:
  
plot(x,y, xlab="brain3", ylab="heart1")

# If you want red dots instead of black open circles, you can also specify this:
  
plot(x,y, xlab="brain3", ylab="heart1", type="p", pch=20, col="red")

# There are many more ways to specify size, shape, color of dots, draw lines, label 
# your dots etc. Start by checking out the help pages for plot() and then follow the 
# links to get more information.

# There are also many other types of plots that can easily be created with R. For instance 
# a box plot:
  
boxplot(y)
boxplot(y, xlab="heart1", ylab="Expression level")

# Some more graphics examples will follow later.

# You can also plot functions. In this case you have to specify the range in which 
# you want to plot the function using the arguments *from* and *to*. Let's, as an 
# example, plot a normal distribution and a sine function:

plot(dnorm, from=-5, to=5)
plot(sin, from = -2*pi, to = 2*pi)

# For publications you sometimes want to have a figure consisting of multiple panels 
# with plots. To create for instance a figure with 4 plots, you first need to set up 
# your plotting parameters, e.g. specifying that you want them to be plotted in 
# a 2x2 layout. You use the function par() for this:

par(mfrow = c(2,2), pty = "s")

# The pty argument just specifies that the plotting should be squared and independent 
# of the device size. Then you can draw your four plots into this area:

plot(x,y, xlab="brain3", ylab="heart1", type="p", 
pch=20, col="red")
boxplot(y, xlab="heart1", ylab="Expression level")
plot(sin, from = -2*pi, to = 2*pi)
plot(1:10, 1:10)

# Very fancy stuff can be done with gplot and ggplot. Have a look at it when you have time: 
# http://docs.ggplot2.org/current/


Exercise III
~~~~~~~~~~~~

# Some people prefer to have slides with a black background. Then it would look nicer, 
# if also your plot had a black background and if your axes and axis labels would be in white. 
# Hint: First, to change the color of your graphics background use par (bg="black"). 
# Then check the help pages for par() to see how you can change the color for your axes 
# and axis labels.

# For further reading and reference have a look at the language definition of R, which 
# can be found here: http://cran.r-project.org/doc/manuals/R-lang.html


------------------------
Gene Expression Analysis
------------------------

*by Vijaykumar Yogesh Muley, 
*Institute of Neurobiology, Universidad Nacional Autónoma de México, Juriquilla Campus, Queretaro



RNA-Seq
-------
# We are going to work with a dataset consisting of 12 transcriptome data from four 
# tissues of mouse, sequenced by Illumina HiSeq 4000 (taken from Dong J et al. 2012). 
# The reads were mapped to the mouse genome. 
# For each gene the number of reads were counted. The counts are in the file 
# "mouse_exp.txt". We will use this dataset to try out two different R libraries for 
# identifying genes that are differentially expressed between the different tissues  and compare 
# the results. 
# Note, that in this tutorial we will use the terms "genes" and "differentially expressed 
# genes". You can use the same RNA-seq library also for counts on transcripts, exons, ESTs, tags etc. 
# and hence identify differentially expressed transcripts, exons, ESTs, tags etc. 

# We first read in our data with the read counts stored in mouse_exp.txt file, using the function read.table(). 

counts=read.table("mouse_exp.txt", head=T, 
                   sep="\t", quote="", stringsAsFactor=FALSE, row.names="gene")

dim(counts)
head(counts)

# The first column contains the gene names. The following columns contain discrete values 
# (counts), representing the expression levels for the genes. 

# The first 3 samples come from brain tissue, the next nine from heart, liver and lung. 

# We need to create a vector that can be used to assign our data to the respective 
# tissue for differential expression analysis
  
groups=c(1,1,1,0,0,0,0,0,0,0,0,0)

# Let's plot them to see how similar our 12 samples are:

boxplot(counts)
boxplot(log2(counts+0.5))

# I added 0.5 to avoid log2 transformation of zero counts to NaNs
# What do you think? Are the samples of the two groups comparable with each other?

# Now we have prepared our input data and can use the libraries.

Remove non-informative genes
~~~~~

# Lets count genes having more than non-zero counts in less than three samples
# before that we check the total number of genes

dim(counts)
sum_non_zeros <- apply(counts!=0, 1, sum)  
counts <- counts[sum_non_zeros>=3,]
dim(counts)

# How many genes have been removed ?

DESeq
~~~~~
  
# DESeq is based on a model using the negative binomial distribution.
# Start by loading the library:
  
library(DESeq)

# Using the function newCountDataSet() we can create an object that can be used by DeSeq 
# to perform the analysis. It will hold the count information for each gene for each sample:
  
counts_DESeq_obj=newCountDataSet(counts, groups)

# Have a look at the object:
  
counts_DESeq_obj

# You got 22573 features - these stand for the number of genes assessed - and 12 samples. 
# You also got the sample names stored in the object.

# The function counts() gives us access to the matrix with our expression counts:
  
head(counts(counts_DESeq_obj))

# Let's plot them to see how similar our 12 samples are:

boxplot(counts(counts_DESeq_obj))
boxplot(log2(counts(counts_DESeq_obj)+0.5))

# I added 0.5 to avoid log2 transformation of zero counts to NaNs
# What do you think? Are the samples of the two groups comparable with each other?

# For the analysis we need to estimate the effective library size to normalize for. 
# Imagine, a gene has the same number of counts in two samples. But the library size 
# (total number of reads) was twice as high in the first sample. Then we would conclude 
#that the gene was higher expressed in the second sample. You can estimate the size 
# factors from the count data using the function estimateSizeFactors():

counts_DESeq_obj=estimateSizeFactors(counts_DESeq_obj)

# Have a look at the size factors:

sizeFactors(counts_DESeq_obj)

# Next we need to estimate for each condition a function that allows to predict the 
# dispersion (= variance across samples). The core assumption of this method is that the 
# mean is a good predictor of the variance, i.e., that genes with a similar expression 
# level also have similar variance across the samples:

counts_DESeq_obj=estimateDispersions(counts_DESeq_obj)

# Remember our boxplot? Let's see if the normalization helped to make the samples more even:
  
boxplot(log2(counts(counts_DESeq_obj, normalized=TRUE)+0.5))

# Now we can test for differential expression of the genes between the two groups 
# by calling the function nbinomTest(). We provide the counts_DESeq_obj and the names of 
# the groups of our samples to this function. This might take a few minutes:
  
DESeq_DEgenes=nbinomTest(counts_DESeq_obj, "0", "1")  

# Have a look at the variable DESeq_DEgenes you produced:
  
head(DESeq_DEgenes)

# It is a data frame with p values (raw and adjusted), mean values, fold changes, and 
# other results. The Mean is the mean of all 6 samples. Then you got the means for each 
# group (here of the brain and of the other tissue samples). The fold change is calculated 
# dividing the mean of the second group by the mean of the first group, the p-values 
# come from the binomial test. Also adjusted p-values are provided and were calculated 
# using the method of Benjamini and Hochberg, which controls for the False Discovery Rate 
# (FDR).

# You can visualize your result by plotting the log2 fold changes against the base means 
# and coloring in red those genes that are significant at 5% FDR:
  
plot(DESeq_DEgenes$baseMean,DESeq_DEgenes$log2FoldChange,log="x", pch=20, 
     cex=.1,col = ifelse( DESeq_DEgenes$padj < 0.05, "red", "black" ) )

# You can see that you got quite a lot of significantly differentially expressed genes. 
# Genes with higher mean expression seems to lead to higher FC.

# To extract the significantly differentially expressed genes we can use one of the 
# next two statements. They do the same thing, just in the first case we use the column 
# name, in the second case the column number for indexing.:
  
signDESeq_DEgenes=DESeq_DEgenes[DESeq_DEgenes$padj<0.05,]
signDESeq_DEgenes=DESeq_DEgenes[DESeq_DEgenes[,8]<0.05,]




edgeR
~~~~~
# edgeR uses empirical Bayes estimation and exact tests based on the negative binomial 
# distribution to call differentially expressed genes.:
  
library(edgeR)

# For the edgeR library we use the function DGEList() to create an object with the sample 
# and read count information::
  
count_edgeR_obj=DGEList(counts=counts, group=groups)

# Here, when you create this object of class DGEList, R already calculates the library size 
# for you.:
  
count_edgeR_obj

# You see how samples got assigned to groups (species), what the library sizes are, the 
# counts for each gene etc.

# edgeR uses the quantile-adjusted conditional maximum likelihood (qCML) method to estimate 
# the dispersion(s) before comparing two groups of samples. It first determines the common 
# dispersion and then the dispersion for each gene. For the gene-wise dispersion it 
# implements an empirical Bayes strategy for squeezing the gene-wise dispersions towards 
# the common dispersion. This may take a few seconds:
  
count_edgeR_obj=estimateCommonDisp(count_edgeR_obj)
count_edgeR_obj=estimateTagwiseDisp(count_edgeR_obj)  

# The edgeR test for differential expression is similar to a Fisher's exact test and is 
# based on the qCML method. This may take a few seconds:

edgeR_DEgenes=exactTest(count_edgeR_obj)  
edgeR_DEgenes

# The edgeR_DEgenes object contains multiple elements. The first one is the table with 
# the fold changes (logFC), logCPM, and p-values for each gene. By default, the p-values 
# in this table are adjusted for multiple testing by the Benjamini-Hochberg method.

# To see the top differentially expressed genes use the function topTags():

topTags(edgeR_DEgenes)

# To get access to the table of the edgeR_DEgenes object for later analysis use 
# edgeR_DEgenes$table:

edgeR_DEgenesTable=edgeR_DEgenes$table
head(edgeR_DEgenesTable)

# You can visualize your result by plotting the log2 fold changes against the base means 
# and coloring in red those genes that are significant at 5% FDR:

plot(edgeR_DEgenesTable$logCPM,edgeR_DEgenesTable$logFC,log="x", pch=20, 
     cex=.1,col = ifelse( edgeR_DEgenesTable$PValue < 0.05, "red", "black" ) )


# Then a table with the significant genes can be extracted:

signedgeR_DEgenes=edgeR_DEgenesTable[edgeR_DEgenesTable[,3]<0.05,] 



Exercises IV
~~~~~~~~~~~~

# 1. Compare the lists of genes that are called differentially expressed between the two 
# groups by the two methods (DESeq, edgeR). 
# How many genes are DE with each method? 
# How many genes are DE with both methods (overlap)? 
# Hint: Check where the gene name information is stored in the DESeq and in the edgeR tables 
# with significant DE genes. Create a vector called "signDEboth" that contains the gene 
# names ("gene") of the genes that are DE with both methods.

# 2. Draw a venn diagramm with the overlap. 
# Hint: To do this, you first need to create a list from two vectors. The first vector 
# contains # the gene names of the DESeq significant genes and the second vector the gene 
# names of the  edgR significant DE genes. Check the help for ?list on how to do this. 
# Then you need to load the library gplots. The function venn() of that library let's 
# you draw a venn diagramm.


Gene functions
~~~~~~~~~~~~~~
  
# To draw some biological insights it is now useful to find out more about the function of 
# the differentially expressed genes. Most genes were DE with both methods. So let's take 
# these genes and analyze their functions. We can do this using the library topGO. This 
# library allows us to test for enriched Gene Ontology (GO) groups among our genes of interest.
# We also need another library, org.Hs.eg.db, that allows us to assign GO functions to 
# each gene.

library(topGO)
library("org.Mm.eg.db")

# We are going to perform a Fishers exact test, which tests for each GO group whether we 
# have more genes of that group in our list of interesting genes compared to what we would 
# expect from the whole set of expressed genes. To do so, we will create a "topGOdata" 
#object, which needs to have the data in a particular format: We need a vector of all 
# gene names and we need a vector consisting of the factors 0 and 1, indicating whether 
# a gene is among our vector of DE genes.

allGenes=rownames(counts)
geneList=factor(as.integer(allGenes%in% rownames(signedgeR_DEgenes)))
head(geneList)
names(geneList)=allGenes
head(geneList)

# Now we can create our topGodata object, which we call GOdata. This may take a few minutes:

GOdata=new("topGOdata", description="DE-Genes", 
ontology="BP", allGenes=geneList, geneSel=signDEboth, 
annot=annFUN.org, mapping="org.Mm.eg.db", ID="symbol")

# Have a look at the GOdata object: 

GOdata

# It tells us again which ontology we chose, how many genes we had in total (#12050), 
# and how many significant DE genes we had (#1192). Note, that not for all genes GO 
# information is available (# of feasible genes). Finally it tells us about the 
# sturcture of the graph that is avaliable for our analysis.

# We can now run the Fisher exact test for each GO group. The computer might take a 
# few seconds for doing this:

resultFisher=runTest(GOdata, algorithm="classic", statistic="fisher")

# Have a look at the resulting object:

resultFisher

# We see that 860 GO groups are significant with p<0.01.

# Now let's collect the significant terms and then create a nice table with the significant 
# GO groups. Again, this might take a few seconds:
  
signTerms=resultFisher@geneData["SigTerms"]
signGroups=GenTable(GOdata, classicFisher=resultFisher, topNodes=signTerms)
head(signGroups)

# We need to correct our p-values for multiple testing. Here we use the Benjamini-Hochberg 
# correction:
  
qval=p.adjust(signGroups$classicFisher, method="BH")
signGroups=cbind(signGroups,qval)

# How many significant groups do we have?

dim(signGroups[signGroups[,"qval"]<0.05,])

# Ok, let's see them:

signGroups[signGroups[,"qval"]<0.05,][1:20,]

# Which biological functions do you see? What could this mean?

# Lets try now main method of topGO

resultFisher=runTest(GOdata, algorithm="elim", statistic="fisher")
resultFisher
# We see that 266 GO groups are significant with p<0.01.
signTerms=resultFisher@geneData["SigTerms"]
signGroups=GenTable(GOdata, classicFisher=resultFisher, topNodes=signTerms)
head(signGroups)
qval=p.adjust(signGroups$classicFisher, method="BH")
signGroups=cbind(signGroups,qval)
dim(signGroups[signGroups[,"qval"]<0.05,])
# Ok, let's see them:
signGroups[signGroups[,"qval"]<0.05,][1:20,]

# Which biological functions do you see? What could this mean?
