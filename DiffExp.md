# Wragling Count matrix

> Noysi raw data: increase signal-to-noise-ratio

Infovis:

1. Library size ie. sum of total reads in the sample
2. Frequency of reads counts
3. Most genes are expressed in low/high amount, many outliers
4. Counts on log2-scale, median differs among samples
5. Percent of genes vs average of counts to figure out N% of genes than have average expression less than 1 (or other 

>  We can remove these percent genes having very low reads counts



Plot expression levels before, log2(normalized) and quantil normalization.



# Infovis Count data

- Check correlation bewtween samples
  - Scatterplots can be used
- Hierarchical clustering 
  - A simple euclidian distance based on hierarchical clustering should reveals trends in the data.
  - Correlation between groups (plot all pairs between samples)
- Perform principal component analysis (PCA)
  - Due to count matrix is a Multidimensinal dataset

### PCA (variation between components)

- The samples can be clustered in two dimension using Multi-dimensional scaling (MDS) plots.
- The distance between each pair of samples can be interpreted as the leading log-fold change between the samples for the genes that best distinguish that pair of samples.
- The difference between groups are much larger than those within groups, meaning that there are likely to be statistically significant differences between the groups.



# Differential expression analysis

For each gene in each sample we have a measure of abundance. ie. Number of reads mapping across gene. We wanto to know wheter there is a statistically significant difference in abundance of genes between treatments / groups / genotypes. **identifying genes ...

Transcriptonal outputs has a direct consequence on the phenotype. Studying expression differences allows us to identify genes responsible for phenotypic differences.

## Statisticall significants

### Hypothesis testing

> H testing for testing differential expression is posed as follow:

- Ho : The distribution of expression values in treatment and control samples do not differ
- Ha : It differs between treatment and control

> Statistical test report the p-value and logFC:

- The  Pvalue is the probability of obtain a resut equal to or more extreme than what was actually obsered (under null hypothesis) .

  - > Use edgeR to measure Pvalues of genes.

- The null hypothesis is rejected if the pvalue is less than or equal to a samll fixed but arbitrary pre-defined threshold value a, which is referred to as the level of significance.

- logFC is the mean difference between log2 transformed expression values between treatment and control.

- Positive FC reflects up-regulation of a gene (ex. condition1 in contrast to control), whereas negative score reflects down-regulation.

## DEseq and edgeR for DiffExp

## 

| DESeq                                                        | EdgeR                                                        |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| Analyzes the median of the ratios across all samples         | Picks one sample as the reference. Analyzes the weighted mean of the ratios of each sample to that reference. |
| If median is not 1 the read counts for that sample are corrected so that it become 1 | If mean is not 1 the libsize for that sample is corrected so that it become 1 |
| Normalized read counts are calculated from raw read counts divided by the correction factor. | Normalized read counts are recalculated with the corrected libsize |
| `estimateSizeFactors()`                                      | `calcNormFactors()`                                          |
|                                                              |                                                              |

### Adjust p-value to reflect multiple testing

Most common statistical teste were good in the times before large scale datasets. Ex. P = 0.01 means of you were repeating the same test 100 times you would expect the same outcome only ones ( 1 in 100). ie. You would not expect the same outcome even if you repeated the test 100 times

We are testing thousands of genes at the same time . If you had 10,000 genes and none were different, you would expect 100 to have p < 0.01 by chance alone.

**it is necessary to adjust pValues to reflect multiple testing.** Some methods are:

- Bonferroni correction
  - Multiple pvalue by number of tests done
  - very conservative. greatly increase false negative rate
- Benjamin & Hochberg correction:
  - Less conservative, willing to accept some number of false positive in order to decrease false negative rate.
- Calculate False Discovery Rate:
  - Try multiple p-value cutoffs and use what you're what willing to accept.

### Visualizing significat (FDR) differential Expression



- heatmaps 
- volcano plots
- boxplot

## Biology of DEGs

...

### Pathway analysis

- GO
- KEGG
- GSEA
- DEG networks

# Final words

...



Open 

# Reference

Vijaykumar Yogesh Muley (2018). Differential Expression analysis. III Summer School in bioinformatics. UNAM Juriquilla Mexico