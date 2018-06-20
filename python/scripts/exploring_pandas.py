#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
#
# @uthor:      acph - dragopoot@gmail.com
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#
# This file contains code that helps to explore pandas module
#
# WARNING: This file is not intend to be executed as a script, so if you execute
# it from a command line, it may not function as you expect.
#------------------------------
# Prepare ipython for interactive mode
%matplotlib qt

# first import the modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


##########################
# 1. LOADING A DATAFRAME #
##########################

# Pandas can load plain text files in many formats. Index_col argument
# tell pandas that use the first column as index.
data = pd.read_table('data/data2.txt', index_col=0, parse_dates=True)

# The DataFrame is a data structure with labeled columns and rows(index)
# to easily access the data.
data.columns
data.index

# A column in a DataFrame is a Series
data['craspi01']
print(type(data['craspi01']))

# Also has methods to conveniently plot data
data.plot()
data.plot.box()
data.plot(subplots=True)

# ... and to describe the data.
c_stats = data.describe()
# some operations returns a new dataframe
c_stats.T['mean']

# We can also do operations on the rows by transposing the matrix
r_stats = data.T.describe()

# .loc and .iloc methods allows the access to rows
c_stats.loc['mean']
c_stats.iloc[1]

# We can use the series with other functions
means = r_stats.loc['mean']
stds = r_stats.loc['std']
xs = range(len(means))
plt.errorbar(xs, means, yerr=stds, fmt='ob--', alpha=0.5)
# ... or maybe prettier
plt.plot(xs, means)
plt.fill_between(xs, means - stds, means + stds, alpha=0.3)

# But usually there is a way to do it on pandas
r_stats.T.plot(y='mean', yerr='std', fmt='bo--', alpha=0.5)
# r_stats.T.plot(y='mean', yerr='std', fmt='bx--', alpha=0.5)
# r_stats.T.plot(y=['mean', 'max'], yerr='std', fmt='bo--', alpha=0.5)

# We can visualize data in many ways
plt.imshow(data.T, cmap='gray')
plt.colorbar(label='Temperature', orientation='horizontal')

# We can obtain smart slices
data.loc['2017-01-06 12:50':'2017-01-06 13:50']
data.loc['2017-01-06 12:50':'2017-01-06 13:50'].plot()

# Finally, what about data normalization?
# Now, we will plot the data relative to the initial value
data.iloc[0]

n_data = data / data.iloc[0]
n_data.plot()

plt.imshow(n_data.T, cmap='gray')
plt.colorbar(label='Temp (normed)', orientation='horizontal')


###############################
# 2. A MORE COMPLEX DATAFRAME #
###############################

# KEGG is a public and specialized database for metabolic data.
kegg = pd.read_table('data/data3.txt', index_col=0)

# Each column (or data Series) may have a different data type
kegg.dtypes

# How many missing data do we have?
# First create a mask of boolean values
miss = kegg.isnull()
# ... then, sum the values. Remember that:
# True == 1
# False == 0
miss.sum()

# What organisms has these null values in nucleotide column?
mask = kegg['nucleotides'].isnull()
nuc_null = kegg[mask]['tax1']
nuc_null.value_counts()
# in proportion?
nuc_null.value_counts() / kegg['tax1'].value_counts() * 100
# well... was that unexpected?

# Although we have a lot of missing data, we can still work
# with it.

# Now, explore the taxa included in the data.
# First, how many prokaryotes and eukaryotes?
tax1 = kegg['tax1'].value_counts()
tax1.plot.pie(autopct='%.1f')

# Second, what eukaryotes taxa do we have? and how many
# organism by taxon?

euk_mask = kegg['tax1'] == 'Eukaryotes'
euk = kegg[euk_mask]

euk_taxa = euk['tax2'].value_counts()
euk_taxa.plot.pie(autopct='%.1f')

# ... and prokaryotes?
prok_mask = kegg['tax1'] == 'Prokaryotes'
prok = kegg[prok_mask]

prok_taxa = prok['tax2'].value_counts()
prok_taxa.plot.pie(autopct='%.1f')

# It is like we can go further to tax3 in Prokaryotes
# ... maybe later ;)

################################
# 3. RELATIONS BETWEEN COLUMNS #
################################

# Let's think... What would be the relation between the
# genome size in nucleotides and the number of ORFs?
plt.figure()
plt.subplot(131)
kegg.plot.scatter(x='nucleotides', y='orfs')

# Well, some points are in a little awkward position
# Remember the euk and prok dataframes?
plt.subplot(132)
euk.plot.scatter(x='nucleotides', y='orfs')
plt.subplot(133)
prok.plot.scatter(x='nucleotides', y='orfs')

# Fairly good...
# Let's see a correlation. .corr() method performs pairwise
# correlations
euk.corr()

prok.corr()

# Look!!! something interesting with enzymes in prokaryotes
# let's explore...
kegg.plot.scatter(x='orfs', y='enzymes')

# Mmm!!! Still awkward points.
# However, there is something fishy.
# We can create another column, that contains the
# proportion of orfs described as enzyme

kegg['proportion'] = kegg['enzymes'] / kegg['orfs']
kegg.plot.scatter(x='orfs', y='proportion')
kegg.corr()                     # correlations?

# Nice plot, and nice negative correlation!!! 

# Wait ... Who the hell has more than 50000 ORFs?
alot_mask = kegg['orfs'] >= 50000
kegg[alot_mask][['tax3', 'spname', 'orfs','proportion']]

# :P well ..... a nice plot always is well received

# we already have an Eukaryotes mask
# remember that we modified the original kegg DataFrame
euk = kegg[euk_mask]

# Prokaryotes are in reality two natural groups
# Bacteria and Archaea
bmask = kegg['tax2'] == 'Bacteria'
amask = kegg['tax2'] == 'Archaea'

bact = kegg[bmask]
arch = kegg[amask]

# We will create a custom plot with matplotlib

plt.figure()
plt.scatter(bact['orfs'], bact['proportion'],
            alpha=0.3, label='Bacteria', s=20)
plt.scatter(arch['orfs'], arch['proportion'],
            alpha=0.3, label='Archaea', s=20)
plt.scatter(euk['orfs'], euk['proportion'],
            alpha=0.3, label='Eukarya', s=20)

leg = plt.legend()
plt.xlabel('ORFs', fontweight='bold')
plt.ylabel('Enzyme proportion', fontweight='bold')

plt.savefig('ecVSorf.png', dpi=300)

######################
# ####               #
# I hope you had fun #
# ####               #
######################

# What do you think about using seaborn?
import seaborn as sns
...
