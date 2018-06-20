#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#------------------------------
#
# @uthor:      acph - dragopoot@gmail.com
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
#
# This file contains code that helps to explore numpy and matplotlib modules
#
# WARNING: This file is not intend to be executed as a script, so if you execute
# it from a command line, it may not function as you expect.
#------------------------------

# Prepare ipython for interactive mode
%matplotlib qt

# first import the modules
import numpy as np
import matplotlib.pyplot as plt


####################################
# 4. MORE PLOTTING AND A LITTLE OF #
#    SCIPY STATS                   #
####################################

from scipy import stats

# -----
# ----- What about a simple correlation example?
# -----
# scipy stats submodule has a lot of commonly used
# statistics functions

# lets create some data
x = np.arange(26)

pcor = stats.pearsonr(x, x)
print('r = ', pcor[0])
print('p-value = ', pcor[1])

# now, some noise to the data
r_factor = random.normal(loc=0, scale=5, size=len(x))
y = x + r_factor

# - Pearson correlation
pcor = stats.pearsonr(x, y)

# Least squares with numpy polyfit
fit = np.polyfit(x, y, 1)
m, b = fit[0], fit[1]


def lin_func(x):
    return m*x + b

# To plot the fitted line we need to create the corrsponding
# x,y points
x2 = np.linspace(0, 25, 1000)   # a lot of points (1000)
y2 = lin_func(x2)               # Corresponding y values

# Lets plot the data

plt.scatter(x, y, color='g', s=20, label='Data')
plt.plot(x2, y2, color='r', label='Linear fit')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()

annotation = 'Linear fit : $y = {:.2}x + {:.2}$\n $r = {:.2}$\n $p = {:.2}$'
annotation = annotation.format( m, b, pcor[0], pcor[1])
ann = plt.annotate(annotation, xy=(0, 18), fontsize='small')

# ----
# ---- What about a t-test?
# ----
# Lets create 2 independent samples from random normal distributions

# First visualize the distributions that we are going to sample
# Means 18 and 22 (~20% of difference between means)
# Same standard deviation = 2

plt.subplot(211)
plt.hist(random.normal(18, 2, size=10000), bins=20)
plt.title('Distribution 1')
plt.xlim((10, 30))

plt.subplot(212)
plt.hist(random.normal(22, 2, size=10000), bins=20)
plt.title('Distribution 2')
plt.xlim((10, 30))

########################
# Now lets sample data #
########################
n1 = 10
n2 = 10

sample1 = random.normal(18, 2, size=n1)
sample2 = random.normal(22, 2, size=n2)

# Mean value
mean1 = sample1.mean()
mean2 = sample2.mean()

# Standard deviation
std1 = sample1.std()
std2 = sample2.std()

# Standard error of the mean
sem1 = std1 / np.sqrt(n1)
sem2 = std2 / np.sqrt(n2)

# Now, plot the data
# boxplot
plt.subplot(311)
plt.boxplot([sample1, sample2])
plt.title('Boxplot')
plt.ylabel('Arbitrary units')
plt.xticks([1, 2], ['Sample 1', 'Sample 2'])

# point plot with standard deviation of the mean
plt.subplot(312)
plt.errorbar([1, 2], [mean1, mean2], yerr=[
             std1, std2], label='Mean +- std', fmt='o')
plt.title('Point plot (std)')
plt.ylabel('Arbitrary units')
plt.xticks([1, 2], ['Sample 1', 'Sample 2'])
plt.xlim((0.5, 2.5))

# bar plot with standard error of the mean
plt.subplot(313)
plt.bar([1, 2], [mean1, mean2], yerr=[sem1, sem2],
        label='Mean +- s.e.m', width=0.4)
plt.title('Bar plot (s.e.m.)')
plt.ylabel('Arbitrary units')
plt.xticks([1, 2], ['Sample 1', 'Sample 2'])
plt.xlim((0, 3))
plt.ylim((15, 24))

plt.tight_layout()


# Are there statistically significant differences between the measn?

# -- Hypothesis testing --
# What about a simple t-test
# Are the data normally distributed?

def is_normal(x, pvalue=0.05):
    """Returns True if the data is normally
       distributed with a given p-value threshold

    - x : numpy 1d array or list of integers|floats
    - pvalue : threshold to reject null hypothesis
    """
    s, pv = stats.shapiro(x)
    if pv > pvalue:
        return True
    else:
        return False


print('sample1', is_normal(sample1))
print('sample2', is_normal(sample2))


# if both samples are normally distributed you can
# use student's t
tstat = stats.ttest_ind(sample1, sample2)
print('T test p-value:', tstat.pvalue)

# Is there a statistically significant difference between the
# sample means (t >= 0.05)?
# If not, how can you simulate samples that are statistically different?

### -- Exercise

# What would happen to the samples parameters
# (std, s.e.m.) if you change the sample size to:
# 1) n1 = 3, n2 = 3
# 2) n1 = 20, n2 = 20
# Take special attention to the s.e.m. and to the boxplot
# save the plots and take note of the differences.

# Is there a statistically significant difference between the
# sample means?
# If not, how can you simulate samples that are statistically different?


###################################################
###################################################
###################################################
###################################################
###################################################

### -- Another exercise
# t-test brute force simulation

# Try different sample sizes
# maintaining the same mean and stds!!!!
samples = {3: 0,
           10: 0,
           20: 0}

# making 1000 simulations
for i in range(1000):
    # testing each sample sizes
    for j in [3, 10, 20]:
        # samples
        s1 = random.normal(18, 2, size=j)
        s2 = random.normal(22, 2, size=j)
        # testing
        ttest = stats.ttest_ind(s1, s2)
        # Reject Ho?
        if ttest.pvalue <= 0.05:
            samples[j] += 1

samples
###################################################
