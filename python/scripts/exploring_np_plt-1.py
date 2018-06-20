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
# iPython contains a lot of convenient functions called magic functions
%magic
# Magic functions only work in ipython and Jupyter notebooks. This means
# that the previous code will raise an error if executed in python
# interpreter or in a script.


# first import the modules
import numpy as np
import matplotlib.pyplot as plt


#####################
# 1. A BASIC PLOT   #
#####################

# You can create a numpy array from an existing data type like a list
mylist = [1, 2, 3, 4, 5]
x = np.array(mylist)

# Also, the functions included in numpy module returns a numpy array.
y = np.linspace(1, 5, 5)

# You can plot these data easily using matplotlib (plt)
plt.plot(x, y)
# In a script context, matplotlib does not display the plot by default.
# If you want to display the plot you need to use the show function
# as follows :
# plt.show()

# ... to store the plot
plt.savefig('plot1.tif', dpi=100)

# Close figure - This is necessary to avoid undesired mixed plots
# or to avoid RAM memory saturation.
plt.close()


#########################################
# 2. MORE PLOTS - OPERATIONS WITH NUMPY #
#    ARRAYS                             #
#########################################

# A numpy array is a Python object and is implemented to easily
# apply many operations to it.

x = np.linspace(0, 50, 100)
z = x[::-1]
y = x ** 2
y2 = x / x.max()
y3 = 2 ** x

# Now the plots

plt.plot(x, z)
# plt.show()
# plt.close()

plt.plot(x, y)  # power function
# plt.show()
# plt.close()

plt.plot(x, y2)  # normalization or scaling

plt.plot(x, y3)  # exponential function

# Really straightforward ... what happen if...
plt.plot(x, x * z)

# Now lets make a more complete plot

plt.plot(x, y, '--b', label='$y = x^2$')
plt.xlabel('x')
plt.ylabel('y')
leg = plt.legend()

plt.title('First plot')
plt.tight_layout()

#plt.savefig('first_plot.png', dpi=200)


##############################
# 3. NUMPY RANDOM SUBPACKAGE #
#             AND            #
#    MORE MATPLOTLIB FIGURES #
##############################

# numpy includes some subpackages like random data generation
# np.random

# Create some random normal distributed data:
n_random = np.random.randn(1000)

# Plot as a histogram
plt.hist(n_random)

# To few bars? you can change the bins parameters using a
# list of values or a scalar.
plt.hist(n_random, bins=25)

# Plot as a boxplot
plt.boxplot(n_random)

# You can create also another types of random data, for
# example: uniform distributed data or gumbel distributed data
u_random = np.random.uniform(size=1000)
g_random = np.random.gumbel(size=1000)

# mmm, to much data... and to much plots...
# maybe a simple boxplot may show us all the data together
plt.boxplot([n_random, u_random, g_random])

# it seems that the plot lacks the name corresponding of each
# boxplot. We can assign it using plt.xticks function
plt.xticks([1, 2, 3], ['Normal', 'Uniform', 'Gumbel'])
plt.xlabel('Distribution')

# It Looks nice, still, it might be confusing, particularly in the
# uniform distribution. Maybe in this case the boxplot is not properly
# used. So we may add all the histogram in one figure.

# Let's create a compound figure
# histograms
plt.figure()


## -- BEGIN Advanced -- ##
# put all the plots in one figure!!

# Figure
plt.figure(figsize=(7, 7)) # hacemos un grid de 7*7
# empezamos a acomodar cada uno de los graficos en el grid creado. describiendo:
# numero de filas, numero de columnas, posicion donde se imprime el plot 
plt.subplot(211)
plt.boxplot([n_random, u_random, g_random]) # metemos la figura en las coordenadas del grid
# y editamos la primer figura
plt.xticks([1, 2, 3], ['Normal', 'Uniform', 'Gumbel'],
           fontsize='large')
plt.title('Distributions')

# Histograms
# normal
plt.subplot(234)
plt.hist(n_random, bins=40)
plt.ylabel('Count')

# uniform
plt.subplot(235)
plt.hist(u_random, bins=40)

# gumbel
plt.subplot(236)
plt.hist(g_random, bins=40)
plt.tight_layout()

## --END Advanced -- ##
