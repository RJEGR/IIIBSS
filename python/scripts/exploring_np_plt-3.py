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

########################
# 5. NOW SOME MATRICES #
########################

# How about some noise in 2d?
noise = np.random.rand(100, 100)
# shape of the array
noise.shape

# You can visualize a matrix (2d) with the following functions
plt.matshow(noise)
plt.colorbar()
#
plt.imshow(noise)
# And you can see the histogram if you flatten the matrix
# What disribution is it?
plt.hist(noise.flatten())

# lets create a convenience function to visualize the matrix data


def show_matrix(matrix):
    plt.figure()
    plt.subplot(121)
    plt.imshow(matrix)
    plt.colorbar()
    plt.subplot(122)
    plt.hist(matrix.flatten())


# No plot the uniform noise
show_matrix(noise)
plt.suptitle('Uniform noise')

# Now, normal distributed noise
n_noise = np.random.randn(100, 100)
show_matrix(n_noise)

# You can create a mask for the 2d array. For example,
# we are going to replace every value < 0 to 0
mask = n_noise < 0
# Check mask shape
mask.shape == n_noise.shape
# Visualizing the mask
plt.imshow(mask)

# replace values
n_noise[mask] = 0
show_matrix(n_noise)

# ---
# -- Now, we are going to use a matrix stored in a file
data = np.loadtxt('data/data1.txt', dtype=np.int)
data.shape

# ... :D
show_matrix(data)

# Use a more suitable argument selection for the image
plt.imshow(data, cmap='gray', vmin=0, vmax=255)

# ...
# The data represent a picture, that is cool.
# We can take an slice of the picture
fragment = data[100:850, 750:1600]
plt.imshow(fragment, cmap='gray', vmin=0, vmax=255)

# The picture is a little dark, we can create new pictures
# with more brightness.

pic1 = data * 2
pic2 = data + 100 # increase the luminosity of the pixels

plt.subplot(121)
plt.imshow(pic1, cmap='gray', vmin=0, vmax=255)
plt.subplot(122)
plt.imshow(pic2, cmap='gray', vmin=0, vmax=255)

# --- What about a nice figure with the results and including
# the histograms?

# Figure with 3 columns and 2 rows
fig = plt.figure()
# original image
plt.subplot(231)
plt.imshow(data, cmap='gray', vmin=0, vmax=255)
plt.axis('off')
# histogram
plt.subplot(234)
plt.hist(data.flatten())
plt.xlim((0, 255))

# multiplied
plt.subplot(232)
plt.imshow(pic1, cmap='gray', vmin=0, vmax=255)
plt.axis('off')
# histogram
plt.subplot(235)
plt.hist(pic1.flatten())
plt.xlim((0, 255))

# + 100 brightness
plt.subplot(233)
plt.imshow(pic2, cmap='gray', vmin=0, vmax=255)
plt.axis('off')
# histogram
plt.subplot(236)
plt.hist(pic2.flatten())
plt.xlim((0, 255))

plt.tight_layout()
