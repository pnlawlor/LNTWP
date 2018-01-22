## Linear-Nonlinear-Time-Warp-Poisson models of neural activity

This repository contains code and sample data related to the following manuscript:

Lawlor P.N., Perich, M.G., Miller, L.E., & Kording, K.P. *Linear-Nonlinear-Time-Warp-Poisson models of neural activity*

which is currently under review but is on Bioarxiv here:

TODO

If you use the code/model, please cite the above paper. If you use the data (to be posted on CRCNS), please cite the above paper and TODO. If you have questions or comments, please email patrick *hyphen* lawlor *at* northwestern *dot* edu

***
### Overview

This project presents a model of trial-to-trial temporal variability in neural data. The intuition comes from movement planning: even when two reaching movements appear to be similar, they could have been planned differently. We can plan movements and execute them immediately, or execute them at a later time. Neural activity related to that planning should therefore not be time-locked to the movement execution. More generally, neural activity and events in the physical world need not be tightly locked. And furthermore, this temporal variability should be similar across neurons performing the same task (e.g., planning).

To build such a model, we combine the Linear-Nonlinear-Poisson (LNP) process and Dynamic Time Warping (DTW). We fit it to data from macaque premotor cortex and find that the effect is considerable.

***
### Contents

This repository contains all of the code used in the manuscript as well as a sample data file to demonstrate how the code works. The main pieces of the code are:


Scripts:
* *main.m* - The main script where the user sets the parameters for the analysis and chooses the data. This script calls all of the other necessary scripts. It also contains code to visualize and analyze the results of the model.
* *load_data.m* - This script loads the data files.
* *initialize.m* - This script initializes some variables and selects the neurons to be analyzed.
* *process.m* - This script generates the model covariates from the loaded data.
* *fit_model.m* - This script actually fits the model.
* *organize_results.m* - This script collects the results of different model fits and plots the summary stats.

The repository also contains directories which contain even more scripts and functions. Most of these are unlikely to be very useful to the user, but feel free to explore.

Requirements:
* *Matlab* - I used v2015b to write the code.
* *Matlab Statistics and Machine Learning toolbox* - for use of *lassoglm* which is used in the model fitting.
* A reasonably recent computer to run the code. The basic model takes about 30 minutes to fit on my machine (64 bit i7 2.3 GHz, 20 GB RAM). I've parallelized the code where possible; most processors can do 2 or 4 virtual cores which helps a bit. It can run considerably faster on a cluster with many cores.

***
### Getting started

1. Download the code and sample data

2. For the simplest possible use, run *main.m* without changing any of the inputs. This reproduces some of the figures from the manuscript (Figs. 7-9).

3. Tinker with the user inputs!
