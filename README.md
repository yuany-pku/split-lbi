MATLAB version of split lbi in the paper: 
Chendi Huang, Xinwei Sun, Jiechao Xiong and Yuan Yao.
Split LBI: An Iterative Regularization Path with Structural Sparsity.
Advances In Neural Information Processing Systems. 2016: 3369-3377.
=======================================
Version: 0.1
Date: 26/03/2017
Authors: Chendi Huang, Xinwei Sun, Jiechao Xiong and Yuan Yao
License: GPL-2

This package is an implementation of the split lbi in MATLAB.

Requirements
------------

- MATLAB 7.8 or later
- CVX(cvx version: 2.1)

Installation
------------
The following steps will help you to install the package.
1 Unzip this archive anywhere on your machine, say to '<path>’. Then in 
MATLAB Commond Window, type:
	addpath('<path>/splitlbi_matlab’);

2 If the cvx file is in the ‘<cvx_path>’, then type:
        addpath(‘<cvx_path>/cvx’);
(You can SKIP this step If you already add the directory of cvx in MATLAB search path)

3 Type:
        cvx_setup;
    
To test your installation, you can run one of the demo files in the ‘First examples’
subdirectory in the Command Window.  This file contain examples of linear and 
logit models of split lbi with grouped and ungrouped sparsity type, which is a good 
place to start with.

Documentation
------------
Main codes of splitlbi are in the file of +splitlbi, which includes the main entry for 
split lbi: +splitlbi/splitlbi.m. For an overview of all functions in the +splitlbi, run

    help splitlbi

The file +lbi is Linearized Bregman Iteration(LBI) for l1 sparsity, without variable 
splitting term. For an overview of all functions in the +lbi, run

    help lbi

For details of LBI, please reference: https://arxiv.org/abs/1406.7728

The file ‘exp_simulation’ reproduce the simulation experiments, such as IRR curve, 
AUC and regularization path plot in the paper. (Note that the figures in the 
paper are plotted by R studio). 

The file ‘examples’ contains basketball and university applications in the paper.