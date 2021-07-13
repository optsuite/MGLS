# MGLS
MGLS is a package in MATLAB for multigrid/mutilevel optimization.


MGLS solves an infinite dimensional minimization problem
   min F(X),
using a multigrid algorithm.

The purpose of this code "MGLS" is to demonstrate the algorithms proposed in 

       @article{WenGoldfarb2007b,
               AUTHOR = {Wen, Zaiwen and Goldfarb, Donald},
               TITLE = {A line search multigrid method for large-scale nonlinear optimization},
               JOURNAL = {SIAM J. Optim.},
               FJOURNAL = {SIAM Journal on Optimization},
               VOLUME = {20},
               YEAR = {2009},
               NUMBER = {3},
               PAGES = {1478--1503},
       }


Please refer to the example driver file "TestBatchMGDriver.m" and "GetProb.m" if the
 document is not clear.


# Installation 
   - download MGLS. This should create a directory with
       /MGLS, /MGLS/MG-OPT, /MGLS/MG-Lineq, /MGLS/Model, /MGLS/utils 
   - Add these folders to your MATLAB path

# Quick Start:
 `>> TestBatchMGDriver`


# The Test Problems
 A collection of test problems is stored in "Model"

# Acknowledgement
 
Some parts of the linear multigrid methods (MG-Lineq) for solving Ax = b are modified from the package MGLab by James Bordner and Faisal Saied 

# The Authors
We hope that MGLS is useful for your application.  If you have any bug reports or comments, please feel free to  the toolbox author:

* Zaiwen Wen, wenzw at pku.edu.cn

# Copyright
------------------------------------------------------------------------
  Copyright (C) 2017, Jiang Hu, Zaiwen Wen

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without  even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>