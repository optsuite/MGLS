% MGLS: A MATLAB Demo for multigrid/mutilevel optimization
%
% -------------------------------------------------------------------------
%  This package is still under development
% -------------------------------------------------------------------------
%
% 1. Introduction
%
% Thank you for downloading MGLS! MGLS is a package for
% solving an infinite dimensional minimization problem
%   min F(X),
% using a multigrid algorithm
%
% The purpose of this code "MGLS" is to demonstrate the algorithms proposed in 
%
%       @TECHREPORT{WenGoldfarb2007b,
%           author = {Wen, Zaiwen and Goldfarb, Donald},
%           title = {Line search Multigrid Methods for Large-Scale NonConvex Optimization},
%           institution = {Dept of IEOR, Columbia University},
%           year = {2007}
%       }
%
%       @TECHREPORT{WenGoldfarb2007a,
%           author = {Wen, Zaiwen and Goldfarb, Donald},
%           title = {A Line search Multigrid Method for Large-Scale Convex Optimization},
%           institution = {Dept of IEOR, Columbia University},
%           year = {2007}
%        }
% 
%
% *. Please refer to the example driver file "TestBatchMGDriver.m" and "GetProb.m" if the
% document is not clear.
%
%           
% -------------------------------------------------------------------------
% 2. Installation and Setup
% 
%   3.1) Remove any old version of MGLS
%   3.2) unzip MGLS.zip. This should create the structure
%       /MGLS, /MGLS/MG-OPT, /MGLS/MG-Lineq, /MGLS/Model, /MGLS/utils 
%   3.3) Add these folders to your MATLAB path
%
% -------------------------------------------------------------------------
% 3. Quick Start:
% 
% >> TestBatchMGDriver
% 
%
% -------------------------------------------------------------------------
% 4. The Test Problems
% A collection of test problems is stored in "Model"
%
%
% -------------------------------------------------------------------------
% 5. Thanks
% 
% Some parts of the linear multigrid methods (MG-Lineq) for solving Ax = b are
% modified from the package MGLab by James Bordner and Faisal Saied 
%
% 6. License
% 
% See License.txt for the license of this program
%
% -------------------------------------------------------------------------
% 7. The Authors
%
% We hope that MGLS is useful for your application.  If you have
% any bug reports or comments, please feel free to email one of the
% toolbox authors:
%
%   Zaiwen Wen, zw2109@columbia.edu
%
% Enjoy!
% Zaiwen
%
%

