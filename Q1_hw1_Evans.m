%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 1 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 1

% Date intialized: Sept. 12th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Question 1

% Testing for (non)linearity. Suppose, for each of the systems below, you observe the indicated
% input/output pairs of vectors (or scalars). Determine whether each system could possibly
% be a linear system. If so, provide an example of a matrix that is consistent with the observed
% input/output pairs, and state whether you think that matrix is unique (i.e., the only matrix
% that is consistent with the observations). If not, explain why.


%%
%*System 1: x = 1 -> [4, 6]
%          y = 2.5 -> [10, 14]

%This is not a linear system because it does not show additivity or homogenity. The
%scaled output of the first option (multiplied by 2.5) would be [10,15] 
%which does not equal the given output of the second vector. 


%%
%* System 2: x = [6, 3] -> [12, 12]
%           y = [-2, -1] -> [-6, -6]

% This system is not linear because when the two inputs are scaled to be
% the same (second function multiplied by -3), the outputs do not match 
% y = 6, 3 -> 18,18 (does not match x system [12,12])

%%
%* System 3: x = [1, 2] -> [5, -1]
%           y = [1, -1] -> [1, 4]
%           z = [3, 0] -> [7, 8]

% This system is not linear becuse when the y vector is mutiplied by a
% constant 2 to make x+y=z for the input, the output does not match when 
% multiplied by the same constant 7,7 ~= 7,8

%%
%* System 4: x = [2, 4] -> 0
%           y = [-2, 1] -> 3

x = [2;4]; %input 1
y = [-2;1]; %input 2
testMat = [-1.2;.6]; %solved for the system by setting x and y both to zero
output1 = x'*testMat
output2 = y'*testMat

disp('the output of the matrix testMat (-1.2;.6) is equal to the system output and is stable with regards to homogeneity and additivity')

%This is a linear system
%%
%* System 5: x = 0 -> [1, 2]

%This is not a linear system, any matrix multiplied by zero will be a zero
%matrix, thus a non-zero intiger output is not possible. 
