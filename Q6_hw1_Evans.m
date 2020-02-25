%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 1 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 6

% Date intialized: Sept. 12th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Null and Range spaces
% Load the file mtxExamples.mat into your MATLAB world.You?ll find a set of 
% matrices named mtxN, with N = 1, 2.... For each matrix, use the SVD to: 
load('mtxExamples.mat');

%%
% (a) determine if there are non-trivial (i.e., non-zero) vectors in the 
% input space that the matrix maps to zero (i.e., determine if there's a nullspace). 
% If so, write a MATLAB expression that generates a random example of such a vector, 
% and verify that the matrix maps it to the zero vector; 

%%
%* MTX1 Null Space Test
[U1 S1 V1] = svd(mtx1);
disp('This matrix has no null space')

%%
%* MTX2 Null Space Test
[U2 S2 V2] = svd(mtx2)
disp('This matrix has 2x rows of null space')
testVec2 = randn(3,1);
testVec2(1) = 0;

inputVec2 = V2*testVec2; %get the input vector by multiplying the test vector by V2 (inverse of V2 transpose)
output2 = U2*S2*V2'*inputVec2; %output is 0 as all non-zero components map to zero

nulltest2 = mtx2*inputVec2 
disp('test was successful as the matrix maps to zero when multiplied by the input vector')

%% 
%* MTX3 Null Space Test
[U3 S3 V3] = svd(mtx3);
disp('This matrix has 1x row of null space')
testVec3 = randn(3,1);
testVec3(1:2) = 0;

inputVec3 = V3*testVec3; %get the input vector by multiplying the test vector by V3 (inverse of V3 transpose) 
output3 = U3*S3*V3'*inputVec3;

nulltest3 = mtx3*inputVec3
disp('test was successful as the matrix maps to zero when multiplied by the input vector')

%%
%* MTX4 Null Space Test
[U4 S4 V4] = svd(mtx4);
disp('This matrix has 1x row of null space')
testVec4 = randn(2,1);
testVec4(1) = 0;

inputVec4 = V4*testVec4; %Input vector found by multiplying the test vector by V4 (inverse of V4 transpose)
output4 = U4*S4*V4'*inputVec4;

nulltest4 = mtx4*inputVec4
disp('test was successful as the matrix maps to zero when multiplied by the input vector')


%%
% (b) write a MATLAB expression to generate a random vector y that lies in the 
% range space of the matrix, and then verify that it's in the range space by 
% finding an input vector, x, such that Mx = y.

%%
%* MTX1 Range Space Test
%inverse of MTX1
S1I = S1;
S1I(1,1) = 1/S1I(1,1); %change the values of S to their inverse
S1I(2,2) = 1/S1I(2,2);
M1Inv = V1*S1I'*U1'; %calculate M# using the pseudo inverse of S

%Compare Mx to y
y1 = U1(:,1)*randn(1); %pull the range space values from U and raise by a random scaler 
x1 = M1Inv*y1; %multiply by M# to get the value for x
rangetest1 = mtx1*x1; %test for equality 

disp(['The value of y1 (' num2str(y1') '), which is derived from the range space'])
disp(['is equal to the value of rangetest1 (' num2str(rangetest1') ')'])
disp('which is the product of mtx1 and x1. Input vector x1 is found')
disp('by multiplying the inverse of mtx1 by y1.')


%%
%* MTX2 Range Space Test

% Calculate inverse of MTX2
S2I = S2;
S2I(1) = 1/S2I(1);
M2Inv = V2*S2I'*U2'; %Calculate M#

%Compare Mx to y
y2 = U2(:,1)*randn(1); %get y by taking the range space values from U and applying a random scaler 
x2 = M2Inv*y2;
rangetest2 = mtx2*x2; %test for equality 

disp(['The value of y2 (' num2str(y2') '), which is derived from the range space'])
disp(['is equal to the value of rangetest2 (' num2str(rangetest2') ')'])
disp('which is the product of mtx2 and x2. Input vector x2 is found')
disp('by multiplying the inverse of mtx2 by y2.')

%%
%* MTX3 Range Space Test

% Calculate inverse of MTX3
S3I = S3;
S3I(1,1) = 1/S3I(1,1);
S3I(2,2) = 1/S3I(2,2);
M3Inv = V3*S3I'*U3'; %find M#

%Compare Mx to y
y3 = U3(:,1)*randn(1);
x3 = M3Inv*y3;
rangetest3 = mtx3*x3;

disp(['The value of y3 (' num2str(y3') '), which is derived from the range space'])
disp(['is equal to the value of rangetest3 (' num2str(rangetest3') ')'])
disp('which is the product of mtx3 and x3. Input vector x3 is found')
disp('by multiplying the inverse of mtx3 by 3.')

%%
%* MTX4 Range Space Test

% Calculate inverse of MTX4
S4I = S4;
S4I(1,1) = 1/S4I(1,1);
M4Inv = V4*S4I'*U4'; %find M#

%Compare Mx to y
y4 = U4(:,1)*randn(1);
x4 = M4Inv*y4;
rangetest4 = mtx4*x4;

disp(['The value of y4 (' num2str(y4') '), which is derived from the range space'])
disp(['is equal to the value of rangetest4 (' num2str(rangetest4') ')'])
disp('which is the product of mtx4 and x4. Input vector x4 is found')
disp('by multiplying the inverse of mtx4 by y4.')

