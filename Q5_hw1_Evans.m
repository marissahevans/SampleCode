%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 1 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 5

% Date intialized: Sept. 12th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc
%% Gram-Schmidt. 

% A classic method for constructing an orthonormal basis is known as Gram- Schmidt 
% orthogonalization. First, one generates an arbitrary unit vector 
% (e.g., by normalizing a vector created with randn). Each subsequent basis 
% vector is created by generating another arbitrary vector, subtracting off 
% the projections of that vector along each of the previously created basis 
% vectors, and normalizing the remaining vector.

%Basic test of Gram-Schmidt priciples 
vec1 = randn(3,1);
vec2 = randn(3,1);
lenV1 = sqrt(vec1'*vec1);

unit1 = vec1/lenV1;

projectV = (vec2'*unit1)*unit1;
unit2 = (vec2 - projectV)/sqrt(vec2'*vec2);

orthoCheck = unit1'*unit2;

%%
%* Write a MATLAB function gramSchmidt that takes a single argument, N, specifying 
% the dimensionality of the basis. It should then generate an N × N matrix 
% whose columns contain a set of orthogonal normalized unit vectors. Try your function for N = 3, and 
% plot the basis vectors (you can use MATLAB's rotate3d to interactively examine these). 
gramMat = gramSchmidt(3);

a = gramMat(:,1);
b = gramMat(:,2);
c = gramMat(:,3);
d = zeros(3,1);

%3D plot of orthogonal vectors produced with GramSchmidt
figure
quiver3(d,d,d,a,b,c)
xlabel('X')
ylabel('Y')
zlabel('Z')
title('3D Gram-Schmidt ')


%%
%* Check your function numerically by calling it for an N larger than 3 and 
% verifying that the resulting matrix is orthonormal. 

gramMat2 = gramSchmidt(7);
orthoCheck = gramMat2'*gramMat2
disp('The matrix is orthogonal because multiplying by the transpose results in the identity')




