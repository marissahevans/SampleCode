%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 2 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 2

% Date intialized: Sept. 27th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Polynomial regression.
% Load the file regress1.mat into your MATLAB environment.
load('regress1.mat')

%%
%* Plot variable Y as a function of X. 

figure
plot(x,y, 'o')
axis equal
axis square
xline(0);
yline(0);
title('Regression Data')
legend('data')
box off
set(gca, 'TickDir', 'out')
xlabel('X')
ylabel('Y')

%%
%* Find a least-squares fit of the data with polynomials of
% order 0 (a constant), 1 (a line, parameterized by intercept and slope), 2, 3, 4, and 5. 
%[Compute this using svd and basic linear algebra manipulations that you've learned in class!] 

%Poly 0
[beta0, error0] = betaMathTools(x,y,0);
errorMat(1,1) = error0;
errorMat(2,1) = 0;

%Poly 1
[beta1, error1] = betaMathTools(x,y,1);
errorMat(1,2) = error1;
errorMat(2,2) = 1;

%Poly 2
[beta2, error2] = betaMathTools(x,y,2);
errorMat(1,3) = error2;
errorMat(2,3) = 2;

%Poly 3
[beta3, error3] = betaMathTools(x,y,3);

errorMat(1,4) = error3;
errorMat(2,4) = 3;

%Poly 4
[beta4, error4] = betaMathTools(x,y,4);
errorMat(1,5) = error4;
errorMat(2,5) = 4;

%Poly 5
[beta5, error5] = betaMathTools(x,y,5);
errorMat(1,6) = error5;
errorMat(2,6) = 5;
    
%%
%* On a separate graph, plot the squared error as a function of the order of the polynomial. Which
% fit do you think is ?best?? Explain.

figure
plot(errorMat(2,:), errorMat(1,:), 'LineWidth', 2)
title(['Plot of squared error as a function of polynomial order'])
legend('Error')
xlabel('# Polynomials');
ylabel('Squared Error');
set(gca, 'XLim', [-.5 5.5])
set(gca, 'XTick', [0 1 2 3 4 5])
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 14)
box off

%The fit at 3 polynomials is the best as it grants the largest reduction of
%error. Going further less and less error is reduced for each dimension
%added. 
