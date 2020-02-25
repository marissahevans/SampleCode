%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 4 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 3

% Date intialized: Oct. 31th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Multi-dimensional Gaussians.

%% PART A

%%
%* Write a function samples = ndRandn(mean, cov, num) that generates a set of
% samples drawn from an N-dimensional Gaussian distribution with the specified mean
% (an N-vector) and covariance (an NxN matrix). The parameter num should be optional
% (defaulting to 1) and should specify the number of samples to return. The returned 
% value should be a matrix with num rows each containing a sample of N elements. 

% (Hint: use the MATLAB function randn to generate samples from an Ndimensional Gaussian 
% with zero mean and identity covariance matrix, and then transform these to achieve the 
% desired mean/cov. Recall that the covariance of Y = MX is E(Y Y') = MCXM' where CX 
% is the covariance of X). 


% Please use mean µ = [4, 5] with CX = [9, ?5; ?5, 6] to 
% sample and scatterplot 1,000 points to verify your function work as intended.

mu = [4,5];
CY = [9, -5; -5, 6];
num = 1000;

[samples] = ndRandn(mu, CY, num);

figure
scatter(samples(:,1), samples(:,2))
title(['Random Sample Transformed by mu and CX'])
xlabel('X')
ylabel('Y')
box off
set(gca, 'TickDir', 'out')

%testing my function against MATLAB function:
R = mvnrnd(mu,CY,num);
figure
scatter(samples(:,1), samples(:,2))
hold on 
scatter(R(:,1), R(:,2))
title(['Verification ndRandn matches Matlabs mvnrnd'])
legend('Sample from my function', 'Sample from Matlab function')
xlabel('X')
ylabel('Y')
box off
set(gca, 'TickDir', 'out')

%% PART B

%%
%* Now consider the marginal distribution of a generalized 2-D Gaussian with mean µ
% and covariance ? in which samples are projected onto a unit vector u? to obtain a 1-D
% distribution. Write a mathematical expression for the mean, µ?, and variance, ??2, of this
% marginal distribution as a function of u?.

% 

% mean = u^'*mean(x) || mean equals the unit vector times the mean vector of
% the 2d Gaussian. In matab: mu * unitVecs'

% sigma^2 = u^'*covaiance*u^ || variance equals the unit vector times the covariance 
% times the unit vector. In matlab: diag(unitVecs*CY*unitVecs')


%% 
%* Check it for a set of 48 unit vectors spaced evenly around the unit circle. 

%Create 48 unit vectors evenly around a circle 

n=48;
tet=linspace(-pi,pi,n);
unitVecs(:,1)=cos(tet);
unitVecs(:,2)=sin(tet);

%%
%* For each of these, compare the mean and variance predicted 
% from your mathematical expression to the sample mean and variance estimated
% by projecting your 1,000 samples from part (a) onto u?. 

%Predicted mean and variance 
predMean = mu * unitVecs';
predVar = diag(unitVecs*CY*unitVecs');

%Sample mean and variance
projectedData = samples*unitVecs';
calcMean = mean(projectedData);
calcVar = var(projectedData);

%%
%* Stem plot the mathematically computed mean and the sample mean (on the same plot), 
% and also plot the mathematical variance and the sample variance.
numVecs = 1:48;

figure
subplot(2,1,1)
stem(numVecs,predMean)
hold on
stem(numVecs,calcMean)
title(['Sample Mean vs Mathmatical Mean'])
legend('Mathmatical Mean','Sample Mean')
xlabel('Unit Vectors')
ylabel('Mean')
box off
set(gca, 'TickDir', 'out')

subplot(2,1,2)
stem(numVecs,predVar)
hold on
stem(numVecs,calcVar)
title(['Sample Variance vs Mathmatical Variance'])
legend('Mathmatical Variance','Sample Variance')
xlabel('Unit Vectors')
ylabel('Variance')
box off
set(gca, 'TickDir', 'out')

% The values are almost identical however due to the random sampling there
% is a small about of mismatch. The larger the sample tested the closer the
% two values would be. 
%% PART C

%%
%* Now scatterplot 1,000 new samples of a 2-dimensional Gaussian using µ and CX in
% part (a). Measure the sample mean and covariance of your data points, comparing to
% the values that you requested when calling the function.

mu = [4,5]
CY = [9, -5; -5, 6]
num = 1000;

[samples2] = ndRandn(mu, CY, num);

checkMean = mean(samples2)
checkCov = cov(samples2)

for ii = 1:1000
[samples4] = ndRandn(mu, CY, num);

checkMean = mean(samples4);
checkCov = cov(samples4);
mean1000(ii,:) = mu - checkMean;
cov1000(ii,:) = CY(1,2) - checkCov(1,2);
end 

figure
subplot(1,2,1)
histogram(mean1000)
title(['Sample Mean Difference over 1000x Samples'])
xlabel('Mean Difference')
ylabel('# of Samples')
box off
set(gca, 'TickDir', 'out')

subplot(1,2,2)
histogram(cov1000)
title(['Sample Cov Difference over 1000x Samples'])
xlabel('Cov Difference')
ylabel('# of Samples')
box off
set(gca, 'TickDir', 'out')

%given the random sample the mean and covariance are within the expected range for 1000
%samples, however will not fall specifically on the value entered unless
%the sample approaches infinity. By calculating the difference between the
%sample mean and the requested population mean we can see in the histogram
%that over 100 samples the difference is normally distributed around zero.
%This is also the case when comparing the difference in covariance between
%the samples and the requested population. 

%%
%*  Plot an ellipse on top of the
% scatterplot by generating unit vectors equi-spaced around the circle, and transforming
% them with a matrix as in part (a) to have the same mean and covariance as the data. 

[u,s,v] = svd(CY*2); %multiplied by 2 because unitVecs' starting covairance was .5 instead of 1
M = sqrt(s)*v';
Y = unitVecs*M;

transEllipse = Y+ mu;

figure
scatter(samples2(:,1), samples2(:,2))
hold on
plot(transEllipse(:,1),transEllipse(:,2), 'LineWidth',2);
title(['Random Sample Transformed by mu and CX'])
legend('Transformed Random Data', 'Transformed Unit Circle')
xlabel('X')
ylabel('Y')
box off
set(gca, 'TickDir', 'out')


%%
%* Try this on three additional random data sets with different means and 
% covariance matrices. Does this ellipse capture the shape of the data?
figure

for ii = 1:3
mu3 = randi([-20,20],1,2);
CY3 = randi([-20,20],2,2);
[samples3] = ndRandn(mu3, CY3, num);

[u,s,v] = svd(CY3*2);
M = sqrt(s)*v';
Y = unitVecs*M;

transEllipse = Y+ mu3;

subplot(3,1,ii)
scatter(samples3(:,1), samples3(:,2))
hold on
plot(transEllipse(:,1),transEllipse(:,2), 'LineWidth',2);
title(['Random white noise transformed by random mu and cov'])
xlabel('X')
ylabel('Y')
box off
set(gca, 'TickDir', 'out')

end 

% The ellipse captures the shape of the data and will appear smaller than
% the total parameter of the data cloud because the data cloud is a
% gaussian distribution and will continue away from the mean with lower 
% proability the further it gets. The circle itself represents the highest 
% density of the gaussian distribution. 

%% PART D

%%
%* How would you, mathematically, compute the direction (unit vector) that maximizes
% the variance of the marginal distribution? Compute this direction and verify that it is
% consistent with your plot.


% The principle component represents the greatest variance so a unit vector laying
% in that direction would capture the fullest amount of the variance in the
% sample. 
[u,s,v] = svd(CY);
unit1(1,:) = v(:,1);

figure
scatter(samples2(:,1), samples2(:,2),'b');
hold on 

%Plot vector for maximum variance 
xlim = get(gca,'XLim');
m = (unit1(1,2)/unit1(1,1));
y1 = m*xlim(1);
y2 = m*xlim(2);
L = line([xlim(1) xlim(2)],[y1 y2],'LineWidth',2);
set(L, 'Color','r')
hold off
axis square
title(['Vector Representing Greatest Variance'])
legend('Transformed Random Data', 'Line of Greatest Variance')
xlabel('X')
ylabel('Y')
box off
set(gca, 'TickDir', 'out')

% In the plot the vector direction is consistant with the axis of greatest
% variance. 