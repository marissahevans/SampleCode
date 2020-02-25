%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 5 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 1

% Date intialized: Nov. 20th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Dueling estimators. 
% In this problem, we use simulation to compare three estimators of the
% mean of a Normal (Gaussian) distribution.

%% PART A

%%
%* First consider the average, which minimizes the sum of squared deviations, and is also
% the Maximum Likelihood estimator. Generate 10,000 samples, each of size 10, from the
% Normal(0,1) distribution (a 10x10000 matrix). Compute the average of each of the 10,000
% samples. Plot a histogram of the resulting estimates (use 50 bins, and set the plot range to
% [-2.3,2.3]). What shape should the histogram have (explain why)?

samples = randn(10000,10);
meanSamples = nan(10000,1);

for ii = 1: length(samples)
    meanSamples(ii) = mean(samples(ii,:));
end

figure
histogram(meanSamples,50)
xlim([-2.3,2.3])
title('Mean of Samples')
xlabel('Mean')
ylabel('# Samples')
box off
set(gca, 'TickDir', 'out')

% The histogram should be bell shaped like a normal distribution because we
% are sampling from a normal distribution so the mean should follow a
% similar structure. Since the distribution we're sampling from has a mean
% of 0 and an SD of 1 our histogram should follow the same metric. 

%%
%*  What is the (theoretical) variance of the average of 10 values drawn 
% from a univariate Gaussian (derive this)? Is the empirical variance of 
% your 10,000 estimates close to this?

% Var(A1 + A2 ... An) = Var(A1) + Var(A2)... + Var(An)
% Var(alpha A) = alpha^2 * Var(A)
% (1/N) * Var(A) = Var(A)/N ----> theoreticalVar(samples) = var(samples)/10
% thorVar = 1/Nsamples

for ii = 1:10000
   sampleVar(ii) = var(samples(ii,:))/10;
end

empiricalVar = mean(sampleVar)

% As shown in the equation derived above the theoretical variance should be
% .1 as there are 10 values. Our empirical variance is very close to this
% value as it is the average of multiple sets of ten values drawn randomly 
% from a univariant guassian. 

%% PART B

%%
%* Now consider the median, which minimizes the sum of absolute deviations. Compute
% the median of each of the 10,000 samples, and again plot a histogram. What shape does this
% one have? 

for ii = 1:10000
medianSamp(ii) = median(samples(ii,:));
end

figure
histogram(medianSamp,50)
xlim([-2.3,2.3])
title('Median of Samples')
xlabel('Median')
ylabel('# Samples')
box off
set(gca, 'TickDir', 'out')

%The shape of this distribution is also bell shaped and reflects a normal
%distribution. It is very close in shape to the distribution of the means,
%which makes sense because since the data was drawn from a normal
%distribution with mean 0 and SD 1 the average and center of the data
%would likely be very close together on most trials. 

%%
%* Compare it to a normal distribution using the function normplot, which plots the
% quantiles of a sample of data versus the normal quantiles (known as a Q-Q plot: if data are
% normally distributed, the points shuld fall nearly on a straight line.) Does the distribution of
% estimated values deviate significantly from a Normal distribution? Specifically, compare the
% Q-Q plot for the median estimator to that for the mean from part (a).

figure
subplot(1,2,1)
normplot(medianSamp)
xlim([-2.3,2.3])
title('Normplot Median')
xlabel('Median')
ylabel('# Samples')
box off
set(gca, 'TickDir', 'out')

subplot(1,2,2)
normplot(meanSamples)
xlim([-2.3,2.3])
title('Normplot Mean')
xlabel('Mean')
ylabel('# Samples')
box off
set(gca, 'TickDir', 'out')

% The plot does not deviate from the normal distribution. Both the mean and the
% median are good estimators for the data as shown by the QQ plot's red
% line. 

%% PART C

%%
%* Finally, consider an estimator that computes the average of the minimum and maximum
% over the sample (as shown in class, this one minimizes the L-norm). Again, compute this
% estimate for each of your 10,000 samples, plot the histogram, and examine and comment on
% the Q-Q plot, just as in part (b).

for ii = 1:10000
    
   maxSample = max(samples(ii,:));
   minSample = min(samples(ii,:));
   
   extremaSamples(ii) = (maxSample + minSample)/2;
end

figure
subplot(1,2,1)
histogram(extremaSamples,50)
xlim([-2.3,2.3])
title('Extrema Norm')
xlabel('Average')
ylabel('# Samples')
box off
set(gca, 'TickDir', 'out')

subplot(1,2,2)
normplot(extremaSamples)
xlim([-2.3,2.3])
title('Normplot Extrema')
xlabel('Average')
ylabel('# Samples')
box off
set(gca, 'TickDir', 'out')

% The extrema, or center of data based on the furthest outlying values, is
% also represented as a normal distribution in this case. Because the data
% is drawn from a normal distribution (m 0 sd 1) the edges of the data
% should be evenly dispursed across trials, but centering at around the
% same point as the mean and the median. Looking at the QQ plot the
% majority of the distribution is normal however at the edges it is skewed
% a little off normal (just barely). 

%% PART D

%%
%* All three of these estimators are unbiased (because of the symmetry of the distribution),
% so we can use variance as the sole criterion for quality. Generate a new set of 10,000 samples, 
% this time of dimension 256. Apply each estimator to sub-matrices of samples of size
% {8, 16, 32, 64, 128, 256}, and compute the variance of each estimator for each. 

samples2 = randn(256, 10000);

subsamples = [8 16 32 64 128 256];

for ii = 1: length(subsamples)
    jj = subsamples(ii);
    mean2(ii) = var(mean(samples2(1:jj,:)));
    median2(ii) = var(median(samples2(1:jj,:)));
    extrema2(ii) = var((max(samples2(1:jj,:)) + min(samples2(1:jj,:)))/2);
    theorVar(ii) = 1/(jj);
end


%%
%* Plot these (on a single log-log plot) along with a line showing the 
% theoretically-computed variance of the average estimator. 

figure
loglog(subsamples, mean2, 'LineWidth',2)
hold on
loglog(subsamples,median2, 'LineWidth',2)
loglog(subsamples,extrema2, 'LineWidth',2)
loglog(subsamples,theorVar,'r--','LineWidth',2)
yline(extrema2(6),'b--');
xline(8,'--'); xline(16,'--'); xline(32,'--');
xline(64,'--'); xline(128,'--'); xline(256,'--');
xticks([8 16 32 64 128 256])
set(gca,'xscale','log')
set(gca,'yscale','log')
title('Plot of all Estimators')
xlabel('Sample Size (log)')
ylabel('Variance (log)')
legend('Variance of the Means', 'Variance of the Medians', 'Variance of the Extrema','Theoretical Var','Location','bestoutside')
box off
set(gca, 'TickDir', 'out')

%%
%* % Does the variance of all three estimators converge at the same rate (1/N)?
% How much larger is the variance of the median estimator than the average estimator? How
% large a sample would you need for the average and median estimators to achieve the same
% variance as the average-extrema estimator (from part (c)) on samples of size 256?

% The estimators do not converge at the same rate, the extrema is very
% different and the median has a very slight almost negligable difference
% from the average. 

mDif = median2 - mean2

mRatio = median2./mean2

% The median estimator on average starts at around .04 larger than the
% average estimator and they begin to converge with a final difference of
% .002, however when plotted on a log scale this convergence is less
% noticible. 

% From looking at the plot the mean would need a sample size of around 12
% to have the same variance as the extrema at a sample size of 256

% The ratio between the median and the mean is usually about 1.5 meaning
% the variance in the median is usually 50% more than in the mean. 
