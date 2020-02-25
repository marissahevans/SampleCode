%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 5 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 4

% Date intialized: Nov. 20th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Simulating a 2AFC experiment. 
% Consider a two-alternative forced choice (2AFC) psychophysical experiment in
% which a subject sees two stimulus arrays of some intensity on a
% trial and must say which one contains the target. (One and only one contains the target.)
% Her probability of being correct on a trial is:

% pc(I) = 1/2 + 1/2?(I; µ, ?)

% where ?(I; µ, ?) is the cumulative distribution function of the Gaussian (normcdf in matlab)
% with mean µ and standard deviation ? evaluated at I. The function pc(I) is known as the
% psychometric function. (Minor note, somewhat subtle: This setup only makes sense if I is on
% a logarithmic scale, e.g., I = k log C, where C is stimulus contrast.)

%% PART A

%%
%* (a) Plot two psychometric functions, for {µ, ?} equal to {5, 2} and {4, 3}. (Use I = [1 : 10]).
% Describe the difference between these. 

mu1 = 5;
mu2 = 4;
sd1 = 2;
sd2 = 3;
I = [1:10];

psychFun1 = .5 + .5*(normcdf(I,mu1,sd1));
psychFun2 = .5 + .5*(normcdf(I,mu2,sd2));

figure
plot(I, psychFun1, 'LineWidth',2)
hold on
plot(I, psychFun2,'LineWidth',2)
title('Plot of 2 psychometric functions')
xlabel('Stimulus Intensity (log)')
ylabel('Probability Correct Resp')
set(gca,'xscale','log')
legend('Psychometric Function 1','Psychometric Function 2','Location','bestoutside')
box off
set(gca, 'TickDir', 'out')

%Function 1 has a more 's' like curve and has the lower SD but it is shifted
%to the right more than funtion 2 because it has a higher mean, which also 
% give it a steeper slope. Function 2 is straighter but has a more gentle
% slope with a higher y intercept. 

%%
%* If you increase µ, how does the curve change?
% If you increase sigma, how does the curve change? (If you are not sure, make more plots
% with different parameter values.) 

figure
for ii = 1:4

cdfGaus1 = normcdf(I,mu1+ii,sd1);
cdfGaus2 = normcdf(I,mu2+ii,sd2);

psychFun1 = .5 + .5*cdfGaus1;
psychFun2 = .5 + .5*cdfGaus2;

subplot(2,4,ii)
plot(I, psychFun1, 'LineWidth',2)
hold on
plot(I, psychFun2,'LineWidth',2)
title(['mean +' num2str(ii)])
xlabel('Stimulus Intensity')
ylabel('Probability Correct Resp')
ylim([.5 1])
set(gca,'xscale','log')
box off
set(gca, 'TickDir', 'out')

end

for ii = 1:4

cdfGaus1 = normcdf(I,mu1,sd1+ii);
cdfGaus2 = normcdf(I,mu2,sd2+ii);

psychFun1 = .5 + .5*cdfGaus1;
psychFun2 = .5 + .5*cdfGaus2;

subplot(2,4,ii+4)
plot(I, psychFun1, 'LineWidth',2)
hold on
plot(I, psychFun2,'LineWidth',2)
title(['SD + ' num2str(ii)])
xlabel('Stimulus Intensity')
ylabel('Probability Correct Resp')
ylim([.5 1])
set(gca,'xscale','log')
box off
set(gca, 'TickDir', 'out')
end

% As the mean increases the curve becomes less of an 's' curve and more of
% an exponential increase curve. The point where the slope changes from a
% gradual increase to a steep increase is moved further in the direction of
% th new mean. When the mean is increased only very intense stimuli become
% easy to reliably dectect and the performance plateau at the top of the
% curve is illiminated. 

% As the standard deviation increases the curves flatten out to become
% stright lines removing any pleateaus in the response accuracy rate. The y
% intercept of the curve is increased as well showing that there are fewer
% very poor responses as well as fewer very good responses. Most of the
% performance is centered in the middle of the accuracy range which would
% suggest the participant is not effected much by changes in the stimuli. 
%%
%* What is the range of pc(I)? Explain why this range is
% appropriate.

%The range of pc(I) is .5-1, which makes sense because it is a measure of
%the probability of the participant correctly identifying the stimulus at a
%range of intensity levels above chance. Because this is a forced choice
%experiment the minimum performance level should be at 50% because if a
%participant was reliably performing at a lower accuracy rate than that it
%would appear that they are sensitive to the simulus but did not understand
%the experiment instructions. 

%% PART B

%%
%* Write a function C=simpsych(mu,sigma,I,T) which takes two vectors (I,T) of the
% same length, containing a list of intensities and the number of trials for each intensity,
% respectively, simulates draws from pc(I), and returns a vector, C, of the same length as
% I and T, which contains the number of trials correct out of T, at each intensity I.

%% PART C

%% 
%* Illustrate the use of simpsych with T=ones(1,7)*100 and I=1:7 for µ = 4 and ? = 1.
% Plot C ./ T vs I (as points) and plot the psychometric function pc(I) (as a curve) on
% the same graph.

T = ones(1,7)*100;
I = 1:7;
mu = 4;
sigma = 1;

C = simpsych(mu,sigma,I,T);

pcI = .5 + .5*(normcdf(I,mu,sigma));

figure
plot(pcI, 'LineWidth',2)
hold on
plot(I,C./T,'ro','LineWidth',2)
xlabel('Stim Intensity')
ylabel('Prob Correct Resp')
title('Plot of 100 trials per intensity')
legend('pc(I)','simpsych points','Location','bestoutside')
set(gca,'xscale','log')
box off
set(gca, 'TickDir', 'out')


%% PART D

%%
%* Do the same with T=ones(1,7)*10 and plot the results (including the psychometric
% function). What is the difference between this and the plot of the previous question?

T = ones(1,7)*10;
I = 1:7;
mu = 4;
sigma = 1;

C = simpsych(mu,sigma,I,T);

pcI = .5 + .5*(normcdf(I,mu,sigma));

figure
plot(pcI, 'LineWidth',2)
hold on
plot(I,C./T,'ro','LineWidth',2)
xlabel('Stim Intensity')
ylabel('Prob Correct Resp')
title('Plot of 10 trials per intensity')
legend('pc(I)','simpsych points','Location','bestoutside')
set(gca,'xscale','log')
box off
set(gca, 'TickDir', 'out')

% When using fewer trial values (10 in this case) it's impossible to get
% the test points exactly on the parameter line especially in the center of
% the data set as there is much more room for variablility. When we use 100
% trials per intensity level our simulated participant scores much closer
% to the expected value of the psychometric curve on average. When the
% stimulus in intense the values will be much more relable because there is
% inherently less variability in the possible responses. 