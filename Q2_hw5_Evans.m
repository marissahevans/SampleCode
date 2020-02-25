%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 5 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 2

% Date intialized: Nov. 20th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Bayesian inference of binomial proportions. 
% Poldrack (2006) published an influential attack on the practice of 'reverse inference' 
% in fMRI studies, i.e. inferring that a cognitive process was engaged on the basis of 
% activation in some area. For instance, if Broca's area was found to be activated 
% using standard fMRI statistical-contrast techniques, researchers might
% infer that the subjects were using language. In a search of the literature, Poldrack found that
% Broca's area was reported activated in 103 out of 869 fMRI contrasts involving engagement
% of language, but this area was also active in 199 out of 2353 contrasts not involving language.

% P(Lang|Broca) = P(Broca|Lang) * P(Lang) / P(Broca)

% P(Broca) = P(Broca|Lang) * P(Lang) + P(Broca|NoLang) * P(NoLang)

%% PART A

%%
%* Assume that the conditional probability of activation given language, as well as that of
% activation given no language, each follow a Bernoulli distribution (i.e., like coin-flipping),
% with parameters xl and xnl. Compute the likelihoods of these parameters, given Poldrack's
% observed frequencies of activation. Compute these functions at the values x=[0:.001:1] and
% plot them as a bar chart.

x1 = [0:.001:1];
x2 = [0:.001:1];
k1 = 103;
n1 = 869;
k2 = 199;
n2 = 2353;

likeLang = nchoosek(n1,k1)*(x1.^k1.*(1-x1).^(n1-k1));
likeNoLang = nchoosek(n2,k2)*(x2.^k2.*(1-x2).^(n2-k2));

figure
bar(x1,likeLang)
hold on
bar(x2,likeNoLang)
title('Likelihood of Broca Activation given Language or not')
xlabel('X')
ylabel('Likelihood')
set(gca, 'TickDir', 'out')
legend('Likelihood Language','Likeihood No Language')
box off

%% PART B

%%
%* Find the value of x that maximizes each discretized likelihood function.
% Compare these to the exact maximum likelihood estimates given by the 
% formula for the ML estimator of a Bernoulli probability.

[~, I1] = max(likeLang)
[~, I2] = max(likeNoLang)

x1Ind = x1(I1)
x2Ind = x2(I2)

ML1 = 103/869
ML2 = 199/2353

%The values are very similar with the ML estimator being slightly smaller
%on both accounts, but the closest available value of the parameter. 

%% PART C

%%
%* Using the likelihood functions computed for discrete x, compute and plot the discrete
% posterior distributions P(x | data) and the associated cumulative distributions P(X <= x |
% data) for both processes. For this, assume a uniform prior P(x) ~ 1 and note that it will
% be necessary to compute (rather than ignore) the normalizing constant for Bayes' rule. 
% Use the cumulative distributions to compute (discrete approximations to) 
% upper and lower 95% confidence bounds on each proportion. 

N1 = sum(likeLang);
N2 = sum(likeNoLang);

posterior1 = likeLang.* 1 / N1;
posterior2 = likeNoLang.* 1 /N2;

%Cumulative Dist. 
for ii = 1:length(posterior1)
    cdf1(ii) = sum(posterior1(1:ii));
    cdf2(ii) = sum(posterior2(1:ii));
end

figure
subplot(2,1,1)
bar(x1,posterior1)
hold on
bar(x2,posterior2)
title('PDF Dist')
xlabel('X')
ylabel('Posterior')
legend('Language','No Language','Location','bestoutside')
box off
set(gca, 'TickDir', 'out')

subplot(2,1,2)
plot(x1,cdf1, 'LineWidth',2)
hold on
plot(x2,cdf2, 'LineWidth',2)
yline(.975, '--','LineWidth',2);
yline(.025, '--','LineWidth',2);
title('CDF Dist')
xlabel('X')
ylabel('Cumulative Posterior')
legend('Language','No Language','Location','bestoutside')
box off
set(gca, 'TickDir', 'out')

% For the No Language conditon the lower & upper 95% CI are: .074 - .096
% For the Language condition the lower & upper 95% CI are: .099 - .14

% These values are estimated based on the CI lines on the CDF plot. 

%% PART D

%%
%* Are these frequencies different from one another? Consider the joint posterior distribution
% over xl and xnl, the Bernoulli probability parameters for the language and non-language
% contrasts. Given that these two frequencies are independent, the (discrete) joint distribution
% is given by the outer product of the two marginals. Plot it (with imagesc).


jointDist = posterior1' * posterior2;
figure
imagesc(x1, x2, jointDist)
colorbar
title('Heatmap of Joint Distribution')
box off
set(gca, 'TickDir', 'out')
xlabel('X')
ylabel('Y')


%Compute (by summing the appropriate entries in the joint distribution) 
%the posterior probabilities that xl > xnl and, conversely, that xl <= xnl.

testMat = zeros(size(jointDist));
for ii = 1:length(testMat)
    testMat(1:ii,ii) = 1;
end 

langLessOrEqual = sum(jointDist.*testMat,'all')
langGreater = 1 - langLessOrEqual

%These frequencies are different when looking at the output of joint dist
% > and <= than xnl. 

%% PART E

%%
%* Is this difference sufficient to support reverse inference? Compute the probability P(language |
% activation). This is the probability that observing activation in Broca's area implies engagement 
% of language processes. To do this use the estimates from part (b) as the relevant conditional 
% probabilities, and assuming the prior that a contrast engages language,
% P(language) = 0.5. 

norm = (k1 + k2) / (n1 + n2) 

posterior1a = (x1Ind* .5) /norm %given Broca's activation language was used
posterior2b = (x2Ind* .5) /norm %given Broca's activation NO language was used


% Poldrack's critique said that we cannot simply conclude that activation
% in a given area indicates that a cognitive process was engaged without computing the posterior
% probability. Is this critique correct? To answer this, compare the Bayes factor (probability
% of language vs. not language) after taking Poldrack's data into account, compared to before
% having done so. 

bayesFactor = posterior1a/posterior2b

% The difference is not significnat enough to support reverse inference. 
% Given the low base factor it is unlikely that language is important for
% Broca's activation. A Bayes' factor under 3.2 is not considered
% substantaial and our current value is only 1.4. There is not substantial
% enough evidence in this case to conclude that Broca's area is specific to
% language activation. 
