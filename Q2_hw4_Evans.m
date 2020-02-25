%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 4 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 2

% Date intialized: Oct. 31th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Poisson neurons. 
%The Poisson distribution is commonly used to model neural spike counts:

% p(k) = u^k *exp(1)^?u/factorial(k)

%where k is the spike count (over some specified time interval), and � is the expected number
% of spikes over that interval.

%% PART A

%%
%*  We would like to know what the Poisson distribution looks like. Set the 
% expected number of spikes to � = 6 spikes/interval then create a vector p 
% of length 21, whose elements contain the probabilities of Poisson spike 
% counts for k = [0...20]. Since we?re clipping the range at a maximum value 
% of 20, you?ll need to normalize the vector so itsums to one (the distribution 
% given above is normalized over the range from 0 to infinity) to make the 
% vector p represent a valid probability distribution. Plot p in a bar plot
% and mark the mean firing rate. Is it equal to �?

p = nan*ones(21,1);
u = 6;
k = 0:20;
for kk = 0:20
    p(kk+1,1) = u^kk*exp(1)^-u/factorial(kk);
end 

p = p/sum(p); %normalize 

spikeMean = sum(p.*k');

figure
bar(k,p)
title('PDF of Poisson Distribution, � = 6')
xlabel('spikes/interval')
ylabel('probability')
xline(spikeMean,'--r', 'LineWidth',2);
legend('Data', 'Mean')
box off
set(gca, 'TickDir', 'out')

% The caluclated mean firing rate is equal to � (6 spikes/interval)

%% PART B

%% 
%* Generate samples from the Poisson distribution where each sample 
% represents the number of spike count ranging from 0 to 20. To simplify 
% the problem, use a clipped Poisson vector p to write a function 
% samples = randp(p, num) that generates num samples from the probability 
% distribution function (PDF) specified by p. [Hint: use the rand 
% function, which generates real values over the interval [0...1], and 
% partition this interval into portions proportional in size to the probabilities in p]. 
% Test your function by drawing 1,000 samples from the Poisson distribution in (a), plotting a 
% histogram of how many times each value is sampled, and comparing this to the frequencies 
% predicted by p. 

numSamples = 1000

samples = randp(p,numSamples,1); %Testing function with a sample of 1000

for ii = 0:20
    sampleProb(ii+1,1) = (sum(samples == ii))/numSamples;
end

sampleProb(:,2) = p;
disp(sampleProb)

%Looking at the probabilities of p (column 2) compared to the proportion of
%spikes which occured in each interval (column 1) they are very close, (confirmed
%visually by the plot as well). As the total spike count increases the frequency values
%converge toward an equality. 

%%
%* Verify qualitatively that the answer gets closer (converges) as you 
% increase the number of samples (try 10 raised to powers [2, 3, 4, 5]).

numSamp = [10^2, 10^3, 10^4, 10^5];
testPVec(:,1) = p;

for jj = 1:length(numSamp)
    sampleNum = numSamp(jj);
    samples = randp(p,sampleNum,1);
    for ii = 0:length(p)-1
        testProb = length(find(samples == ii));
        testPVec(ii+1,jj+1) = testProb/sampleNum;
    end
end

disp(testPVec)

testSimilarity = round(p,2) == round(testPVec(:,5),2)
%As the number of samples increases to 10^5, the two probabilities are
%almost all equal until the tenths place, and will continue to converge toward
%each ohter as the spike count increases. 

%% PART C

%% 
%*  Imagine you?re recording with an electrode from two neurons simultaneously, whose
% spikes have very similar waveforms (and thus can?t be distinguished by the spike sorting software). 
% Create a probability vector, q, for the second neuron, assuming a mean
% rate of 4 spikes/interval. 

q = nan*ones(21,1);
u = 4;
k = 0:20;
for kk = 0:20
    q(kk+1,1) = u^kk*exp(1)^-u/factorial(kk);
end 
q = q/sum(q); %normalize 

%%
%* What is the PDF of the observed spike counts, which will be
% the sum of spike counts from the two neurons derived from p and q? [Hint: the output
% vector should have length m + n ? 1 when m and n are the lengths of the two input
% PDFs. This is because the maximum spike count will be bigger than the maximum of
% each respective individual neuron.]

% Verify your answer by comparing it to the histogram of 1,000 samples generated by
% summing two calls to randp.


PDFcombined = conv(q,p);

k = 0:length(p)+length(q)-1;
numSamples = 1000;

sample1 = randp(p, numSamples, 0);
sample2 = randp(q, numSamples, 0);

combined = sample1 + sample2;


figure
yyaxis left
bar(PDFcombined)
title(['PDF of observed spike counts and Combined Poisson mu = 4+6 Distribution'])
xlabel('Spikes/Interval')
ylabel('Probability')
box off
set(gca, 'TickDir', 'out')

yyaxis right
histogram(combined,k)
xlabel('Spikes/Interval')
ylabel('# of Samples')
box off
set(gca, 'TickDir', 'out')
legend('PDF', 'Poisson Dist')


% The PDF of the distribution and the combined sample spike outputs have
% relatively the same size, however the PDF is more standardized given we
% do not have an infinite number of samples. This is especially noticible
% in the tails as the probability of having sample data present there is
% much lower

%% PART D

%%
%*  Now imagine you are recording from a neuron with mean rate 10 spikes/interval (the
% sum of the rates from the neurons above). 

n = nan*ones(21,1);
u = 10;
k = 0:20;
for kk = 0:20
    n(kk+1,1) = u^kk*exp(1)^-u/factorial(kk);
end 
n = n/sum(n); %normalize 

%%
%* Plot the distribution of spike counts for this
% neuron, in comparison with the distribution of the sum of the previous two neurons.


sample4 = randp(n,numSamples,0);

figure
subplot(1,2,1)
histogram(combined)
title(['Poisson Distribution with mean 4 + 6'])
xlabel('Spikes/Interval')
ylabel('# of Samples')
box off
set(gca, 'TickDir', 'out')
subplot(1,2,2)
histogram(sample4)
title(['Poisson Distribution with mean 10'])
xlabel('Spikes/Interval')
ylabel('# of Samples')
box off
set(gca, 'TickDir', 'out')


%%
%* Based on the results of these two experiments, if we record a new spike train, can you
% tell whether the spikes you have recorded came from one or two neurons just by looking
% at their distribution of spike counts? Comment about the reason why based on the
% intuition behind Poisson distribution

% As u increases in the Poisson distribution the shape becomes wider and
% flatter. The two distributions look very similar (never exactly because
% we're using random samples) because when mean=10 the density is distributed 
% over more values leading to a wider curve than if the mean=2 for example. 
% When summing the two means (4 + 6 in this case) taking Poissons with 
% PDF's of x1, x2,... xn and means of u1, u2,...un, then x1 + x2 +...+ Xn 
% is a Poisson distribution too, with mean parameter u1 + u2 +...+ un.

% The primary difference in our graphs as we are using a truncated poisson
% distribution is that the dispersion for u=10 will never exceed a value of
% 20, while the distributions of 4+6 are able to extend past 20 if both
% have values present in those areas. 
