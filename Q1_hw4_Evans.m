%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 4 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 1

% Date intialized: Oct. 31th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Bayes? rule and eye color. 

% A male and female chimpanzee have blue and brown eyes, respectively. The 
% brown-eyed allele can be denoted as a capital B, whereas the blue-eyed allele
% can be represented as a lowercase b. Assume a simple genetic model in which the gene for
% brown eyes is always dominant (so that the trait of blue eyes can only arise from 
% two blueeyed genes, but the trait of brown eyes can arise from two brown-eyed genes, or one of each).
% You can also assume: i) the probability of the mother being BB is 50% and the probability of
% her being Bb is 50%; and ii) the a priori probability that each of the four gene configurations
% is equally probable. For each question, provide the math, and explain your reasoning.

%% PART A

%%
%* Suppose you observe that they have a single child with brown eyes. What
% is the probability that the female chimp has a blue-eyed gene?

% male = bb
% female = BB or Bb

% possible combinations: BBxbb = Bb, Bb, Bb, Bb & Bbxbb = Bb, Bb, bb, bb

% There are 8x total combinations, 6x lead to brown eyes and 2x lead to
% blue. 

% Since only 2x of those 6x combinations resulted from the mother having a
% blue gene the probability is 1/3 

% P(B) = event of child having brown eyes
    % 50% of the time the outcome will ALWAYS be brown 
    % 50% of the time the outcome will be brown 50% of the time
% P(A) = event of mother having blue gene
% P(B|A) = event of brown eyes occuring if blue gene present

% P(B|A)*P(A) / P(B)

probBrown = .5 * .5 + .5 * 1; % P(B) = event of child having brown eyes
probBlueGene = .5; % P(A) = event of mother having blue gene
eventBrownIfBlue = .5; % P(B|A) = event of brown eyes occuring if blue gene present

probBlueGene = eventBrownIfBlue*probBlueGene/probBrown;

disp(probBlueGene)

% P(A|B) = .5^2 / 3/4


%% PART B

%%
%* Suppose you observe that they have a second child with brown eyes. Now what is the
% probability?

% P(C) = event of 2x children having brown eyes
    % for each child 50% of the time the outcome will ALWAYS be brown 
    % for each child 50% of the time the outcome will be brown 50% of the time
% P(A) = event of mother having blue gene
% P(C|A) = event of brown eyes occuring in both children if blue gene present

% P(C|A)*P(A) / P(C)

probBrown2 = (.5 * .5 + .5 * 1)^2; % P(C) event of 2 children having brown eyes
probBlueGene = .5; % P(A) event of mother having blue gene
eventBrownIfBlue2 = .5^2; % P(C|A) event of brown eye occuring twice if blue gene present

probBlueGene = eventBrownIfBlue2*probBlueGene/probBrown2;

disp(probBlueGene)

% P(A|C) = .5^3 / 3/4^2


%% PART C

%%
%*  Generalizing, suppose they have N children with brown eyes... express the probability,
% as a function of N.

% P(N) = event of Nx children having brown eyes
    % for each child 50% of the time the outcome will ALWAYS be brown 
    % for each child 50% of the time the outcome will be brown 50% of the time
% P(A) = event of mother having blue gene
% P(N|A) = event of brown eyes occuring in Nx children if blue gene present

% P(N|A)*P(A) / P(C)

N = 10; % N = number of children

probBrown3 = (.5 * .5 + .5 * 1)^N; % P(N) event of N children having brown eyes
probBlueGene = .5; % P(A) event of mother having blue gene
eventBrownIfBlue3 = .5^N; % P(N|A) event of brown eye occuring N times if blue gene present

probBlueGene = eventBrownIfBlue3*probBlueGene/probBrown3;

disp(probBlueGene)

% P(A|N) = .5^N+1 / 3/4^N    Expression in terms of 'N' 
