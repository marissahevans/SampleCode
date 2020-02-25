%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 2 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 1

% Date intialized: Sept. 27th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Trichromacy
% Load the file colMatch.mat in your MATLAB environment. 
% This file contains matrices and vectors related to the color matching 
% experiment presented in class. In particular, the variable P is an N × 3 
% matrix containing wavelength spectra for three 'primary' lights, that could 
% be used in a color-matching experiment. For these problems
% N = 31, corresponding to samples of the visible wavelength spectrum from 
% 400nm to 700nm in increments of 10nm. The function humanColorMatcher.p 
% simulates a normal human observer in a color matching experiment. The input 
% variable light should contain the wavelength spectrum of a test
% light (a 31-dimensional column vector). The variable primaries should contain 
% the wavelength spectra of a set of primary lights (typically, a 31 × 3 matrix, 
% as for matrix P described above). The function returns a 3-vector containing 
% the observer's 'knob settings' - the intensities of each of the primaries 
% that, when mixed together, appear identical to the test light. The function 
% can also be called with more than one test light (by passing a matrix
% whose columns contain 31-dimensional test lights), in which case it returns a matrix whose
% columns are the knob settings corresponding to each test light.

load('colMatch.mat');

%% PART A
% (a) Create a test light with an arbitrary wavelength spectrum, by generating 
% a random column vector with 31 positive components (use rand). Use humanColorMatcher to
% 'run an experiment', asking the 'human' to set the intensities of the three primaries in
% P to match the appearance of the test light. Compute the 31-dimensional wavelength
% spectrum of this combination of primaries, plot it together with the original 
% light spectrum, and explain why the two spectra are so different, even though they appear the
% same to the human.

lightRand = rand(31,1); % 31x1 column test vector
knobSetting = humanColorMatcher(lightRand, P); %3 value knob output

testWave = P*knobSetting; %31x1 light output based on primary knobs

figure
plot(lightRand, 'LineWidth', 2)
hold on 
plot(testWave, 'LineWidth',2)
title('Plot wave w/ Original Spectrum')
xlabel('Wavelength Spectrum')
ylabel('Output Level')
legend('Original Light Specturm', 'Spectrum of Primaries')
box off
set(gca, 'TickDir', 'out')

disp('The original spectrum has much more variablilty that is smoothed out')
disp('when computed only using the 3x knob settings. This shows that there is')
disp('much more variability in the color spectrum than is perceived by the human')
disp('eye using only three cones. Using these three light sensors humans can')
disp('repicalte their view of the world in accordance with how it is neurally processed')

%% PART B
% (b) Now characterize the human observer as a linear system that maps 31-dimensional
% lights to 3-dimensional knob settings. Specifially, run a set of experiments to estimate
% the contents of a 3 × 31 color-matching matrix M that can predict the human responses.
% Verify on a few random test lights that this matrix exactly predicts the responses of the
% function humanColorMatcher.

disp('We are using a loop to move though each of the 31 points on the spectrum')
disp('individually do calculate the exact "human" output for each point. Then')
disp('these 3x value outputs are concatonated into a matrix to use for analysis')
disp('of additional light sources.')

for ii = 1:31
testLight = zeros(31,1); % 31x1 column vector
testLight(ii) = 1;
humanOutput = humanColorMatcher(testLight, P); %3 value knob output

M(:,ii) = humanOutput; %makes a matrix M of the human predictors 
end

disp('Our matrix creates the same output as the human color matcher (confirmed')
disp('by the "verify" checks below) therefore our manipulation was a success.')

%Testing with random number light vector inputs
testLight2 = rand(31,1);
check2 = M*testLight2;
testOut2 = humanColorMatcher(testLight2, P);
verify2 = round(check2,4) == round(testOut2,4)

testLight3 = rand(31,1);
check3 = M*testLight3;
testOut3 = humanColorMatcher(testLight3, P);
verify3 = round(check3,4) == round(testOut3,4)

testLight4 = rand(31,1);
check4 = M*testLight4;
testOut4 = humanColorMatcher(testLight4, P);
verify4 = round(check4,4) == round(testOut4,4)

%% PART C
% (c) The variable Cones contains (in the rows) approximate spectral sensitivities of the
% three color photoreceptors (cones) in the human eye: Cones(1,:) is for the 
% L (longwavelength, or red) cones, Cones(2,:) the M (green) cones, and Cones(3,:) the
% S (blue) cones. Applying the matrix Cones to any light ~l yields a 3-vector containing
% the average number of photons absorbed by that cone (try plot(Cones') to visualize them!). 

figure
plot(Cones(1,:),'r','LineWidth',2)
hold on
plot(Cones(2,:),'g','LineWidth',2)
plot(Cones(3,:),'b','LineWidth', 2)
title('Cone Responses')
xlabel('Wavelength Spectrum')
ylabel('Output Level')
legend('Long Wave', 'Medium Wave', 'Short Wave')
box off
set(gca, 'TickDir', 'out')


%%
%* Verify that the cones provide a physiological explanation for the 
% matching experiment, in that the cone absorptions are equal for any pair of lights 
% that are perceptually matched. First, do this informally, by checking that randomly generated lights and
% their corresponding 3-primary matching lights produce equal cone absorptions.

lightRand = rand(31,1); %random light wavelength
colorPrime = P*M*lightRand; %light wavelength outputs based on primes

coneLightTest = Cones*lightRand; %cone absorbtion for random light
conePrimeTest = Cones*colorPrime;%cone absorbtion for primaries 

coneTest = round(coneLightTest,4) == round(conePrimeTest,4)

disp('The two cone responses from the random light and the percepturally matched')
disp('light are equal as shown by "coneTest" logical above. This shows that the')
disp('cones provide a psysiological explanation for the matching experiment')


%%
%* Then, provide a few lines of matlab code that provide a more mathematical demonstration,
% along with an extended comment explaining your reasoning using concepts from linear
% algebra. [Hints for two possible approaches: (1) write math/code that computes cone
% responses for any test light and then computes the weighted combination of primaries
% that would produce the same cone responses - show that this is numerically the same
% as the color-matching matrix; (2) convince yourself, and explain why, it is sufficient to
% show that M and Cones have the same nullspace. Then use SVD to demonstrate that
% this is true!]

%RESPONSE
disp('Because there are only 3x color receptiors in the eye any additional')
disp('values drop to the null space as the data is not sent further into the')
disp('brain. The three receptors cover the visible spectrum (as shown in the')
disp('plot above) so the weighted response of a combination between all 3x of')
disp('them will allow for a point to be reached anywhere on the visible')
disp('spectrum. Because we only recieve signals from 3x inputs we are able to')
disp('duplicate these signals using only 3x lights that match the colors the')
disp('Cones are designed to pick up. Since we are only using 3x knobs the other')
disp('information from the test light is dropped to the null space as well.') 

[UCone SCone VCone] = svd(Cones);
[UM SM VM] = svd(M);

disp('The SCone Matrix')
disp(SCone')
disp('The SM Matrix')
disp(SM')

disp('Both SCone and SM both are 3x31 diagonal matricies and they only have values in the')
disp('first three columns and zeros in the rest of the matrix represnting the')
disp('null space. Only the light signal from the top 3x rows of the V matrix')
disp('will make it further into the calculation of the light perception. The')
disp('values absorbed by each cone will be what appears in the range space, just')
disp('as the values that are represnted in the knobs will be in the range space')
disp('for the M matrix as well. ')


%% PART D
% (d) The function altHumanColorMatcher(light,primaries) simulates a color-deficient
% human observer in a standard color matching experiment. 

%%
%* (i) for a random test light,
% compare the knob settings for this observer with those of the normal human. Do this for
% several runs of altHumanColorMatcher(light,primaries). How do they differ?

lightRand = rand(31,1); %random light wavelength
altOutput = altHumanColorMatcher(lightRand, P) %dificent human response to light
regOutput = humanColorMatcher(lightRand, P) %regular human response to light

disp('The color deficient human is unable to see one of either medium or short')
disp('wavelengths because the the long wavelengths are fairly consistent but')
disp('the medium and short are over manipulated showing they are using one to')
disp('compensate for the other ')


%%
%* (ii) Compute cone absorptions for the test light, and for the mixture of three matching
% primaries (by applying the Cones matrix). Do this for both the normal human observer,
% and for multiple runs of the abnormal observer. Try this for several different test lights.
% How do the cone responses of the normal and abnormal observers differ? Can you offer
% a diagnosis of the underlying cause of color deficiency in the abnormal observer?

figure
coneOpt = [1;2;3;];

for jj = 1:6
lightRand = rand(31,1);
regOutput = humanColorMatcher(lightRand, P);
regConeAbs = Cones * P * regOutput;
%regConeAbs = regConeAbs*P;

for ii = 1:10
    altOutput = altHumanColorMatcher(lightRand, P);
    altConeAbs(ii,:) = Cones * P * altOutput;
    %altConesAbs(ii,:) = altConeAbs(ii,:)*P;
end
subplot(2,3,jj)
plot(coneOpt,altConeAbs', 'g')
hold on
plot(coneOpt,regConeAbs')
axis([.5 3.5 0 12])
xlabel('Cone L/M/S')
ylabel('Cone Response')
title(['Test Light ' num2str(jj)])
box off
set(gca, 'TickDir', 'out')
end

disp('By looking at the plots we can tell that the participant is color')
disp('deficient in the medium wavelength spectrum. We can see this because they over')
disp('compinsate with short and long to account for the detected luminance change of the')
disp('stimulus but are regularly under representing the green knob on the color')
disp('spectrum.')

disp('GREEN is the deficient human and PURPLE is the normal human (putting a legend on each subplot made it hard to read on my laptop')