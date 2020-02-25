%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 3 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 1

% Date intialized: Oct. 15th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% LSI system characterization. 

% You are trying to experimentally characterize three auditory
% neurons, in terms of their responses to sounds. For purposes of this problem, the responses
% of these neurons are embodied in compiled matlab functions unknownSystemX.p, (with
% X=1, 2, 3), each of which takes an input column vector of length N = 64 whose elements
% represent sound pressure over time. The response of each is a column vector (of the same
% length) representing the mean spike count over time. Your task is to examine them to see if
% they: i) behave like they are linear; and/or ii) shift-invariant; with/without iii) circular (i.e.
% periodic) boundary handling . For each neuron,

%% PART A

%%
%* "Kick the tires" by measuring the response to an impulse in the first position of an
% input vector. 

testSound = zeros(64,1);
testSound(1) = 1;

sysTest(:,1) = unknownSystem1(testSound);
sysTest(:,2) = unknownSystem2(testSound);
sysTest(:,3) = unknownSystem3(testSound);

disp(sysTest)


%%
%* Check that this impulse response is shift-invariant by comparing to the
% response to an impulse at positions n = 2, 4, 8. 

set = [2 4 8];

for ii = 1:3
    sampNum = set(ii);
    testSound = zeros(64,1);
    testSound(sampNum) = 1;
    
    sysTest1(:,ii) = unknownSystem1(testSound);
    sysTest2(:,ii) = unknownSystem2(testSound);
    sysTest3(:,ii) = unknownSystem3(testSound);
    
end 

%The impulse response for system 1 and 2 are shift invariant because the output is always the
%same just shifted by a few spots relative to the position of n. The
%response for system 3 however increases with n so the power is different
%when shifted. Confirmed below:

verifyLin(:,1) = abs(sum(sysTest1(:,1)) - sum(sysTest1(:,3))) < 1e-8; %shift invariant
verifyLin(:,2) = abs(sum(sysTest2(:,1)) - sum(sysTest2(:,3))) < 1e-8; %shift invariant
verifyLin(:,3) = abs(sum(sysTest3(:,1)) - sum(sysTest3(:,3))) < 1e-8; %NOT shift invariant

disp('System 1 : System 2 : System 3')
disp('System 3 is not shift invariant')
disp(verifyLin)

%%
%* Use different n to determine how the system handles inputs near the boundary. 

impulse = zeros(64,1);

impulse(1) =1;
output1 = unknownSystem1(impulse);
output2 = unknownSystem2(impulse);
output3 = unknownSystem3(impulse);

impulse63 = circshift(impulse, 63);
output63_1 = unknownSystem1(impulse63);
output63_2 = unknownSystem2(impulse63);
output63_3 = unknownSystem3(impulse63);

verify(:,1) = abs(output63_1 - circshift(output1,63)) < 1e-8; %perodic boundary
verify(:,2) = abs(output63_2 - circshift(output2,63)) < 1e-8; %perodic boundary
verify(:,3) = abs(output63_3 - circshift(output3,63)) < 1e-8; %NOT a perodic boundary

disp('System 1 : System 2 : System 3')
disp('System 3 does not exhibit circular boundary handling')
disp(verify)

% Systems 1 and 2 use periodic boundary handling because inputs near the boundary (such as 63) 
% roll over to the top (or rather drop off the bottom and are replaced with
% new values at the top). System 3 does not have perodic boundary handling 
% because it lacks shift invariance. This is verified by comparing the shifted output
% against the orginal output, shifted. 

%%
%* Also check that the response to a sum of any two of these impulses is equal 
% to the sum of their individual responses. Be sure to describe your findings.

[StestSys1, StestSys2, StestSys3] = linCheckSum(20) 

%The function above checks the linearity property of additivity of 3x unknown 
%systems which are 64 rows in length. Input is how many random points of the 
%system you'd like to compare (10 points leads to 9 comparisons). 

%For systems 2 and 3 when randomly selecting any two response points, the 
%sum of the response to those two points run through the system separatly 
%is the same as the sum of both those responses run though the system 
%simultaniously, this is because the resting point of these neurons is 0 so
%when summed the total remains the same because you're just adding zero. In
%system 1 the responses do not sum to the same value because the resting
%value of this system is 1 meaning there is still a signal present. 

%% PART B

% If the previous tests succeeded, examine the response of the system to sinusoids with
% frequencies {2pi/N, 4pi/N, 8pi/N, 16pi/N}, and random phases, and check whether the
% outputs are sinusoids of the same frequency (i.e., verify that the output vector lies 
% completely in the subspace containing all the sinusoids of that frequency). [note: Make the
% input stimuli positive, by adding one to each sinusoid, and the responses should then
% be positive (mean spike counts)].

freqSet = [1 2 4 8];

for ii = 1:4
    
    k = freqSet(ii); %frequency
    n = 64; %total number of points
    x = 0:n-1; %length of signal domain
    A = 1; %amplitude
    phi = rand; %phase
    signalTest = A*cos(2*pi*k/n * x - atan(phi))+1;
    
    outputWaves(:,ii) = signalTest';
end

for ii = 1:4
waveTest2(:,ii) = unknownSystem2(outputWaves(:,ii));
end 

figure
subplot(2,1,1)
plot(x,outputWaves,'Linewidth',2)
title('Original Waves')
xlabel('time (s)')
ylabel('amplitude')
legend('Frequency 1', 'Frequency 2', 'Frequency 4', 'Frequency 8', 'Location', 'bestoutside')
box off


subplot(2,1,2)
plot(x, waveTest2,'Linewidth',2)
title('System 2 Transform')
xlabel('time (s)')
ylabel('amplitude')
legend('Frequency 1', 'Frequency 2', 'Frequency 4', 'Frequency 8', 'Location', 'bestoutside')
box off


% mean(waveTest1)
meanWave2 = mean(waveTest2) %mean is equal for all frequencies 

% The transformation of unknownSystem2 leaves the frequency intact but
% makes changes to the phase and amplitude with lower frequencies'
% amplitude being impacted more significantly than the higher frequences.
% The whole set seems to have had a constant added to it as well because
% the new waveforms are centered at 12.4896 instead of 1. 
% This is reflected by the mean spikes of each wave as they remains equal 
% to the others before and after the transformation, although the mean itself increases. 

%For graphing purposes removing the DC component we have from the +1 we
%added to a sinusoids to make them all positive. It makes it easier to view
%the component ampltudes. 
waveTest2 = waveTest2 - mean(waveTest2); 
outputWaves = outputWaves - mean(outputWaves);

freqTest1 = abs(fft(outputWaves));
freqTest2 = abs(fft(waveTest2));

figure
for ii = 1:4
    subplot(2,2,ii)
    stem(freqTest1(:,ii), 'Linewidth',1)
    hold on
    stem(freqTest2(:,ii), 'Linewidth', 1)
    title(['Transform of Sinusoid ' num2str(ii)])
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
    legend('Original Waves', 'Transformed Waves')
    box off
end

% This amplitude plot shows that the frequencies before and after the
% transformation are aligned with eachother, thus no change to the
% frequency occured. 

%% PART C

% If the previous tests succeeded, verify that the change in amplitude and phase from
% input to output is predicted by the amplitude (abs) and phase (angle) of the corresponding 
% terms of the Fourier transform of the impulse response. If not, explain which
% property (linearity, or shift-invariance, or both) seems to be missing from the system. If
% so, do you think that the combination of all tests guarantees that the system is linear and
% shift-invariant? What combination of tests would provide such a guarantee?


figure
for ii = 1:4

rVec(:,ii) = fft(waveTest2(:,ii));

lenWave = length(x);                         
fShift = (-lenWave/2:lenWave/2-1)*(50/lenWave);
yShift = fftshift(rVec(:,ii));

tol = 1e-6;
yShift(abs(yShift) < tol) = 0;

theta = angle(yShift);

fftAbs(:,ii) = abs(yShift);
testPhi(:,ii) = theta;

subplot(2,1,1)
stem(fShift,fftAbs(:,ii),'Linewidth', 1)
title('Amplitude Plot')
ylabel('Amplitude')
xlabel('Frequency (Hz)')
hold on

subplot(2,1,2)
stem(fShift,theta,'Linewidth', 1)
xlabel ('Frequency (Hz)')
ylabel ('Phase / \pi')
title('Phase Plot')
hold on

end

subplot(2,1,1)
legend('Frequency 1', 'Frequency 2', 'Frequency 4', 'Frequency 8', 'Location', 'bestoutside')

subplot(2,1,2)
legend('Frequency 1', 'Frequency 2', 'Frequency 4', 'Frequency 8', 'Location', 'bestoutside')



%Determining the amplitude and phase of the output based on the amplitude
%and phase of the impluse vector and the input waves mathmatically. 

fftImpulse2 = fft(sysTest(:,2));
fftImpulse2 = fftshift(fftImpulse2);

tol = 1e-6;
fftImpulse2(abs(fftImpulse2) < tol) = 0;
impPhase2 = angle(fftImpulse2);
imp2Amp = abs(fftImpulse2);

% Finding amplitude of impulse reponse * input
for ii = 1:4
testAmp2 = fft(outputWaves(:,ii));
testAmp2 = fftshift(testAmp2);

tol = 1e-6;
testAmp2(abs(testAmp2) < tol) = 0;

testAmp(:,ii) = abs(testAmp2) .* imp2Amp;

end 


%Finding phase of impulse response + input
testPhase = zeros(64,4);
for ii = 1:4
    testPhase2 = fft(outputWaves(:,ii));
    testPhase2 = fftshift(testPhase2);
    
    tol = 1e-6;
    testPhase2(abs(testPhase2) < tol) = 0;
    phi2(:,ii) = angle(testPhase2);
    for jj = 1:length(testPhase2)
        if testPhase2(jj) ~= 0 %phase is always present but we only care about when amplitude is more than 0
            testPhase(jj,ii) = phi2(jj,ii) +  impPhase2(jj,:);
        end
    end
end

%Testing amplitude comparisons
ampCheck = round(fftAbs,4) == round(testAmp,4)

%Testing phse comparisons
phiCheck = round(testPhi,4) == round(testPhase,4)

% We can see here though these comparisons that the amplitude and phase of
% the output can be predicted by combining the amplitude and phase of the
% impulse response to the system and the input waves. 

% The change in amplitude and phase from the input to the output is shown
% in the fourier transform of the impulse response. In the original wave
% form all the waves are of equal amplitude with varying phases, viewing
% the amplitude plot it is apparent that after the tranformation
% frequency 1's amplitude is increased the most, followed by frequency 2.
% Both the frequencies 4 and 8 are increased very slightly however they
% reamin relatively similar to eachother. 

% Reagarding the phase small shifts occur in all four frequencies however
% frequency 2 has the largest shift, as visible in the plot, followed by
% frequency 1, then 4 then 8 with the smallest change from the original.

% UnknownSystem3 failed the test of shift invariance because the output
% increased as the input time increased, this system also does not
% have perodic boundaries. 

% UnknownSystem1 failed the test of linearity because it does not hold the 
% property of aditivity and thus superposition, e.g. when to separate inputs
% are added together they do not equal the output of those two inputs run
% simultaniously. 

% UnknownSystem2 was the only system to pass both tests and to continue on. 

% Although system2 passed all of these tests it does not absolutly
% mean beyond a shadow of a doubt that the system is linear. The tests we
% did above were assuming resonable signals, and we have no way to measure
% every possible signal imaginalble or sounds so loud they would destroy
% the listener. Since we are measuring a discrete signal we cannot be 100%
% sure how the signal responds between those measured points. Additionally
% when testing the systems we only used impulse vectors however the system
% could be set up to handle inputs of particular amplitudes differently
% than others, or another variance we did not explcitly test. 
% When testing linearity above we did not test for homogeneity in addition 
% to additivity which could have lead to a false assumption about the linearity
% of systems 2 and 3. 

% There is no combination of tests that would provide a guarentee of true
% linearity, however taking into account the type of stimulus being used as
% an input we can reasonably assume that system 2 behehaves like it is
% linear within the standard input range. 

