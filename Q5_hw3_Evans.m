%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 3 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 5

% Date intialized: Oct. 15th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Sampling and aliasing. 
% Load the file myMeasurements.mat into matlab. It contains a
% vector, sig, containing voltage values measured from an EEG electrode, sampled at 100Hz.
% Plot sig as a function of vector time (time, in seconds), using the flag 'ko-' in matlab's plot
% command so you can see the samples.

load('myMeasurements.mat')
time = [0:119]/100;

figure
plot(time, sig, 'ko-')
xlabel('time (s)')
ylabel('voltage')
box off
set(gca, 'TickDir', 'out')
title('EEG Electrode Sampled at 100Hz')
hold on

%% PART A

%%
%*  This signal is only a small portion of the full data, and you don't want to store all those
% values. Create a subsampled version of the signal, which contains every fourth value.
% Plot this, against the corresponding entries of the time vector, on top of the original
% data (use matlab's hold function, and plot with flag 'r*-'). How does this reduced
% version of the data look, compared to the original? Does it provide a good summary
% of the original measurements? Is the subsampling operation linear? Shift-invariant?
% Explain.

indexVec1 = 4:4:120;

for ii = 1:length(indexVec1)
    jj = indexVec1(ii);
    subsampledSig1(ii) = sig(jj);
    subsampledTime1(ii) = time(jj);
end

plot(subsampledTime1, subsampledSig1, 'r*-', 'Linewidth', 2)
legend('Full Dataset', 'Subsampled') 
title('Comparison of Full Data and Subsampled Data')
box off
set(gca, 'TickDir', 'out')


% The reduced version of the data does not look like the original, it cuts
% off all the high voltage spikes and only stores the refractions where the
% volage dips.

% This subsampled system does not provide a good summary of the original
% because the subsampled structure is misleading, and underrepresents the data shape.  

%%
%* Test of Additivity 
% The subsampling opperation is likely linear because the properity of 
% additivity is upheld, as shown below in the plot. 
% The combindation of two outputs from different signals put through the
% same system is equal to two signals combined then put through the system.
% This shows additivity. 

indexVec2 = 4:4:120;
randSig = randn(1,120);

for ii = 1:length(indexVec2)
    jj = indexVec2(ii);
    subsampledRandSig(ii) = randSig(jj);
    subsampledSig2(ii) = sig(jj);
    subsampledTime2(ii) = time(jj);
end

combinedVec = [subsampledSig2'+ subsampledRandSig'];

figure
plot(subsampledTime2, combinedVec, 'r*-', 'Linewidth',2)
hold on

addedVec = [sig' + randSig'];

for ii = 1:length(indexVec2)
    jj = indexVec2(ii);
    subsampledSigA(ii) = addedVec(jj);
    subsampledTimeA(ii) = time(jj);
end

plot(subsampledTimeA, subsampledSigA, 'b*--','Linewidth',2)
legend('Combined Outputs', 'Output of Combined Inputs') 
title('Test of Linearity Additive Property')
xlabel('time (s)')
ylabel('voltage')
box off
set(gca, 'TickDir', 'out')

%%
%* Test of Homogeneity 
% The subsampling opperation is linear because it also obeys the rules of
% homogenity. The normal output multiplied by a constant after it goes though
% the system is equal to the output of the input multiplied by the same
% constant before it goes through the subsampling system. 


indexVec3 = 4:4:120;

sig5 = sig*5;
time5 = time*5;

for ii = 1:length(indexVec3)
    jj = indexVec3(ii);
    subsampledSigH(ii) = sig(jj);
    subsampledTimeH(ii) = time(jj);
end

for ii = 1:length(indexVec3)
    jj = indexVec3(ii);
    subsampledSig5(ii) = sig5(jj);
    subsampledTime5(ii) = time5(jj);
end


figure
plot(subsampledTimeH*5, subsampledSigH*5, 'r*-', 'Linewidth', 2)
hold on
plot(subsampledTime5, subsampledSig5, 'b*--', 'Linewidth', 2)
legend('Original Output *5', 'Input* 5') 
title('Test of Homogeneity')
xlabel('time (s)')
ylabel('voltage')
box off
set(gca, 'TickDir', 'out')


%%
%* Test of Shift Invariance 
% The subsampling opperation is not shift invariant, if you start at a
% different point the outcome will be different (ie if you start at 1
% instead of 4 (so each point selected is shifted by 1) the output uses all
% high voltage points instead of low voltage points, even when still counting by 4). 
% The plot below shows the outputs of the shifted inputs 

indexVec2 = 1:4:120;

for ii = 1:length(indexVec2)
    jj = indexVec2(ii);
    subsampledSig3(ii) = sig(jj);
    subsampledTime3(ii) = time(jj);
end

figure
plot(subsampledTime1, subsampledSig1, 'r*-', 'Linewidth', 2)
hold on
plot(subsampledTime3, subsampledSig3, 'b*-', 'Linewidth', 2)
legend('Input 1', 'Input 2') 
title('Test of Shift Invariance')
xlabel('time (s)')
ylabel('voltage')
box off
set(gca, 'TickDir', 'out')

%% PART B

%%
%* Examine your EEG result in the frequency domain. First plot the magnitude (amplitude) 
% of the Fourier transform of the original signal, over the range [-N/2,(N/2) - 1]
% (use fftshift). By convention, the 'Delta' band corresponds to frequencies less than
% 4Hz, 'Theta' band is 4-7Hz, 'Alpha' band 8-15Hz, and 'Beta' is 16-31Hz. For these
% data, which band shows the strongest signal?


fftEEG = fft(sig);

lenWave = length(sig);   
sampFreq = 100; %Hz
fShift1 = (-lenWave/2:lenWave/2-1)*sampFreq/lenWave;
yShift1 = fftshift(fftEEG);

figure
stem(fShift1,abs(yShift1),'Linewidth',1)
ylabel('Amplitude')
xlabel('Frequency (Hz)')
title('Amplitude of EEG in Frequency Domain')
box off
set(gca, 'TickDir', 'out')

% The stronget frequency signal response is 30Hz which falls in the Beta range, the
% second highest is 10Hz which is part of the Alpha range. 


%% PART C

%%
%* Plot the Fourier magnitude for signals downsampled by factors of 2, 3, and 4, after
% upsampling them back to full size (i.e., make a full-size signal filled with zeros, and
% set every kth sample equal to one of the subsampled values, for subsampling by factor
% k). What is the relationship between these plots and the original frequency plot. What
% has happened to the frequency components of the original signal? Does the strongest
% signal band change?

resampledOutput = zeros(120,3);

sampleVec = 2:4;
figure
subplot(2,2,1)
stem(fShift1,abs(yShift1),'Linewidth',1)
title('Full Data Set')
ylabel('Amplitude')
xlabel('Frequency (Hz)')
box off
set(gca, 'TickDir', 'out')

for ii = 1:3
    pp = sampleVec(ii);
    
    indexVec = pp:pp:120;

for qq = 1:length(indexVec)
    hh = indexVec(qq);
    subsampledSig(qq) = sig(hh);
  
end

resampledOutput(indexVec,ii) = subsampledSig;

fftEEG = fft(resampledOutput(:,ii));

lenWave = length(sig);                         
fShift = (-lenWave/2:lenWave/2-1)*sampFreq/lenWave;
yShift = fftshift(fftEEG);

subplot(2,2,ii+1)
stem(fShift,abs(yShift),'Linewidth',1)
ylabel('Amplitude')
xlabel('Frequency (Hz)')
title(['Sampling Rate of ' num2str(pp)])
box off
set(gca, 'TickDir', 'out')
axis([-50,50,0,150])

clear subsampledSig;
end

% When downsampling rate increases the spread of frequencies becomes greater.
% The original highest amplitude frequencies (10 and 30) are still present
% at all sampling rates, and appear with a stable proportion to one another.

% The overall maximum amplitude drops as the downsampling rate increases,
% however more frequencies are represented, splitting the original total
% amplitude between them. 

% Because of this as the downsampling rate increases more activity is seen
% across the different signal bands with equal amplitudes being represented
% in Delta, Theta, Alpha, Beta and Gamma. The signal appears to repeat with
% more cycles as the downsampling rate increases. 

% Each time the sampling rate is halved the frequencies are folded across
% eachother adding mirrored allias of equiviolent amplitude. This is
% because the sampling rate is a lower frequency than the original signal
% so it is picking up values outside each discrete wavelengths. Since there
% is an overlap between the groups of the copies and the original frequencies 
% it will be more difficult to filter them out to return to k=0 unless the
% Niquist criteron is satisifed. 


