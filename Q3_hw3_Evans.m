%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 3 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 3

% Date intialized: Oct. 15th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all
close all
clc

%% Gabor filter.
% The response properties of neurons in primary visual cortex (area V1) are often
% described using linear filters. We?ll examine a one-dimensional cross-section of the most
% common choice, known as a 'Gabor filter' (named after Electrical Engineer/Physicist Denis
% Gabor, who developed it in 1946 for use in signal processing).

%% PART A

%%
%* Create a one-dimensional linear filter that is a product of a Gaussian
% and a sinusoid. The filter should contain 25 samples, and the Gaussian should be centered on the middle (13th)
% sample. Plot the filter to verify that it looks like what you'd expect. Plot the amplitude
% of the Fourier transform of this filter, sampled at 64 locations (MATLAB's fft function
% takes an optional additional argument). What kind of filter is this? Why does it have
% this shape, and how is the shape related to the choice of parameters (sigma, w)? Specifically,
% how does the Fourier amplitude change if you alter each of these parameters?

sigma = 3.5;
w = 2*pi*10/64;
mu = 13;
x = (0:1:24);  % Plotting range

gaborFilt = exp(-((x-mu).^2/(2*sigma^2))) .*cos(w*x);

figure
subplot(2,1,1)
plot(x, gaborFilt,'Linewidth',2)
xlabel('Sample #')
ylabel('Signal Amplitude')
title('Visualization of Gabor Filter Freq. 10')
box off
set(gca, 'TickDir', 'out')

fourierG = fft(gaborFilt, 64);

lenWave = length(fourierG);
fShift = (-lenWave/2:lenWave/2-1);
yShift = fftshift(fourierG);

subplot(2,1,2)
stem(fShift, abs(yShift),'Linewidth',1)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Amplitude of Frequency Componets')
box off
set(gca, 'TickDir', 'out')

% This is a band pass filter. It has this shape because it is selective to
% frequencies near the original input of the sinusoid (10 in this example),
% and not to frequencies higher or lower than it. 

%%
%* Looking at the DTF of the components themselves we can see that when the 
% gaussian and the sinousiod are combined the structure of the filter
% remains similar to the fourier components of the two inputs. 

x = -32:1:31;
c = cos(w*x);
g = exp(-((x-mu).^2/(2*sigma^2)));

fftC = fft(c);
fftG = fft(g);

cShift = fftshift(fftC);
gShift = fftshift(fftG);

figure
stem(x,abs(cShift), 'Linewidth',1)
hold on
stem(x,abs(gShift),'Linewidth',1)
legend('Cosine Component','Gaussian Component')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('DFT of Separate Components of Filter')
box off
set(gca, 'TickDir', 'out')

%%
%* As sigma is incerased from zero the signal begins to take shape around 3,
% then as it is increased more beyond that the noise starts to spread. When
% graphed the lower the sigma the more filter resembles the gaussian
% component of the filter and as sigma increases the more input it takes
% from the sinousoidal component. 

sigSet = [1,3,5,7,9];
x = (0:1:24);
figure
for ii = 1:5
    jj = sigSet(ii);
    w = 2*pi*10/64;
    gaborFilt1 = exp(-((x-mu).^2/(2*jj^2))) .*cos(w*x);
    fourierTest = fft(gaborFilt1, 64);
    
    lenWave = length(fourierTest);
    fShift = (-lenWave/2:lenWave/2-1);
    yShift = fftshift(fourierTest);
    subplot(5,1,ii)
    stem(fShift, abs(yShift), 'Linewidth',1)
    title(['Amplitude of Sigma ' num2str(jj)])
    box off
    set(gca, 'TickDir', 'out')
    ylim([0 10])
end

%%
%* When w is close to 1 the two high amplitude frequencies are close
% together, as the frequency of the sinousoid is increased the two primary
% frequency components of the filter move apart until they are
% furthest at 32 (64/2), they then start moving closer together again until
% they overlap at 64 since the second half of the spectrum is a mirror of
% the first half, just with the signs flipped on the imaginary components. 
% Each point on the spectrum is covered while moving from 1 to 32.

freqSet = [0, 8, 16, 24, 32];
figure
for ii = 1:5
    jj = freqSet(ii);
    w = 2*pi*jj/64;
    gaborFilt2 = exp(-((x-mu).^2/(2*sigma^2))) .*cos(w*x);
    fourierTest = fft(gaborFilt2, 64);
    
    lenWave = length(fourierTest);
    fShift = (-lenWave/2:lenWave/2-1);
    yShift = fftshift(fourierTest);
    subplot(5,1,ii)
    stem(fShift, abs(yShift), 'Linewidth',1)
    title(['Amplitude of Frequency ' num2str(jj)])
    box off
    set(gca, 'TickDir', 'out')
    ylim([0 10])
end
%% PART B

%%
%*  If you were to convolve this filter with sinusoids of different frequencies, which of them
% would produce a response with the largest amplitude? Obtain this answer by reasoning
% about the equation defining the filter (above), and also by finding the maximum of the
% computed Fourier amplitudes (using the max function), and verify that the answers
% are the same. 

% Based on the fourier transform of the filter, amplitude is highest for
% frequencies matching the sinousoid used to make the filter, a frequency
% of 10. This will not be compressed or scaled as it's full input will be
% let through the filter. Any frequencies more or less than 10 will be
% compressed when the filter is applied by the amplutude present at that
% frequency on the filter. 

maxAmp = max(abs(fourierG))
maxAmpFreq = find(abs(fourierG) == max(abs(fourierG)))

%Matlab positon 11 = frequency 10

%%
%* Compute the period of this sinusoid, measured in units of sample spacing
% (hint: this is the inverse of its frequency, in cycles/sample), and verify by eye that this is
% roughly matched to the oscillations in the graph of the filter itself. 
nSamp = 64;
x = 0:1:24;
k = 10;

wavePeriod = (1/k)*nSamp;
index = 0:wavePeriod:30;
figure
plot(x,gaborFilt,'Linewidth',2)
hold on
plot(index,0,'ro','Linewidth',2)
xlabel('Sample #')
ylabel('Amplitude')
title('Period of signal in cycles/sample')
legend ('Gabor Filter', 'f10 Period','Location', 'bestoutside')
box off
set(gca, 'TickDir', 'out')

% The spacing of the wave period is visually congruent with the occilations
% in the filter itself. 

%%
%* Which two sinusoids would produce responses with about 25% of this maximal amplitude?

% The frequencies which would elicit an response with an amplitude one quarter
% the size of the maximum would be 5 and 15 because looking at the
% amplitude plot of the filter those points have an amplitude of
% approximatly 1.1 which is 25% of the maximum amplitude of 4.4. 

quarterAmp = maxAmp/4


%% PART C

%%
%*  Create three unit-amplitude 64-sample sinusoidal signals at the three frequencies (low,
% mid, high) that you found in part (b). Convolve the filter with each, and verify that
% the amplitude of the response is approximately consistent with the answers you gave
% in part (b). (hint: to estimate amplitude, you'll either need to project the response onto
% sine and cosine of the appropriate frequency, or compute the DFT of the response and
% measure the ampitude at the appropriate frequency).

mu = 13;
sigma = 3.5;
w = 2*pi*10/64;
x = (0:1:24);
gaborFilt = exp(-((x-mu).^2/(2*sigma^2))) .*cos(w*x);

nSamp = 64;
x = 0:1:63;
xShift = (-nSamp/2:1:nSamp/2-1);
highFreq = 15;
medFreq = 10;
lowFreq = 5;

highWave = cos(2*pi*highFreq/nSamp * x);
medWave = cos(2*pi*medFreq/nSamp * x);
lowWave = cos(2*pi*lowFreq/nSamp * x);

convHigh = conv(highWave, gaborFilt, 'same');
convMed = conv(medWave, gaborFilt,'same');
convLow = conv(lowWave, gaborFilt,'same');

fftHigh = fft(convHigh);
fftMed = fft(convMed);
fftLow = fft(convLow);

figure
subplot(3,1,1)
stem(xShift,(fftshift(abs(fftHigh))), 'Linewidth',1)
ylim([0 150])
ylabel('Amplitude')
xlabel('Frequency (Hz)')
title('Frequency 15')
box off
set(gca, 'TickDir', 'out')

subplot(3,1,2)
stem(xShift,(fftshift(abs(fftMed))), 'Linewidth',1)
ylim([0 150])
ylabel('Amplitude')
xlabel('Frequency (Hz)')
title('Frequency 10')
box off
set(gca, 'TickDir', 'out')

subplot(3,1,3)
stem(xShift,(fftshift(abs(fftLow))), 'Linewidth',1)
ylim([0 150])
ylabel('Amplitude')
xlabel('Frequency (Hz)')
title('Frequency 5')
box off
set(gca, 'TickDir', 'out')


maxHigh = max(abs(fftHigh));
maxMed = max(abs(fftMed));
maxLow = max(abs(fftLow));

percentLow = 100*(maxLow/maxMed)
percentHigh = 100*(maxHigh/maxMed)

% The selected frequencies were approximatly at the amplitude percentages
% estimated in part B. The difference is due to sampling, since the signal
% is not continuous its likely that we did not sample a point specifically
% at 25%, however our chosen frequencies of 5 and 15 are the two closest
% points to 25% out of the options we have available. 
