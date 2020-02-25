%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 2 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 4

% Date intialized: Sept. 27th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Dimensionality reduction with PCA.
% Professors Hugh Bell and Wee Zell were recording
% extracellular action potentials (i.e. spikes) from cat primary visual cortex late one evening
% when their computer malfunctioned. It had already isolated a set of 400 time windows in
% which voltages had crossed a threshold, indicating the presence of spike. But these traces
% likely arose from multiple cells, with each cell producing a characteristic waveform, and the
% computer failed before sorting the voltage traces to determine how many cells were present,
% and which spikes arose from each cell. The professors come to you (the only 
% math-toolsenabled graduate student still in the building at that hour), asking for help. They provide
% you with a file windowedSpikes.mat containing a 400 × 150 matrix, data, whose rows
% contain the electrode measurements (voltages recorded for each 150 msec window, at 1msec
% intervals). Your task is to determine how many neurons produced these 400 spikes.
load('windowedSpikes.mat');

%%
%* (a) Plot the 400 waveforms superimposed and describe what you see. Be sure to label your
% axes! Using these spike waveform plots, can you devise a way to deduce how many
% neurons produced these spikes? Feel free to include an additional plot containing just
% a subset of the waveforms in order to aid in your explanation.

figure;
plot(data')
title('All Neural Responses')
xlabel('Time (m/s)')
ylabel('Voltage Level')
box off
set(gca, 'TickDir', 'out')

disp('The neurons have the largest activity between 53 and 85 m/s into the')
disp('recording. This will thus be the time point which shows the highest definition')
disp('between neurons')

disp('When randomly plotting 30x different neural signals out of the 400 there appears to be 3x')
disp('distinct patterns that emerge, one with a strong negative dip, one with a prominate postive peak')
disp('and a third with jagged low level persistent activity during the recording time')

disp('There is noise present in these patterns so no two neurons are identical')
disp('however the three distinctive shapes are seen throughout the data, indicating that the')
disp('recording was primarily from three princial sources')

%plot a sample of the neuron with the negative dip
figure;
plot(data(244,:))
title('Sample Signal Response 1 - large negative dip (Sample 244)')
xlabel('time m/s')
ylabel('Voltage Level')
box off
set(gca, 'TickDir', 'out')
axis([-0 150 -100 100])

%plot a sample of the neuron with a prominant positive peak
figure;
plot(data(319,:))
title('Sample Signal Response 2 - large positive peak (Sample 319)')
xlabel('time m/s')
ylabel('Voltage Level')
box off
set(gca, 'TickDir', 'out')
axis([-0 150 -100 100])

%plot a sample of the neuron with jagged low level persistent activity during the recording time
figure;
plot(data(107,:))
title('Sample Signal Response 3 - Persistent Activity (Sample 107)')
xlabel('time m/s')
ylabel('Voltage Level')
box off
set(gca, 'TickDir', 'out')
axis([-0 150 -100 100])

%plotting the set of 30 random example neurons:
r = randi([1 400],30,1);
figure;
sgtitle('Sample 30x Neural Responses')
for ii = 1:30
    subplot(10,3,ii)
    plot(data(r(ii,:),:))
    box off
set(gca, 'TickDir', 'out')
axis([-0 150 -100 100])
end  


%%
%* (b) Perform principal components analysis (PCA) on your data, and plot the eigenvalues
% in descending order (alternatively, compute the SVD of data). It might help to display
% the eigenvalues on a log-scale. Interpret what you see.

% using SVD
[U S V] = svd(data);

S2 = diag(S);
M = V*S'*S*V';

eigenVals = S2.^2;
eigenVectors = V';
logEigen = log(eigenVals);

%Plot the log-scale eigenvalues 
figure;
bar(logEigen)
title('Eigen Values')
box off
set(gca, 'TickDir', 'out')
xlabel('Eigen Index #')
ylabel('EigenValue (log)')

disp('There are 3x neurons which are responsible for the majority of the spikes')
disp('according to the eigen values. This is displayed by the 3x majority')
disp('columns all the way to the left of the function. Because the graph shows')
disp('the log function of the eigen values the disparity between the top 3x')
disp('values and the rest of the group looks less exteme')

%%
%* (c) Project each of the 400 spike waveforms onto the top two principal components of the
% dataset, and plot the resulting values as points in 2 dimensions. Describe what you see.
% Can you deduce how many distinct neurons produced the 400 voltage traces?

principal1 = data*eigenVectors(1,:)';
principal2 = data*eigenVectors(2,:)';

%2D plot of data on pricipal components
figure;
plot(principal1,principal2,'.')
title('Scree Plot')
box off
set(gca, 'TickDir', 'out')
xlabel('PC1')
ylabel('PC2')

   
disp('By looking at the scatter plot of the 400 waveforms on the two principal')
disp('comoments we see three distinct clusters, this supports the eigenvalue')
disp('plot above showing only 3x primary contributing factors.')

%%
%* (d) Now project each waveform onto the top three principal axes, and plot in 3 dimensions
% (you may want to spin it around, using rotate3d in matlab). Are there any significant
% changes you see? Using the 3D plot, can you inform Drs. Bell and Zell how many
% neurons they likely recorded from?

principal1 = data*eigenVectors(1,:)';
principal2 = data*eigenVectors(2,:)';
principal3 = data*eigenVectors(3,:)';

%3D plot of data on pricipal components
figure;
plot3(principal1,principal2,principal3,'.')
title('3D Scree Plot')
box off
set(gca, 'TickDir', 'out')
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
rotate3d on
grid on

disp('In the 3D plot there appears to be 4x primary clusters, however the 4th')
disp('cluster that appears is much less defined than the three seen originally')
disp('in the 2D space. This fourth cluster could account for the noise found in')
disp('the data or a less influencial nearby neuron. ')

%%
%* (e) Extra credit: how would you determine which neuron fired each of the 400 observed
% spikes? Describe your strategy, including appropriate linear algebraic expressions for
% the required calculations.

disp('When taking the dot product of two vectors that are the same, the outcome')
disp('will be larger than for two vectors that are disimilar.')


disp('By taking the dot product of each neural response with each of the three eigen vectors representing ')
disp('the most influencial signals, the pairing with the highest dot product')
disp('will be the most similar. We can then tag this neural response as being')
disp('caused by neuron 1, 2 or 3 based on which combindation resulted in the highest outcome.')
disp('First I"ve normalized the data and made all vectors of the neural responses unit')
disp('length to compare to the unit length eigen vectors.')
disp('Now that the lengths are the same the dot product will be mainly')
disp('influenced by the angle the vector is pointing. We are able to tell which')
disp('neuron was responding in the observed time period by looking at the')
disp('results of these dot products. ')

D = data;

%Standardize data
stData = (D - mean(D,1)) ./std(D,1);

%Make all data vectors unit length
dNorm = stData/norm(stData);
[Un Sn Vn] = svd(dNorm);

eigenVectors = Vn;

sampleout = nan*ones(400,1);

for ii = 1:400
    
    %take the dot product of the data with each eigen vector 
    test(1) = stData(ii,:)*eigenVectors(:,1);
    test(2) = stData(ii,:)*eigenVectors(:,2);
    test(3) = stData(ii,:)*eigenVectors(:,3);
    
    %find the largest dot product
    vec = find(test == max(test));
    
    %place index in matrix
    sampleout(ii) = vec;
    
end

figure;
histogram(sampleout)
title('Approximate date score by neural response')
ylabel('Number of Neural Responses')
xlabel('Responsible Cell (1, 2 or 3)')
set(gca, 'XTick', [1, 2, 3])




