%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% VISUAL CHANGE DETECTION DATA ANALYSIS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Date Created: July 3rd 2019
% Date Last Edit: July 9rd 2019
% Author: Marissa Evans - mhe229@nyu.edu

% DATA DISCLOSURE: All data presented here is the property of Weiji Ma's
% lab and cannot be duplicated or analyzed without explicit permission. 

% Assumptions:
% The 'Data Structure Organizer' script must be run BEFORE this script in
% order for it to function.
% Files 'allParticipantData.m' and 'outputAcrossAprticipants.m' are in the
% same directory as the script.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment description: This data is from a visual change detection
% experiment. During the experiment participants view a display with 4x
% ellipses which either have a rounder shape (low reliability) or a narrow
% shape (high reliability). It is randomly decided which of the ellipses
% are round and which are narrow on a given trail. After a short delay the
% display changes to either 4x ellipses (same shape and orientation as the
% previous screen) or 4x straight lings (on the same orientation as the
% previous screen). On 50% of the trails 1 of the 4 ellipses changes it
% orientation from 1-90 degrees. The role of the participant is to respond
% with a key press if YES there was a change between the two screens, or NO
% the two screens were identical. Each participant completes 1000 trials in
% EACH condition. Within each condition we split the trails into 5x groups
% depending on how many NARROW (high reliability) ellipses were present on
% the screen in a given trial (0,1,2,3,or 4). Due to the fully randomized
% presentation of the stimuli there are an uneven number of trails for each
% participant in each reliability level.

% For this analysis we are mainly looking at the lowest reliability level (0
% narrow ellipses) and the highest relaibiltiy level (4 narrow ellipses) to
% compare the difference in performance between them, with regards to the
% stimulus condition (either ellipses or lines presented on the 2nd
% display).

% Additionally in this analysis we are only working with the metric of
% 'Correct' vs '# of High Reliability Ellipses'. In the data we are also
% planning to do similar analysis on the 'Probability Respond Change' vs
% 'Amount of Rotation'.

% This experiment was orignally designed to be analysed using bayesian
% model testing so fewer predictors were recorded outside of that. 


% What this code does(by section):

% Section 0 - Initalization
% Section 1 - Loader
% Section 2 - Pulling notable data from structure
% Section 3 - Graphing the distributions across conditions
% Section 4 - Classical Statistics
% Section 5 - Bootstrapping data for both condtions to check mean stability
% Section 6 - PCA
% Section 7 - Permutation tests (ellipse high vs low, line high vs low,
%             ellipse high vs line high)
% Section 8 - SVA and Random Forest Models 
% Section 9 - Conclusions based on data analysis

% OUTPUTS: 4x figures and an ANOVA table- *PRESS ANY KEY TO VIEW NEXT
% FIGURE (pause is present between each figure)*


%% 0 INITALIZATION
clear all %Clear workspace
close all %Close all open figures
clc %clear command window

%% 1 LOADER

% Reminder- these files must be created in the 'RUN_FIRST' script, without
% running that first they will not be available to open. 

load('allParticipantData'); %load the structure made in the first script
load('outputAcrossAprticipants'); %load the other structure made in the first script

fS = 20; %Sets font size for figures
lW = 3; %Sets the line width for figures

%% 2 DATA ORGANIZATION AND MATRIX BUILDING

%pulling out the matricies of interest for analysis- the lowest and highest
%reliability levels for the comparison between them is of the most interest
%to us. These are output for the participant response for each trail,
%accross all participants, for each level of reliability. 0 = incorrect
%reponse, 1 = correct response. Columns 2 and 3 read 1 or
%ELLIPSE or 2 or LINE condition, and in the 3rd column is the number of
%high rel (0:4)

%Expaning the ELLIPSE matricies out of the structure
EZeroRelCorrect = outputAccrossSubj.Ellipse.Correct.zeroRel;
EOneRelCorrect = outputAccrossSubj.Ellipse.Correct.oneRel;
ETwoRelCorrect = outputAccrossSubj.Ellipse.Correct.twoRel;
EThreeRelCorrect = outputAccrossSubj.Ellipse.Correct.threeRel;
EFourRelCorrect = outputAccrossSubj.Ellipse.Correct.fourRel;

%Expanding the LINE matricies out of the structure 
LZeroRelCorrect = outputAccrossSubj.Line.Correct.zeroRel;
LOneRelCorrect = outputAccrossSubj.Line.Correct.zeroRel;
LTwoRelCorrect = outputAccrossSubj.Line.Correct.zeroRel;
LThreeRelCorrect = outputAccrossSubj.Line.Correct.zeroRel;
LFourRelCorrect = outputAccrossSubj.Line.Correct.fourRel;

%Concatonating all trials for all participants and conditions 
AllTrialsMat = [outputAccrossSubj.Line.AllRel; outputAccrossSubj.Ellipse.AllRel;];

%Ellipse Condition Means & Standard Error
meansEllipse(1,1) = mean(EZeroRelCorrect(:,1));
meansEllipse(2,1) = mean(EOneRelCorrect(:,1));
meansEllipse(3,1) = mean(ETwoRelCorrect(:,1));
meansEllipse(4,1) = mean(EThreeRelCorrect(:,1));
meansEllipse(5,1) = mean(EFourRelCorrect(:,1));

SEMEllipse(1,1) = std(EZeroRelCorrect(:,1))/sqrt(length(EZeroRelCorrect(:,1)));
SEMEllipse(2,1) = std(EOneRelCorrect(:,1))/sqrt(length(EOneRelCorrect(:,1)));
SEMEllipse(3,1) = std(ETwoRelCorrect(:,1))/sqrt(length(ETwoRelCorrect(:,1)));
SEMEllipse(4,1) = std(EThreeRelCorrect(:,1))/sqrt(length(EThreeRelCorrect(:,1)));
SEMEllipse(5,1) = std(EFourRelCorrect(:,1))/sqrt(length(EFourRelCorrect(:,1)));


%Line Condition Means & Standard Error
meansLine(1,1) = mean(LZeroRelCorrect(:,1));
meansLine(2,1) = mean(LOneRelCorrect(:,1));
meansLine(3,1) = mean(LTwoRelCorrect(:,1));
meansLine(4,1) = mean(LThreeRelCorrect(:,1));
meansLine(5,1) = mean(LFourRelCorrect(:,1));

SEMLine(1,1) = std(LZeroRelCorrect(:,1))/sqrt(length(LZeroRelCorrect(:,1)));
SEMLine(2,1) = std(LOneRelCorrect(:,1))/sqrt(length(LOneRelCorrect(:,1)));
SEMLine(3,1) = std(LTwoRelCorrect(:,1))/sqrt(length(LTwoRelCorrect(:,1)));
SEMLine(4,1) = std(LThreeRelCorrect(:,1))/sqrt(length(LThreeRelCorrect(:,1)));
SEMLine(5,1) = std(LFourRelCorrect(:,1))/sqrt(length(LFourRelCorrect(:,1)));



%% 3 VISUALIZE THE DATA ACROSS CONDITIONS

%Create a figure and add the data
meanFig = figure %open a figure
relNum = 0:4; %set the number of high reliability items to 0,1,2,3,4
ellipseError = errorbar(relNum, meansEllipse, SEMEllipse); %graph the means against the # high rel ellipses with error bars
hold on %wait to add another line
lineError = errorbar(relNum, meansLine, SEMLine); %graph the means against the # high rel ellipses with error bars

%Atributes of the figure
title('Average Correct Response by Reliability Level'); %figure title
ylabel('Mean Correct Response'); %y-axis label
xlabel('Number of High Reliability Ellipses'); %x-axis label
set(gcf, 'Position', [100 100 800 800]); %figure position on screen
set(gca, 'XLim', [-.5  4.5]); %range of x axis
set(gca, 'XTick', [0 1 2 3 4]); %ticks on x axis
set(gca, 'YLim', [.5 .8]); %range of y axis
set(gca, 'YTIck', [.5 .55 .6 .65 .7 .75 .8]); %ticks on y axis
ellipseError.LineWidth = lW; %set the line width for Ellipse
lineError.LineWidth = lW; %set the line width for Line
set(gca,'TickDir', 'out'); %put the ticks outside the plot
set(gca, 'Color', 'w'); set(meanFig, 'Color', 'w'); %set plot and background color to white
set(gca, 'FontSize', fS); %set fontsize
box off %turn off the box around the plot

%Legend
legend('Ellipse Condition', 'Line Condition'); %add a legend
legend('Location', 'SouthEastOutside') %place legend outside the plot to the lower right

shg

pause %wait for button press to continue



%% 4 CLASSICAL STATISTICS

% T-test comparing the lowest number of high reliability stimulus to the
% highest number of high reliability stimulus (ellipse condition)
[sigE, pE, CIE, tstatsE] = ttest2(EZeroRelCorrect, EFourRelCorrect)
%tstat: -5.2593, df: 2250, sd: 0.45594, p: 1.5825e-07 - SIGNIFICANTLY DIFFERENT AT .01

% T-test comparing the lowest number of high reliability stimulus to the
% highest number of high reliability stimulus (line condition)
[sigL, pL, CIL, tstatsL] = ttest2(LZeroRelCorrect, LFourRelCorrect)
%tstat: -2.0122, df: 2267, sd: 0.47533, p: 0.044317 - NOT SIGNIFICANTLY
%DIFFERENT AT .01

% Multifactor ANOVA on the correct answer grouped both by condition and level of high
%reliability ellipses.
[p, anovatab] = anovan(AllTrialsMat(:,1),[AllTrialsMat(:,2),AllTrialsMat(:,3)], 'model','full',...
    'varnames',{'Stim Category','Num High Rel'})
% Significant difference in correct responses between stimulus conditions (F=38.15) and
% reliability levels (F=30.35), however there is NOT a significant
% interaction between stim condition and reliability level (F=1.38). 



%% 5 BOOTSTRAPPING THE DATA FOR BOTH CONDITIONS WITH HIGEST/LOWEST RELIABILITY

% Comparing the LOWEST reliability level and the HIGHEST relaibiltiy level
% for BOTH the ellipse and line conditions.



%% 5a ELLIPSE CONDITION 0 vs 4 NARROW ELLIPSES

eEmpiricalMeanDiff = mean(EZeroRelCorrect) - mean(EFourRelCorrect) %find the emperical difference between the ellipse means
%eEmpiricalMeanDiff = .1

numSamples = 1e5; %100k times
E_resampledMeans = nan(numSamples,2); %preallocate for all resamples means
n1 = length(EZeroRelCorrect); %number of samples equal to zero rel trials
n3 = length(EFourRelCorrect); %number of samples equal to four rel trials

for ii = 1:numSamples
    index = randi(n1, [n1, 1]); %we will use this index to subsample
    resampledEZeroRatings = EZeroRelCorrect(index); %pull ratings out of the zero rel matrix as chosen by index
    E_resampledMeans(ii,1) = mean(resampledEZeroRatings); %put means of resampled zero rel 'correct' scores into matrix
    
    index = randi(n3, [n3, 1]); %we will use this index to subsample
    resampledEFourRatings = EFourRelCorrect(index); %pull ratings out of the four rel matrix as chosen by index
    E_resampledMeans(ii,2) = mean(resampledEFourRatings); %put means of resampled four rel 'correct' scores into matrix
end

E_meanDifferences = E_resampledMeans(:,1) - E_resampledMeans(:,2); %Calculate the resampled mean difference for all 100k runs

% PLOTTING THE NULL DISTRIBUTION 
figure %open a figure
subplot(2,1,1); %add a subplot
histogram(E_meanDifferences,100); %make a histrogram of all mean difference for the ellipse condition scores in bins of 100
set(gca, 'Color', 'w'); set(gcf, 'Color', 'w'); %set figure and background to whtie
title('distribution of resampled Ellipse mean differences'); %title 
xlabel('Difference of resampled means'); %x axis label
ylabel('#'); %y axis label
set(gcf, 'Position', [100 100 800 800]); %position figure
set(gca, 'FontSize', fS); %set font size
set(gca, 'TickDir', 'out') %set tick direction out
set(gca, 'XLim', [-.18 .04])
box off %turn the outer axes off

% ADDING THE EMPERICAL MEAN
h= line([eEmpiricalMeanDiff eEmpiricalMeanDiff], [min(ylim) max(ylim)]); %make a vertical line in the plot from the top to the bottom
set(h,'Color', 'r'); set(h, 'LineWidth',lW); %set the line color and width

% ADDING CONDIDENCE INTERVALS 
sortedResampledMeans = sort(E_meanDifferences); %sort the mean differences in assending order
CIwidth = 99; %using a 99% CI
lowerBound = (100-CIwidth)/2; %calculate the lower bound of the CI
upperBound = 100-lowerBound; %calculate the upper bound of the CI

%find corresponding indices
lowerBoundIndex = round(lowerBound*length(sortedResampledMeans)/100); %Find the placement for the lower bound
upperBoundIndex = round(upperBound*length(sortedResampledMeans)/100); %Find the placement for the upper bound

% reach into the sorted mean differences to find the value of of the upper
% and lower CI
lowerCI = sortedResampledMeans(lowerBoundIndex); %set the lower CI value 
upperCI = sortedResampledMeans(upperBoundIndex); %set the upper CI value

% add the CI lines to the figure
gca
lCI = line([lowerCI lowerCI], [min(ylim) max(ylim)]); set(lCI, 'color','g'); set(lCI, 'LineWidth', 3);
uCI = line([upperCI upperCI], [min(ylim) max(ylim)]); set(uCI, 'color','g'); set(uCI, 'LineWidth', 3);

%Adding a legend
legend('Null Dist','Test Statistic', 'Lower 99% CI', 'Higher 99% CI'); %add a legend
legend('Location', 'SouthEastOutside'); %place legend outside the plot to the lower right

% ANALYSIS: Because the emperical mean fit squarely between the two
% confidence intervals we can assume our mean differnece is stable and not
% likely reached by chance. 



%% 5b LINE CONDITION 0 vs 4 NARROW ELLIPSES

lEmpiricalMeanDiff = mean(LZeroRelCorrect) - mean(LFourRelCorrect) %calculate the emperical mean difference
%lEmpiricalMeanDiff = .04

numSamples = 1e5; %100k times
L_resampledMeans = nan(numSamples,2); %preallocate for all resamples means
n1 = length(LZeroRelCorrect); %number of samples equal to zero rel trials
n3 = length(LFourRelCorrect); %number of samples equal to four rel trials

for ii = 1:numSamples
    index = randi(n1, [n1, 1]); %we will use this index to subsample from the zero narrow ellipses correct score
    resampledLZeroRatings = LZeroRelCorrect(index); %pull ratings out of the zero rel matrix as chosen by index
    L_resampledMeans(ii,1) = mean(resampledLZeroRatings); %put means of resampled zero rel 'correct' scores into matrix
    
    index = randi(n3, [n3, 1]); %we will use this index to subsample from the four narrow ellipses correct score
    resampledLFourRatings = LFourRelCorrect(index); %pull ratings out of the four rel matrix as chosen by index
    L_resampledMeans(ii,2) = mean(resampledLFourRatings); %put means of resampled four rel 'correct' scores into matrix
end

L_meanDifferences = L_resampledMeans(:,1) - L_resampledMeans(:,2); %Calculate the resampled mean difference for all 100k runs

% PLOTTING THE NULL DISTRIBUTION 
subplot(2,1,2); %add to the existing figure in a lower subplot
histogram(L_meanDifferences,100); %plot a histogram of all the mean differences for the line condition, sorting into bins of 100

title('distribution of resampled Line mean differences'); %title
xlabel('Difference of resampled means'); %x axis label
ylabel('#'); %y axis label
set(gcf, 'Position', [100 100 800 800]); %figure postition
set(gca, 'FontSize', fS); %font size on figure
set(gca, 'TickDir', 'out') %put ticks outside the lines
set(gca, 'XLim', [-.18 .04])
box off %turn the outer axes off

h= line([lEmpiricalMeanDiff lEmpiricalMeanDiff], [min(ylim) max(ylim)]); %make a vertical line in the plot from the top to the bottom
set(h,'Color', 'r'); set(h, 'LineWidth',lW); %set the line color and width

% ADDING THE CONFIDENCE INTERVAL TO THE FIGURE
sortedResampledMeans = sort(L_meanDifferences); %sort the mean differences in assending order
CIwidth = 99; % 99% CI
lowerBound = (100-CIwidth)/2; %set the lower bound
upperBound = 100-lowerBound; %set the upper bound

%find corresponding indices
lowerBoundIndex = round(lowerBound*length(sortedResampledMeans)/100); %find which row corresponds to the lower bound
upperBoundIndex = round(upperBound*length(sortedResampledMeans)/100); %find which row corresponds to the upper bound

% reach into the sorted mean differences to find the value of of the upper
% and lower CI
lowerCI = sortedResampledMeans(lowerBoundIndex); %set a value for the lower CI
upperCI = sortedResampledMeans(upperBoundIndex); %set a value for the upper CI

% Add the CI to the figure 
gca
lCI = line([lowerCI lowerCI], [min(ylim) max(ylim)]); set(lCI, 'color','g'); set(lCI, 'LineWidth', 3);
uCI = line([upperCI upperCI], [min(ylim) max(ylim)]); set(uCI, 'color','g'); set(uCI, 'LineWidth', 3);

%Adding a legend
legend('Null Dist','Test Statistic', 'Lower 99% CI', 'Higher 99% CI'); %add a legend
legend('Location', 'SouthEastOutside'); %place legend outside the plot to the lower right

shg

% ANALYSIS: Because the emperical mean fit squarely between the two
% confidence intervals we can assume our mean differnece is stable and not
% likely reached by chance.

pause %wait for button press to continue 



%% 6 PCA ON ALL CONDITIONS

% **NOTE: Given this data a PCA does not make the most sense since there
% are only 2x predictors that go into the output. However I wanted to show
% that I'm able to do it in Matlab, so let's suspend 'science' for a minute
% and just go along with it***


%Look at the spread of the raw data for all participants
figure %open a figure
set(gcf, 'Position', [100 100 1000 1200]); %figure position on the screen
set(gca, 'Color', 'w'); set(gcf, 'Color', 'w'); %set figure and background to white
subplot(3,1,1); %add a subplot
imagesc(AllTrialsMat); colormap(jet); colorbar;shg; %image the raw data
title('Raw Data Visualization'); %title
xlabel('Correct Response / Stim Category / # High Reliabiltiy') %label x axis
set(gca, 'FontSize', fS-6); %set font size
 
%See if any correlations exist between the outcomes and predictors.
subplot(3,1,2); %open a second subplot
imagesc(corrcoef(AllTrialsMat)); colormap(jet); colorbar;shg; %display the correlations between the predictors and output
title('Correlations Between Output and Predictors'); %title
set(gca, 'FontSize', fS-6); %font size
%There is a slight correlation between the number of high reliability
%ellipses and the the correct response.

%Run a PCA on the data
[loadings, origDataInNewDimensions, eigVal] = pca(AllTrialsMat); %run a PCA
varExplained = eigVal./sum(eigVal).*100; %pull out the variance explained 

eigenValMag = eig(corrcoef(AllTrialsMat)); %redoing just hte extraction of the eigen values from the corelation matrix
subplot(3,1,3); %add a 3rd subplot
bar(1:length(eigenValMag),sortrows(eigenValMag,-1)); %create a bar graph of the eigen values
title ('Scree Plot'); %title
line([min(xlim) max(xlim)], [1 1], 'linestyle', '--'); %kaiser criteron line
xlabel('Factors, in order of decreasing Eigenvalue'); %x axis label
ylabel('Eigenvalue'); %y axis label
set(gca, 'FontSize', fS-6); %font size

%Factor 1 is above the kaiser criteron line 

pause %wait for button press to continue 



%% 7 PERMUTATION TESTING ON BOTH CONDITIONS USING HIGH/LOW RELIABILITY

% Comparing the LOWEST reliability level and the HIGHEST relaibiltiy level
% for BOTH the ellipse and line conditions. 



%% 7a ELLIPSE CONDITION PERMUTATION TEST

empiricalE = mean(EZeroRelCorrect) - mean(EFourRelCorrect) %create a test statistic- in this case using the mean difference
% empericalE = -.1

HighLowRel = [EZeroRelCorrect; EFourRelCorrect]; %stack the ratings vertically
n1 = length(EZeroRelCorrect); %1 sample size
n2 = length(HighLowRel); %Combined sample size

N = 1e5; %Do this 100k times
nullDistributionE = nan(N,1); %preallocate to save time

for ii = 1:N
    % Shuffling randomly (drawing without replacment, or randomly permute)
    shuffleIndices = randperm(n2); %gets us a list of indices from 1 to length of HighLowRel
    % reach into the matrix of combined ratings twice
    randomizedScores1 = HighLowRel(shuffleIndices(1:n1)); %first half 
    randomizedScores2 = HighLowRel(shuffleIndices(n1+1:n2)); %second half
    
    nullDistributionE(ii,1) = mean(randomizedScores1)-mean(randomizedScores2); %save the randomized equilivent of the test statistic 
    
end

%PLOTTING THE NULL DISTRIBUTION 
figure %open a figure
set(gca, 'Color', 'w'); set(gcf, 'Color', 'w'); %set the figure and background to white
subplot(2,1,1); %open a subplot
histogram(nullDistributionE,100); %Make a histograph of the null distrion with bins of 100 
title('Null distribution of TestStat E (ellipse condition)'); %title
xlabel('Resampled E, assuming change'); %x axis label
ylabel('Number'); %y axis label
set(gcf, 'Position', [100 100 1000 800]); %figure position
set(gca, 'FontSize', fS); %figure font size
set(gca, 'XLim', [-0.15 0.15]) %range of x axis
set(gca, 'TickDir', 'out') %set ticks outside the box
box off %turn the outer edges off
hold on %wait for additional lines

% CREATING THE CONFIDENCE INTERVALS
sortNullDist = sortrows(nullDistributionE); %put the matrix of null scores in order
highCI = sortNullDist(length(sortNullDist)*.005); %create a variable that marks the row of the lowest .5% of scores
lowCI = sortNullDist(length(sortNullDist)*.995); %create a variable that marks the row of the highest .5% of scores

%Add CI lines to the histogram
line([lowCI, lowCI], ylim, 'LineWidth', 2, 'Color', 'r'); %add a stright line at the lower bound of the 99% CI in red
line([highCI, highCI], ylim, 'LineWidth', 2, 'Color', 'r'); %add a stright line at the higher bound of the 99% CI in red

%Adding the test statistic
EStatLine = line([empiricalE empiricalE], [min(ylim) max(ylim)]);
set(EStatLine, 'Color', 'g', 'Linewidth',2); 

%Adding the legend
legend('Null Dist', 'Lower 99% CI', 'Higher 99% CI','Test Statistic'); %add a legend
legend('Location', 'SouthEastOutside'); %place legend outside the plot to the lower right
shg

EexactP = sum(nullDistributionE<empiricalE)/length(nullDistributionE) %calculate the exact p-value
%exact p-value = 0



%% 7b LINE CONDITION PERMUTATION TEST

empiricalL = mean(LZeroRelCorrect) - mean(LFourRelCorrect)
%empericalL = -.04

HighLowRel = [LZeroRelCorrect; LFourRelCorrect]; %stack the ratings vertically
n1 = length(LZeroRelCorrect); %1 sample size
n2 = length(HighLowRel); %Combined sample size

N = 1e5; %100k times
nullDistributionL = nan(N,1); %preallocate to save time

for ii = 1:N
    % Shuffling randomly (drawing without replacment, or randomly permute)
    shuffleIndices = randperm(n2); %gets us a list of indices from 1 to length of HighLowRel
    % reach into the matrix of combined ratings twice
    randomizedScores1 = HighLowRel(shuffleIndices(1:n1)); %first half
    randomizedScores2 = HighLowRel(shuffleIndices(n1+1:n2)); %second half
    
    nullDistributionL(ii,1) = mean(randomizedScores1)-mean(randomizedScores2); %save the randomized null version of the test statistic
    
end

%PLOTTING THE NULL DISTRIBUTION 
subplot(2,1,2); %open a 2nd subplot
histogram(nullDistributionL,100); %draw a histogram of the null dist. with bins of 100
title('Null distribution of TestStat L (line condition)'); %title
xlabel('Resampled L, assuming change'); %x axis label
ylabel('Number'); %y axis label
set(gcf, 'Position', [100 100 1000 800]); %figure position
set(gca, 'FontSize', fS); %figure font size
set(gca, 'XLim', [-0.15 0.15]) %x axis range
set(gca, 'TickDir', 'out') %place ticks outside the box
box off %turn the outer axes off
hold on; %wait to add another element to the figure

sortNullDist = sortrows(nullDistributionL); %put the matrix of null scores in order
highCI = sortNullDist(length(sortNullDist)*.005); %create a variable that marks the row of the lowest 5% of scores
lowCI = sortNullDist(length(sortNullDist)*.995); %create a variable that marks the row of the highest 5% of scores

%Add CI lines to the histogram
line([lowCI, lowCI], ylim, 'LineWidth', 2, 'Color', 'r'); %add a stright line at the lower bound of the 99% CI in red
line([highCI, highCI], ylim, 'LineWidth', 2, 'Color', 'r'); %add a stright line at the higher bound of the 99% CI in red

%Adding the test statistic
LStatLine = line([empiricalL empiricalL], [min(ylim) max(ylim)]);
set(LStatLine, 'Color', 'g', 'Linewidth',2);

%Adding the legend
legend('Null Dist', 'Lower 99% CI', 'Higher 99% CI','Test Statistic'); %add legend
legend('Location', 'SouthEastOutside'); %place legend outside the plot to the lower right
shg

LexactP = sum(nullDistributionL<empiricalL)/length(nullDistributionL) %Calcualte the exact p value
%LexactP = .01999

pause %wait for button press
%% 7c LINE VS ELLIPSE (4 rel) PERMERATION TEST 

%Looking at the mean difference between the ellipse condition and the line
%condition when there are ALL narrow ellipses

empiricalEL = mean(EFourRelCorrect) - mean(LFourRelCorrect)
%empericalL = .07

HighLowRel = [EFourRelCorrect; LFourRelCorrect]; %stack the ratings vertically
n1 = length(EFourRelCorrect); %1 sample size
n2 = length(HighLowRel); %Combined sample size

N = 1e5; %100k times
nullDistributionEL = nan(N,1); %preallocate to save time

for ii = 1:N
    % Shuffling randomly (drawing without replacment, or randomly permute)
    shuffleIndices = randperm(n2); %gets us a list of indices from 1 to length of HighLowRel
    % reach into the matrix of combined ratings twice
    randomizedScores1 = HighLowRel(shuffleIndices(1:n1)); %first half
    randomizedScores2 = HighLowRel(shuffleIndices(n1+1:n2)); %second half
    
    nullDistributionEL(ii,1) = mean(randomizedScores1)-mean(randomizedScores2); %save the randomized null version of the test statistic
    
end

%PLOTTING THE NULL DISTRIBUTION 
figure; %open a new figure
histogram(nullDistributionEL,100); %draw a histogram of the null dist. with bins of 100
title('Null distribution of TestStat EL (line vs ellipse condition)'); %title
xlabel('Resampled EL, assuming change'); %x axis label
ylabel('Number'); %y axis label
set(gcf, 'Position', [100 100 1000 800]); %figure position
set(gca, 'FontSize', fS); %figure font size
set(gca, 'XLim', [-0.15 0.15]) %x axis range
set(gca, 'TickDir', 'out') %place ticks outside the box
set(gca, 'Color', 'w'); set(gcf, 'Color', 'w'); %set the figure and background to white
box off %turn the outer axes off
hold on; %wait to add another element to the figure

sortNullDist = sortrows(nullDistributionEL); %put the matrix of null scores in order
highCI = sortNullDist(length(sortNullDist)*.005); %create a variable that marks the row of the lowest 5% of scores
lowCI = sortNullDist(length(sortNullDist)*.995); %create a variable that marks the row of the highest 5% of scores

%Add CI lines to the histogram
line([lowCI, lowCI], ylim, 'LineWidth', 2, 'Color', 'r'); %add a stright line at the lower bound of the 99% CI in red
line([highCI, highCI], ylim, 'LineWidth', 2, 'Color', 'r'); %add a stright line at the higher bound of the 99% CI in red

%Adding the test statistic
LStatLine = line([empiricalEL empiricalEL], [min(ylim) max(ylim)]);
set(LStatLine, 'Color', 'g', 'Linewidth',2);

%Adding the legend
legend('Null Dist', 'Lower 99% CI', 'Higher 99% CI','Test Statistic'); %add legend
legend('Location', 'SouthEastOutside'); %place legend outside the plot to the lower right

shg

LexactP = sum(nullDistributionEL>empiricalEL)/length(nullDistributionEL) %Calcualte the exact p value
%LexactP = 0



%% 8 SVM MODEL & RANDOM FOREST MODEL ON DATA FOR ALL PARTICIPANTS

% Doing an 80/20 test/train split. pulls 80% of trials at random into the
% 'train' group and then the other 20% of trials into the 'test' group. 

ll = length(AllTrialsMat); %find length of data
predictors = AllTrialsMat(:,2:3); %isolates only the predictor variables
outcomes = AllTrialsMat(:,1); %isolates the outcome variable

randOrd = randperm(ll); %get a random order of numbers the length of the data matrix
tt = round(ll*.8); %find what trial number is the 80% cut off
ii = randOrd(1:tt); %80% of trials
pp = randOrd(tt+1:ll); %20% of trials

for rr = 1:tt %for trials 1 through the 80% cut off
    numSelect = ii(rr); %find the rr'th row of the randomized set of numbers (equal in length to rr)
    outTrain(rr,:) = outcomes(numSelect,:); %make a variable for the outputs to train on (80% of trials) by selecting them at random from the total trials
    predTrain(rr,:) = predictors(numSelect,:); %Make a matching variable of the predictors to use the same 80% of trials as the outputs
end
for rr = 1:length(pp) %for trials 1 through last 20%
    numSelect = pp(rr); %find the rr'th row of the randomized set of numbers (equal to rr in length)
    outTest(rr,:) = outcomes(numSelect,:); %create a variable of outputs to test on containing the other 20% of trials from the original trial set
    predTest(rr,:) = predictors(numSelect,:); %create a predictor variable to test on using the same 20% of trials 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SVM MODEL 
svmModel = fitcsvm(predTrain, outTrain); %creates the model using the amount of trials declared above
[decision score] = predict(svmModel,predTest); %tests the model using the unseen predictors
comp = [decision outTest];  %compares the model with the unseen outcomes vs decision
SVMmodelAccuracy = (sum(comp(:,1) == comp(:,2))./length(comp)*100) %percentage of model accuracy


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RANDOM FOREST MODEL
numTree =  5; %set the number of trees input by participant

%run the treebagger model
treeModel = TreeBagger(numTree,predTrain,outTrain); %using the number of trials specified above
[treeDecisions score] = predict(treeModel,predTest); %use the model to make a prediction using unseen predictors

%The output of the treebagger are labels, we need to convert to numbers
numericalOutcomes = str2num(char(string(treeDecisions))); %get the numerical output from the prediction
empericalVsPrediction = [outTest numericalOutcomes]; %compare the unseen outcomes to the predicted values
forestModelAccuracy = (sum(empericalVsPrediction(:,1) == empericalVsPrediction(:,2))./length(empericalVsPrediction)*100) %calculate percent model accuracy

%Both model accuracies usually end up in the 66% range. As guessing
%randomly would get an accuracy of 50% in this case as there are two
%options, this is only a 16% improvement on chance. I feel that with a
%model that takes into account other weights of the particular reliability levels
%a more accurate outcome would be achieved. 

%% 9 CONCLUSIONS BASED ON THE DATA ANALYSIS

%The ellipse condition performance is strongly effected by the increase of high
%reliability items, shown by the significant difference between the 0 level
%and the 4 level for that condition. 

%On the otherhand, the line conditon performance does not seem to be
%effected by the increase of high reliability items, as the difference
%between performance at the 0 level and the 4 level are not significantly
%different. 

%The two different conditions themselves result in significantly different
%performance. 