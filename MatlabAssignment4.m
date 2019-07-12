%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Date Created: June 27th 2019
% Date Last Edit: June 237th 2019
% Author: Marissa Evans, mhe229@nyu.edu

% Assumptions: The files 'LOY.mp4', 'LSRP.xlsx' & 'soundSignals.mat' are in
% the same directory as the code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script excecutes signal analysis, PCA and permutation tests. 

% What this code does(by section):

% 0) Initialization for script

%%%%%%%%%%%%%%%%%%%% SECTION 1 - soundSignals DATA %%%%%%%%%%%%%%%%%%%%%%%%
% 1) Load/read the file 'soundSignals.mat' - Initalize for section 1
% 2) Plays each of the 3x sounds in soundSignals.mat separated by pauses,
% then visualizes the signal in the time and frequency domain to find the
% components. 

%%%%%%%%%%%%%%%%%%%%%%%% SECTION 2 - LOY DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load file'LOY.mp4' and initalize for section 2
% 3) Uses the 'audioread' function to get he signal and sampling rate from
% the LOY.mp4 file. Description of sound followed by visualization to match
% article (http://bit.ly/LOY2019S). Signal then replayed with 20% slower
% and faster sampling rates. 

%%%%%%%%%%%%%%%%%%%%% SECTION 3 - LSRP DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load file 'LSRP.xlsx' and initalize for section 3
% 4) Using LSRP.xlsx recode questions 6,14,19,22,24,25,&26 as they are inversly
% coded. Plot of raw data is made, then the missing data is handled. 
% Then a plot of the processed raw data is excuted. 
% After that the total score is calculated for each participant (26 item
% personality test) and a histogram is made of emperical scores with emperical
% mean denoted. Note: 0= female, 1 = male
% 5) Runs a PCA and creates scree plot to assess which factors should be
% retained and are meaningful by either elbow or kaiser criterion.
% Calculates how much total variance is explained by these factors. 
% 6) Resample the SLRP data and perform permutation tests to determine if
% there are potential gender differences in psychopathy. 
% 7) Save the workspace output.

% OUTPUTS:
% Audio files play so be sure speaker sound is audible but not too loud
% Seven figures open with pauses in between. PRESS ANY KEY TO PROCEED TO
% THE NEXT FIGURE

%% 0 INITIALIZATION
clear all %clear workspace
close all %close all open figures
clc %clear the command window

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% SECTION 1 - soundSignals DATA %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD AND VARIABLE INITALIZATION

load('soundSignals.mat'); %Load the soundSignals mat file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2a PLAYING SOUNDS

figure
set(gcf, 'Position', [100 100 600 800]); %position figure
sound(signal1,fs) %play sound 1
plot(signal1) %visualize sound 1
title('Signal 1')
xlabel('Time')
ylabel('Amplitude')
pause %wait for key press
sound(signal2,fs) %play sound 2
plot(signal2) %visiualize sound 2
title('Signal 2')
xlabel('Time')
ylabel('Amplitude')
pause %wait for key press
sound(signal3,fs) %play sound 3
plot(signal3) %visualize sound 3
title('Signal 3')
xlabel('Time')
ylabel('Amplitude')
ylim([-1 1])
pause


%% 2b VISUALIZATIONS

windowShape = hanning(512); %several windows smear over time but relflect the frequency clearly 
overL = (length(windowShape)/2); %maximal overlap
frequencyRange = 0:62:fs/2; %entire time course is not absolutly necessary however incluing for continuity sake. 
samplingFrequency = 1e4; %sampling 10k/s

figure %Open a figure
set(gcf, 'Position', [100 100 600 800]); %position figure

%Spectrogram of sound 1
subplot(3,1,1)
spectrogram(signal1,windowShape, overL, frequencyRange, samplingFrequency, 'yaxis')
title('Signal 1')

%Spectrogram of sound 2
subplot(3,1,2) 
spectrogram(signal2,windowShape, overL, frequencyRange, samplingFrequency, 'yaxis')
title('Signal 2')

%Spectrogram of sound 3
subplot(3,1,3)
spectrogram(signal3,windowShape, overL, frequencyRange, samplingFrequency, 'yaxis')
title('Signal 3')

% All three signals are made up of a 20kHz tone and a 10kHz tone. The 10kHz
% tone remains at a consitant power throughout the signal while the power
% of the 20kHz signal varies in power in different ways across the 3x
% signals. 

pause %wait for key press to show next figure
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% SECTION 2 - LOY DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD AND INITALIZE FOR SECTION 2

[loyDATA, loyFS] = audioread('LOY.mp4'); %read the audio file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3a LOY AUDIO SPECTROGRAM

sound(loyDATA); %sounds like a strange lion roar 

loyWindow = bartlett(1024); %tried several options, this looks the most like the reference image
loyOverL = round(length(loyWindow)/2); %maximal windows
loyFreqRange = 0:16:loyFS/7; %Shorter frequency range 

%Spectrogram of LOY.mp4
figure %open figure
set(gcf, 'Position', [100 100 600 800]); %position figure
spectrogram(loyDATA, loyWindow, loyOverL, loyFreqRange,loyFS, 'yaxis')
colormap(jet)
title('Spectrogram of LOY data')

pause %wait for key press to show next figure
%% 3b CHANGING THE SAMPLING RATE

sound(loyDATA, loyFS*.8) %20% slower sampling rate- sounds like Yanni
pause %wait for key press
sound(loyDATA, loyFS*1.2) %20% higher sampling rate- sounds like Laurel 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SECTION 3 - LSRP DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD AND INITALIZE FOR SECTION 3

lsrpDATA = xlsread('LSRP.xlsx'); %read the Excel file

inverseCode = [6 14 19 22 24 25 26]; %questions with invered coding
fontS = 20; %font size for figures
validQuestions = 1:26; %the columns which reprsent valid quesitons
gender = 27; %denotes the column in which 'gender' resides in the LSRP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4a LSRP DATA ANALYSIS - REVERSE CODING

%view unprocessed raw data
figure %open figure
set(gcf, 'Position', [100 100 600 800]); %position figure
subplot(2,1,1) %add a subplot
imagesc(lsrpDATA); colormap(jet); colorbar; %image data and set color map
title('Unprocessed Raw LSRP Data')
ylabel('Participants')
xlabel('Questions')

%Recode inverse coded questions to correct order
recodeLSRP = lsrpDATA; %copy matrix

%Loop that looks through each column specified and makes applicable changes
for ii = inverseCode
    oneToFive = find(recodeLSRP(:,ii)==1); %find all 1's in column ii
    recodeLSRP(oneToFive, ii) = 5; %replace with 5's
    twoToFour = find(recodeLSRP(:,ii)==2); %find all 2's in column ii
    recodeLSRP(twoToFour, ii) = 4; %replace with 4's
    fourToTwo = find(recodeLSRP(:,ii)==4); %fine all 4's in column ii
    recodeLSRP(fourToTwo, ii) = 2; %replace with 2's
    fiveToOne = find(recodeLSRP(:,ii)==5); %find all 5's in column ii
    recodeLSRP(fiveToOne, ii) = 1; %replace with 1's 
end

%% 4b LSRP - REMOVING NaN's

%removing NaN's participant wise as only 42 of 606 are lost that way. Plus
%it keeps the data lined up accross quesitons without removing whole
%quesitons. 
findNans = sum(isnan(recodeLSRP),2); %counts missing values per row (the number is the DIMENSION!!) 1 = columns, 2=rows
rowNans = find(findNans==0); %only keep rows of participants that answered all quesitons
cleanRecode = recodeLSRP(rowNans,:); %rewrite the scores including only those who answered all questions 

%visualize raw data after inverting and cleaning
subplot(2,1,2)
imagesc(cleanRecode); colormap(jet); colorbar;
title('Cleaned Raw LSRP Data')
ylabel('Participants')
xlabel('Questions')

pause %wait for key press to show next figure
%% 4c LSRP EMPERICAL COMPARISONS

empericalScores = sum(cleanRecode(:,validQuestions)')'; %the total score for each participant- verticly 

empFig = figure; %open a figure
set(empFig, 'Position', [100 100 600 800]); %set figure size

%Put emperical score data on figure
handleHist = histogram(empericalScores,10); %create a histrogram with 10 bins
handleHist.FaceAlpha = .5; %make histogram bars slightly translucent
handleHist.FaceColor = 'b'; %make face color blue
handleHist.EdgeColor = 'k'; %make edge color black

%label axis, title and legend
xlabel('Total LSRP Score'); %x axis label
ylabel('Number of Participants'); %y axis label
title('Total LSRP Score Across Participants'); %plot title

%Physical atributes of the figure
set(gca,'TickDir', 'out'); %put the ticks outside the plot
set(gca, 'Color', 'w'); set(empFig, 'Color', 'w'); %set plot and background color to white
set(gca, 'FontSize', fontS); %set fontsize to indexed size above
box off %turn off the box around the plot
shg %show figure

pause %wait for key press to show next figure
%% 5 RUNNING + PLOTTING PCA

% PLOTTING THE CORRELATIONS 
figure %Open a figure
set(gcf, 'Position', [100 100 600 800]) %resize the figure

%Plot the correlations between the quetions 
subplot(2,1,1)
imagesc(corrcoef(cleanRecode(:,validQuestions))); colormap(jet); colorbar;shg
title('LSRP Question Correlations')

%RUN THE PCA
[eigVec, rotVal, eigVal] = pca(cleanRecode(:,validQuestions));

%calculate the variance explained
varExplained = eigVal./sum(eigVal).*100;

% Extracting eigen values
eigenValMag = eig(corrcoef(cleanRecode(:,validQuestions)));

%Plotting a bar graph of the eigen values with a Kaiser criteron line
subplot(2,1,2)
bar(1:length(eigenValMag),sortrows(eigenValMag,-1))
title ('LSRP Scree Plot')
line([min(xlim) max(xlim)], [1 1], 'linestyle', '--') %kaiser criteron line
ylabel ('Eigen Value')
xlabel('Factors')

% Using the kaiser criteron line the first 6x factors should be meaninfully
% maintained. These factors account for approximatly 52% of the total variance
% when combined. 

% Using the elbow criteron either 1 or 3 factors could be maintained. These
% woudld account for either approximatly 22% or 37% of the total variance when
% combined. 

pause %wait for key press to show next figure
%% 6a RESAMPLING + PERUTATION BY GENDER

%Create groups for compairson
women = find(cleanRecode(:,gender)==0); %find all women
womenAll = cleanRecode(women,validQuestions); %put scores of women in matrix
men = find(cleanRecode(:,gender)==1); %find all men
menAll = cleanRecode(men,validQuestions); %put scores of men in matrix

%Create test statistic - comparing the difference in average ratings on each
%question over all males and females. Test Stat = -4.65
womenScores = sum(womenAll)/length(womenAll); %Collapsing across participants to sum question scores (then dividing by participants
menScores = sum(menAll)/length(menAll); %Collasping across participants to sum question scores (dividied by participants)
testStatM = sum(womenScores-menScores); %using a difference between the group scores as a test statistic

%Combine groups
bothGenders = [womenScores menScores]; %horzcat to stack groups
n1 = length(womenScores); %1st sample size
n2 = length(bothGenders); %Combined sample size

numReps = 1e5; %100k reps
nullDistM = nan(numReps,1); %preallocate the storage matrix

for ii = 1:numReps
    shuffleIndices = randperm(n2); %creates a list of indices from 1 to length of bothGenders (in random order)

    randomizedRatings1 = bothGenders(shuffleIndices(1:n1)); %first half
    randomizedRatings2 = bothGenders(shuffleIndices(n1+1:n2)); %second half
    
    nullDistM(ii,1) = sum(randomizedRatings1-randomizedRatings2); %null version of test stat. 
    
end

%% 6b PLOTTING THE NULL DIST

sortNullDist = sortrows(nullDistM); %put the matrix of null scores in order 
highCI = sortNullDist(length(sortNullDist)*.025); %create a variable that marks the row of the lowest 2.5% of scores 
lowCI = sortNullDist(length(sortNullDist)*.975); %create a variable that marks the row of the highest 2.5% of scores

figure
set(gcf, 'Position', [100 100 1000 800]); %position figure

%Adding the null distribution
histogram(nullDistM,100)
title('Null distribution of Test Statistic M using LSRP')
xlabel('Resampled M, assuming change')
ylabel('#')
set(gca, 'FontSize', 26)
hold on;

%Add CI lines to the histogram
line([lowCI, lowCI], ylim, 'LineWidth', 2, 'Color', 'r'); %add a stright line at the lower bound of the 95% CI in red
line([highCI, highCI], ylim, 'LineWidth', 2, 'Color', 'r'); %add a stright line at the higher bound of the 95% CI in red

%Adding the test statistic 
MStatLine = line([testStatM testStatM], [min(ylim) max(ylim)]); 
set(MStatLine, 'Color', 'g', 'Linewidth',2)
legend('Null Dist', 'Lower 95% CI', 'Higher 95% CI','Test Statistic')
legend('Location', 'SouthEastOutside') %place legend outside the plot to the lower right
shg

%The test statistic falls within the null distribution (confirmed by the
%confidence interval). Therefor the two groups sampled are not
%significantly different 

% Calculate the exact p-value
exactP = sum(nullDistM>testStatM)/length(nullDistM);

% Exact p-value equals .9135 - the LSRP scores that denote psychopathy do
% not occur significantly more frequently in men or women on a question by
% question basis. 