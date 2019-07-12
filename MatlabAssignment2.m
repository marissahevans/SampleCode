%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% DATA RECORDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Date Created: June 14th 2019
% Date Last Edit: June 19th 2019
% Author: Marissa Evans - mhe229@nyu.edu

% Assumptions: The file 'thedress.jpg' is in the same directory as the code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script writes a psychophysics experiment in which participants
% determin if an image of 'the dress' looks bluer in either the original
% image or an image with altered brightness when presented side by side.

% Total trials: 110 (10 per condition)

% What this code does(by section):
% 0) Initialization
% 1) Load the file 'thedress.jpg'
% 2) Make 11 copies of the dress with a varying +/- scale of brightness
% 3) Present the original image side by side with one of the copies in a
% stimulus display screen with a red central fixation cross with the
% original appearing on each side of the display an equal number of times.
% 4) Organization of presentation: timing, response keys, response window,
% RT recording
% 5) Calcuate the mean response of 'bluer' as a function of condition
% 6) Plot a psychometric curve of the responses
% 7) Save the workspace output.


%% 0 INITIALIZATION

clear all %clear workspace
close all %close all open figures
clc %clear command window

brightSet = [-130 -110 -90 -70 -50 50 70 90 110 130 150]; % create a scale of brightness adjustments above and below the original image
numSet = 11; %number of iterations that loops will go though (ie how many conditions there are)
fS = 20; %font size for stimuli
numTrials = 110; %total number of trials
numReps = 1:10; %number of repitions each condition has
stimSize = get(0,'ScreenSize'); %finds the screen size of the current display
stimGrey = 1/255*[238,238,238]; %creates a neutral grey for stimulus background

%fixation cross
crossSize = 10; %total height of fixation cross
crossWidth = 3; %width of cross legs

% Make a 11x10 matrix with a row for each condition (10 stimuli per
% condition, 11 conditions) filled with numbers 1:110
condiAll = reshape(1:110, 10, 11).';

% For psychometric plot
startingNum = 0; %lowest probability
endingNum = 1; %highest proability
hair = .05; %breathing room spacing on plot
fSP = 14; %font size on plot
lW = 3; %line width on plot

%% 1 LOAD FILES

dress = imread('thedress.jpg'); %load the file and assign it to a matrix (will make it a uint8 file) *assuming file is in directory

%% 2 MAKE 11 COPIES OF STIMULI WITH VARYING BRIGHTNESS

dressEdit = cell(11,1); %create a cell array for all dress matricies

for ii = 1:numSet %go through loop 11x times
    for bb = brightSet(ii) %call each iteration of brightness one time
        dressEdit{ii} = dress + bb; %place each dress copy in a different row of cell array, updated version is the original with the brightness change ADDED. 
    end
end

%% 3 MAKE PAIRS OF STIMULI AROUND FIXATION POINT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3a) Pair the images with random left/right placement using an odd/even
% number metric based on a randomized intiger between 1 and 100

% NOTE: because each arrangement is randomized by a coin flip (as opposed to a
% shuffled random) the selections are independent from eachother. 
% Thus there is not always going to be an equal amount of right and left
% arrangements. It is possible to get a set of all right or all left
% (although rather unlikely). 

pairedDresses = cell(11,2); %creates a cell array to hold paired dress images
rng('shuffle') %sets random number seed using clock vector at high precision

for ii = 1:numSet %run through loop 11 times
    randomNumber = randi(100); %finds a random number between 1 and 100
    if rem(randomNumber,2) == 0 %if random number is EVEN
        pairedDresses{ii} = [dress dressEdit{ii}]; %put the original image on the LEFT
        pairedDresses{ii,2} = 1; %output of 1 shows that edited image is on the RIGHT
    elseif rem(randomNumber, 2) ~= 0 %if the random number is ODD
        pairedDresses{ii} = [dressEdit{ii} dress]; %put the original image on the RIGHT
        pairedDresses{ii,2} = 2; % output of 2 shows that edited image is on the LEFT
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3b) add the fixation point to all paired images (size of cross initalized
% in '0'). Cross is not projected onto image but rather pixels in it's area
% are changed to red. 

for ii = 1:numSet %go through loop 11 times
    dressMatSize = size(pairedDresses{ii}); %get size of matrix of paired dresses
    nx = dressMatSize(1,2); %calcuate the x axis size
    ny = dressMatSize(1,1); %calculate the y axis size
    centerPoint = round([nx ny]/2); %find the center point for the image using x and y
    
    %Fitting the size of the fixation cross on the rows and changing those
    %pixels to red for all eleven pairs
    pairedDresses{ii}(centerPoint(2)-crossSize:centerPoint(2)+crossSize,centerPoint(1)-crossWidth:centerPoint(1)+crossWidth,1) = 255;
    pairedDresses{ii}(centerPoint(2)-crossSize:centerPoint(2)+crossSize, centerPoint(1)-crossWidth:centerPoint(1)+crossWidth,2:3) = 0;
    
    %Fitting the fixation cross on the columns and changing the included
    %pixels to red for all eleven pairs
    pairedDresses{ii}(centerPoint(2)-crossWidth:centerPoint(2)+crossWidth, centerPoint(1)-crossSize:centerPoint(1)+crossSize,1) = 255;
    pairedDresses{ii}(centerPoint(2)-crossWidth:centerPoint(2)+crossWidth, centerPoint(1)-crossSize:centerPoint(1)+crossSize,2:3) = 0;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3c) make a blank screen with fixation for pause during response

blankDisp = pairedDresses{1}; %Copy a matrix of the paired dresses at the same size
blankDisp(:,:,:) = 238; %change all pixels to grey

%Fitting the size of the fixation cross on the rows and changing those
%pixels to red
blankDisp(centerPoint(2)-crossSize:centerPoint(2)+crossSize,centerPoint(1)-crossWidth:centerPoint(1)+crossWidth,1) = 255;
blankDisp(centerPoint(2)-crossSize:centerPoint(2)+crossSize, centerPoint(1)-crossWidth:centerPoint(1)+crossWidth,2:3) = 0;

%Fitting the fixation cross on the columns and changing the included
%pixels to red
blankDisp(centerPoint(2)-crossWidth:centerPoint(2)+crossWidth, centerPoint(1)-crossSize:centerPoint(1)+crossSize,1) = 255;
blankDisp(centerPoint(2)-crossWidth:centerPoint(2)+crossWidth, centerPoint(1)-crossSize:centerPoint(1)+crossSize,2:3) = 0;

%% 4 TRIAL PRESENTATIONS, TIMING AND RESPONSE

% 4a) Loops that make 10x copies of each image pair and puts them in a cell
% array with 110 spaces with darker pairs in the lower rows and brighter
% pairs in the higher rows. 

dressStim110 = cell(110,1); %initalize cell array for all 110 trials
dressPlacement = zeros(numTrials,1); %initalize a matrix to store if the edited image was presented on the right or left for all 110 trials

for ii = 1:numSet %run through loop 11 times
    for pp = condiAll(ii,:) %go though all iterations of each condition one at a time
        dressStim110{pp} = pairedDresses{ii}; %place the paired dress in 10x cells in the array
        dressPlacement(pp) = pairedDresses{ii,2}; %store if the paired dress had the edited version on the right or left
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4b) BUILD THE STIMULUS PRESENTATION AND RECORD RESPONSE

% Trial order randomization is achieved by creating a vector of numbers
% 1-110 in a random order, then using that as the call order to pull images
% out of the stimulus matrix. No numbers are repeated and all 110 images
% are presented. 

% Stimulus is displayed for 500 ms followed by a blank screen with fixation
% until participant responds. The reaction time timer starts once the
% participant is able to respond. 

%Output is stored depending on participant response in a separate matrix. 

validKeyPress = 0; %initalize the valid key press
rr = randperm(110)'; %create a random set of numbers from 1 to 110 without repetition
reactionTime = ones(numTrials,5)*999; %initalize a matrix to store all output using 999 since 1 and 0 are output options
dressFig = figure; %initalize the stimulus figure

%physical properties of stimulus figure
set(dressFig, 'Color', stimGrey) %grey background on figure
set(dressFig, 'Position', stimSize) %make the stimulus display at full screen
axis off %turn the axis off

%Welcome page with instructions for task
startHandle = text(0.5,0.8,'For each trial images will appear to the right and left of central fixation');
startHandle.FontSize = fS + 10; %font size is 10pts larger than experimental font size
startHandle.HorizontalAlignment = 'center'; %text is centered on figure

startHandle = text(0.5,0.6,'Please select which of the two images appears more blue.');
startHandle.FontSize = fS + 10; %font size is 10pts larger than experimental font size
startHandle.HorizontalAlignment = 'center'; %text is centered on figure

startHandle = text(0.5, 0.4, 'Press "j" for the left image, and "k" for the right image.');
startHandle.FontSize = fS + 10; %font size is 10pts larger than experimental font size
startHandle.HorizontalAlignment = 'center'; %text is centered on figure

startHandle = text(0.5, 0.2, 'Press any key to begin.');
startHandle.FontSize = fS + 10; %font size is 10pts larger than experimental font size
startHandle.HorizontalAlignment = 'center'; %text is centered on figure
pause %wait for key press
clf %clear current figure

%Task with image displays and response recording

for ii = 1:numTrials %go though all 110 trials
    clear nn %clear randomized trial number variable each time
    nn = rr(ii,1); %pick one of the randomized set of 110 numbers without any repeating.
    validKeyPress = 0; %update valid key press within the loop
    
    %Image display
    image(dressStim110{nn}) %open image of the paired dresses based on random order so presentation is mixed.
    set(gca, 'FontSize', fS) %set font size on the stimulus
    axis off
    
    while validKeyPress == 0 %to start loop must match initalization above
        
        %Stimulus timing
        pause(.5) %display image for 500 ms
        clf %clear the current figure
        
        %Response display (blank) display properties
        image(blankDisp) %display a blank image with the central fixation cross until response
        set(dressFig, 'Color', stimGrey) %grey background on figure
        set(gca, 'FontSize', fS) %set font size on the stimulus
        set(dressFig, 'Position', stimSize) %make stimulus display at full screen
        title('Select which image is more blue, "j" for LEFT and "k" for RIGHT') %title of the figure
        axis off %turn the edges of the figure off
        
        tic; %start timer for reaction time
        pause %wait for for the key press
        
        % IF the response key is pressed too quickly the command window
        % will pop up in front of the stimulus display, just click back on
        % the figure to continue the experiment. 
        
        % If the 'bluer' image is on the left
        if strcmp(dressFig.CurrentCharacter,'j') == 1 %checking if 'j' was pressed (image on left selected)
            reactionTime(ii,2) = toc; %save the reaction time
            reactionTime(ii,1) = ii; %save the trial number
            reactionTime(ii,3) = nn; %save which image was presented
            reactionTime(ii,4) = 2; %output of '2' means that participant selected image on the LEFT
            validKeyPress = 1; %update valid key pressed, goes back to the start of the loop.
            pause(.3) %ITI of 300 ms
        % If the 'bluer' image is on the right
        elseif strcmp(dressFig.CurrentCharacter,'k') ==1 %if j was not pressed, check if 'k' was pressed
            reactionTime(ii,2) = toc; %save the reaction time
            reactionTime(ii,1) = ii; %save the trial number
            reactionTime(ii,3) = nn; %save which image was presented
            reactionTime(ii,4) = 1; %output of '1' means that participant selected image on the RIGHT
            validKeyPress = 1; %update valid key pressed, goes back to the start of the loop.
            pause(.3) %ITI of 300 ms
        % If incorrect key is pressed
        else
            disp('Invalid Entry')
            %no update here for validKeyPressed as input was not correct.
        end
        clf %clear the figure
    end
end
close all %close the open stimulus figure

%% 5 SORTING TRIALS BY CONDITION AND LOGGING RESPONSE BY CONDITION

reactionTime = sortrows(reactionTime,3); %Sorts the output by the condition type (the trial numbers in condiAll are how the conditions are divided)

%Adding another column to reactionTime that records if the original or
%edited image was selected as bluer accounting for the randomized position
%changes.

reactionTime(:,5) = reactionTime(:,4) == dressPlacement(:,1); %1 means the edited version was selected, 0 means the original was selected

probRespBlue = ones(numSet,2)*99; %initalize a matrix for the probability that the edited version was selected

for ii = 1:numSet %go though loop 11 times
    probRespBlue(ii,1) = mean(reactionTime(condiAll(ii,:),5)); %calculate the mean probability of response for each condition over it's ten trials
end


%% 6 PLOT DATA IN PSYCHOMETRIC CURVE

% Plot the probability that the participant responded that the edited
% version was bluer than the original version 

probRespBlue(:,2) = brightSet; %add a column to the matrix for how bright each condition was set at

x=probRespBlue(:,2); %x-axis will be the brightness level
y=probRespBlue(:,1); %y-axis will be the probabily that the edited version was thought to be bluer

%Make Plot
PsychFig = figure; %Open a figure
set(PsychFig, 'Position', [100 100 800 800]) %set figure size
respPlot = plot(x,y); %plot the psychometric function of percieved blueness as a factor of brightness

%label axis, title and legend
xlabel('Brightness Condition') %x axis label
ylabel('Probabiltiy Respond Blue') %y axis label
title('Probability of Blue Response Based on Brightness Condition') %plot title
legend('Participant 1') %add legend
legend('Location', 'SouthEastOutside') %place legend outside the plot to the lower right

%Physical atributes of the plot
set(gca, 'YLim', [startingNum - hair endingNum + hair]) %add a hair of space on either end of the plot
set(gca,'TickDir', 'out') %put the ticks outside the plot
set(respPlot, 'LineWidth', lW) %set the line width on the plot
set(gca, 'Color', 'w'); set(PsychFig, 'Color', 'w') %set plot and background color to white
set(gca, 'FontSize', fSP) %set fontsize to indexed size above
box off %turn off the box around the plot

shg %show figure

%% 7 Save the Output

save('MatlabAssignment2') %saves .mat file to current directory.
