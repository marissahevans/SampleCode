%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% CAFFEINE REACTION TIMES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0 Script Information

% Created: June 6th 2019
% Author contact: Marissa Evans, mhe229@nyu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data file 'Assignment1_data.mat' contains the reaction times (loaded
% below). Should be inside the same directory as the script. 

% This data represents 3x groups of participants, 200mg caffeine, 400mg
% caffeine, and a control group (represented by the data_ variables). 

% The data measures their reaction times in each of the three groups. The stimulus 
% onset time for each trial is jittered (represented by the stimTime_ variables). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script calculates the actual reaction time after stimulus onset and
% removes any trails that included a false start. 
% All 'reaction times' x 'valid trials' are then plotted for each participant 
% on one figure (black=control, blue=200mg caffeine, red=400mg caffeine) 

% Workspace is then saved as a new document. 

%% 1 Load Data

clear all %clear the workspace
close all %close all figures
clc %clear the command window

load('Assignment1_data.mat') %Opens data file (assuming script and data are in the same directory)

%% 2 Initialization 
 startingNum = 0; %first trial
 endingNum = 300; %last trial
 stepSize = 50; %space between ticks on plot
 hair = 15; %breathing room spacing on plot
 fS = 14; %font size on plot

%% 3 Calculate Actual Reaction Times

% New variables represent the actual reaction time once the stimulus onset
% time is removed from the start time. Negative response times mean that
% trial was a false start. 

% 200mg Caffeine Group
RTcaff200 = (data_caffeine_200 - stimTime_caffeine_200); %subtract the stimulus onset time from the response time

% 400mg Caffeine Group
RTcaff400 = (data_caffeine_400 - stimTime_caffeine_400); %subtract the stimulus onset time from the response time

% Control Group
RTcontrol = (data_control - stimTime_control); %subtract the stimulus onset time from the response time


%% 4 Remove All Invalid (false start) Trials

%Initializations for loop below
validStart = 0; %making a start point for valid trials (response after the stimulus)
numTrials = 300; %number of trials
numGroups = 3; %number of groups (200, 400 & control)
kk = 0; %initalize counter
validTrials = zeros(numTrials, numGroups); %create 300x3 matrix to store output

for ii = [RTcaff200 RTcaff400 RTcontrol] %loop through each condition
    kk=kk+1; %increases by 1 each trial to avoid overwritting data
    ii(ii<validStart)=nan; %change all negative numbers to NaN's in response time
    validTrials(:,kk)= ii; %place each response time set in a separare column in the new Valid Trials matrix.
end
          
 %Within validTrials column 1 = 200mg, column 2 = 400mg, & column 3 = control


 %% 5 Plot Valid trials 
 
 %Make Plot
 Plot1 = figure; %Open a figure
 set(Plot1, 'Position', [100 100 800 800]) %set figure size
 RTPlot = plot(validTrials); %plot the valid trials for the three groups
 
 %label axis, title and legend
 xlabel('Trial Number') %x axis label
 ylabel('Reaction Time') %y axis label
 title('Effects of Caffeine on Reaction Time') %plot title
 legend('200mg Caffeine', '400mg Caffeine', 'Control') %add legend
 legend('Location', 'SouthEastOutside') %place legend outside the plot to the lower right
 
 %Physical atributes of the plot
 set(gca, 'XLim', [startingNum - hair endingNum + hair]) %add a hair of space on either end of the plot
 set(gca, 'XTick',(startingNum:stepSize:endingNum)) %set the tick spacing
 set(gca,'TickDir', 'out') %put the ticks outside the plot
 set(gca, 'Color', 'w'); set(Plot1, 'Color', 'w') %set plot and background color to white
 set(gca, 'FontSize', fS) %set fontsize to indexed size above
 set(RTPlot(3,1), 'Color', 'k') %set control group to black
 set(RTPlot(2,1), 'Color', 'r') %set 400mg group to red
 set(RTPlot(1,1), 'Color', 'b') %set 200mg group to blue
 box off %turn off the box around the plot
 
 shg %show figure 

%% 6 Save the Output

save('MatlabAssignment1') %saves .mat file to current directory. 

