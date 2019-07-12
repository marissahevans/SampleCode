%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Date Created: June 21th 2019
% Date Last Edit: June 24th 2019
% Author: Marissa Evans, mhe229@nyu.edu

% Assumptions: The file 'studentGradesAdvancedStats.mat' is in the same 
% directory as the code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script compares the scores of 100 statistics students on a midterm
% and final exam. Some students did not attend all recitations and some
% students were allowed to use their cell phones during class. The analysis
% looks at how these different parameters relate and how they effect the
% student's final grade in the course. 

% What this code does(by section):
% 0) Initialization
% 1) Load the file 'studentGradesAdvancedStats.mat'
% 2) Clean the data and remove NaN's
% 3) Determine if there is a significant correlation between the midterm
% and final exam 
% 4) Plot the midterm vs final exam performance and add regression line -
% discuss how what these outcomes mean. 
% 5) Creates variable 'gradeScore' and run t-tests to see which variables 
% effected gradeScore
% 6) Run a full model 'anovan' (2 group) to see how cell phone use and
% recitation attendance effect gradeScore and what their interaction is
% 7) Calculate the mean diff. between midterm and final, then bootstrap to
% see if means are stable and calculate the 95% CI. 

%% 0 INITALIZATION

clear all
close all
clc

fS = 20; %font size for figures
lW = 3; %line width for figures
numReps = 1e5; %100k repetiions 
numPop = 90; %total number of students after NaN's removed
maxSample = 90; %Sample size drawn from population
%% 1 LOAD DATA FILES

load('studentGradesAdvancedStats.mat')

%% 2 CLEAN AND ORGANIZE DATA 

% Removing NaN's partipant wise, only looking at data from students who did
% not miss any exams. Since only 10x students missed an exam this still
% leaves us with a large subject pool but allows for repeated measures
% within the students themselves. 

DATA = ones(100, 5)*999; %initalize matrix with 999's

DATA(:,1) = attendedRecitation;
DATA(:,2) = usedCell;
DATA(:,3) = midtermScore;
DATA(:,4) = finalScore;

temp = sum(isnan(DATA),2); %counts missing values per row (the number is the DIMENSION!!) 1 = columns, 2=rows
temp2 = find(temp==0); %only keep rows of students that took both exams
DATA = DATA(temp2,:); %rewrite DATA with only scores of studetns who took both exams 

%% 3 CALCULATE CORRELEATION BETWEEN MIDTERM AND FINAL 

examCorr = corr(DATA(:,3), DATA(:,4));

% The two exams are highly and significanly correlated: r = .6, p=.000

%% 4 PLOT MIDTERM vs FINAL PERFORMANCE

% There is a linear trend in comparing the midterm and the final exam
% scores, students who scored higher on the midterm tended to score higher
% on the final as well. 

x=DATA(:,3); %x-axis will be the midterm Exam
y=DATA(:,4); %y-axis will be the final Exam

%Make Plot
examFig = figure; %Open a figure
set(examFig, 'Position', [100 100 800 800]); %set figure size
respPlot = scatter(x,y, 'filled'); %plot a scatterplot of the midterm exams vs the final exams by student
examRegLine = lsline;

%label axis, title and legend
xlabel('Midterm Exam Score'); %x axis label
ylabel('Final Exam Score'); %y axis label
title('Stats Midterm Exam vs Final Exam Scores'); %plot title

%Physical atributes of the plot
examRegLine.Color = 'r';
examRegLine.LineWidth = lW;
set(gca,'TickDir', 'out') %put the ticks outside the plot
set(gca, 'Color', 'w'); set(examFig, 'Color', 'w'); %set plot and background color to white
set(gca, 'FontSize', fS); %set fontsize to indexed size above
box off %turn off the box around the plot

shg %show figure

%% 5 T-TESTS TO SEE WHAT VARIABLES EFFECT FINAL GRADE

gradeScore = (DATA(:,3) + DATA(:,4))/2; %create variable 'gradeScore' which is the average of midterm and final grades
DATA(:,5) = gradeScore; %put the final grade in the 5th column of DATA


% a) Is there an effect of cell phone use on total grade score? 
cellUsed = find(DATA(:,2) == 1); %find all rows where students used cell phone
noCellUsed = find(DATA(:,2) == 0); %find all rows where studented didnt use cell phone

cellUsed = DATA(cellUsed,5); %group the final grades for students who used cell phone
noCellUsed = DATA(noCellUsed,5); %group the final grades for students who didnt use cell phone

[sig, p, CI, tstats] = ttest2(cellUsed, noCellUsed)
        % Yes, there was a significant effect of cell phone use on the the
        % final score; p = < .001, tstat: 7.6374, df: 88, sd: 12.3639, 
        % CI: 14.72 - 25.09
        % Students who did NOT use their cell phone performed statistically better in the
        % course
        
% b) Is there an effect of recitation attendance on total grade score? 
attendedRes = find(DATA(:,1) == 1); %find all rows where students attended classs
notAttended = find(DATA(:,1) == 0); %find all rows where studented didnt attend class

attendedRes = DATA(attendedRes,5); %group the final grades for students who attended class
notAttended = DATA(notAttended,5); %group the final grades for students who didnt attend class

[sig, p, CI, tstats] = ttest2(attendedRes, notAttended)
        % Yes, there was a significant effect of recitation attendance on the the
        % final score; p < .001, tstat: 4.8073, df: 88, sd: 14.1868, 
        % CI: 8.4341 - 20.3215
        % Students who attended recitation recieved statistically higher scores in the class.  

% c) Was one exam harder than the other one?
[sig, p, CI, tstats] = ttest(DATA(:,3), DATA(:,4))
        % Yes, the final exam (mean = 57.43) was significantly more difficult 
        % than the midterm exam (mean = 65.85); p < .001, tstat: 5.0536, df: 89, sd: 15.8107

%% 6 RUN A TWO WAY ANOVA TO SEE THE INTERACTIONS RESPONSIBLE IN THE FINAL GRADE

% n-Group ANOVA comparing how cell use and recitation attendance effect the
% final grade score, and how they interact. 
[p, anovatab] = anovan(DATA(:,5),[DATA(:,1),DATA(:,2)], 'model','full',...
    'varnames',{'attendedRes','cellUse'})

% While both cell phone use(F=18.23, p=.000) and recitation attendace(F=51.17, p=.000) 
% strongly effect the final grade score, they do not have a significant interaction(F=.16, p=.68). 
% This means that cell use and lack of attendance both negatively effected
% the final grade score, however they did not work together exponentially
% to lower the grade even more when someone did both. 

%% 7 BOOTSTRAPPING DATA WITH CI

%find the actual mean difference bwtween the midterm and final score
empiricalMeanDiff = mean(DATA(:,3)) - mean(DATA(:,4));
% mean diff = 8.42

resampledMeans = nan(numReps,2); %preallocate for all resamples means, 1 column per exam, each run is a row
n1 = length(DATA(:,3)); %length of midterm exam results (90)
n2 = length(DATA(:,4)); %length of final exam results (90)
for ii = 1:numReps %go though loop 100k times
    index = randi(n1, [n1, 1]); %make an index of rows to subsample exam scores from (with replacement)
    resampledMidterm = DATA(index,3); %take samples from the 3rd column of DATA (midterm exam)
    resampledMeans(ii,1) = mean(resampledMidterm); %put the mean of the samples drawn in the 1st column of a matrix
    
    index = randi(n2, [n2, 1]); %make an index of rows to subsample exam scores from (with replacement)
    resampledFinal = DATA(index,4); %take samples from the 4th column of DATA (final exam)
    resampledMeans(ii,2) = mean(resampledFinal); %put the mean of the samples drawn in the 2nd column of a matrix
end

% PLOTTING WITH CI
meanDifferences = resampledMeans(:,1) - resampledMeans(:,2); %find the mean differences between the 100k resampled means
meanDifferences = sortrows(meanDifferences); %rewrite variable with rows in ascending order
highCI = meanDifferences(length(meanDifferences)*.025); %create a variable that marks the row of the lowest 2.5% of scores 
lowCI = meanDifferences(length(meanDifferences)*.975); %create a variable that marks the row of the highest 2.5% of scores

bootFig = figure; %open a figure
set(bootFig, 'Position', [100 100 800 800]); %set figure size

%Put data on figure
handleAll = histogram(meanDifferences,100); %create a histrogram with 100 bins
handleAll.FaceAlpha = .8; %make histogram bars slightly translucent
handleAll.FaceColor = 'k'; %make face color back
hold on;

%Add CI lines to the histogram
line([lowCI, lowCI], ylim, 'LineWidth', 2, 'Color', 'r'); %add a stright line at the lower bound of the 95% CI in red
line([highCI, highCI], ylim, 'LineWidth', 2, 'Color', 'r'); %add a stright line at the higher bound of the 95% CI in red

%label axis, title and legend
xlabel('Mean Difference'); %x axis label
ylabel('Number of Samples'); %y axis label
title('Bootstrapping Mean Differences in Exam Scores'); %plot title

%Physical atributes of the figure
set(gca,'TickDir', 'out'); %put the ticks outside the plot
set(gca, 'Color', 'w'); set(bootFig, 'Color', 'w'); %set plot and background color to white
set(gca, 'FontSize', fS); %set fontsize to indexed size above
box off %turn off the box around the plot
shg %show figure

% Because our calcuated emperical mean difference was centered within the
% the 95% CI for the bootstrapped data it is reasonable to assume that this
% mean difference is statistically reliable. Additionally since 0 does not
% fall within the 95% CI we can strongly assume that the two exam score
% grades are statistically different. 