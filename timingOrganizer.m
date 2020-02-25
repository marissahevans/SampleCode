%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TRIAL TIMING EXTRACTION FOR DELAYS, TARGET, & GO CUE %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Author: Marissa Evans; mhe229@nyu.edu
%Initalized: Feburary 25th, 202

% This script is written to extract timing for delay periods to use for
% MVPA analysis as part of the MGRS (memory guided reach and saccade)
% project for Curtis Lab.

% This scipt pulls from timing files created using 'afniMatrix.m' and
% 'SheetBuilder.m', which in turn have their own source files needed to
% create the data they handle. 

% Running this script creates a text file for each trial organized by it's
% delay period in order to be input into the 'do_glm_new' model shell
% script for each participant. 

%% INITIALIZE SCRIPT
clear all 

%select root file
mgrs_root = '/share/datb/mgrs';
subj = 's07';
myOutputDirectory = '/share/datb/mgrs/s07/models'

%% Targets
STarL = sprintf('%s/%s/models/STarL-1.txt',mgrs_root,subj);
STarR = sprintf('%s/%s/models/STarR-1.txt',mgrs_root,subj);
RTarL = sprintf('%s/%s/models/RTarL-1.txt',mgrs_root,subj);
RTarR = sprintf('%s/%s/models/RTarR-1.txt',mgrs_root,subj);

%Concatonate the target files
fid = fopen(STarL);
STarL1 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
fclose(fid);

fid = fopen(STarR);
STarR1 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
fclose(fid);

fid = fopen(RTarR);
RTarR1 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
fclose(fid);

fid = fopen(RTarL);
RTarL1 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
fclose(fid);

STar = sort([STarL1{1,1},STarR1{1,1}],2);
RTar = sort([RTarL1{1,1},RTarR1{1,1}],2);

dlmwrite(fullfile(myOutputDirectory, 'STar.txt'),STar,'delimiter','\t');
dlmwrite(fullfile(myOutputDirectory, 'RTar.txt'),RTar,'delimiter','\t');

clear STarL STarR RTarR RTarL STarL1 STarR1 RTarR1 RTarL1

%% Movements 

SMovL = sprintf('%s/%s/models/SMovL.txt',mgrs_root,subj);
SMovR = sprintf('%s/%s/models/SMovR.txt',mgrs_root,subj);
RMovL = sprintf('%s/%s/models/RMovL.txt',mgrs_root,subj);
RMovR = sprintf('%s/%s/models/RMovR.txt',mgrs_root,subj);

%Concatonate the movement files
fid = fopen(SMovL);
SMovL1 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
fclose(fid);

fid = fopen(SMovR);
SMovR1 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
fclose(fid);

fid = fopen(RMovR);
RMovR1 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
fclose(fid);

fid = fopen(RMovL);
RMovL1 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
fclose(fid);

SMov = sort([SMovL1{1,1},SMovR1{1,1}],2);
RMov = sort([RMovL1{1,1},RMovR1{1,1}],2);

dlmwrite(fullfile(myOutputDirectory, 'SMov.txt'),SMov,'delimiter','\t');
dlmwrite(fullfile(myOutputDirectory, 'RMov.txt'),RMov,'delimiter','\t');

clear SMovL SMovR RMovR RMovL SMovL1 SMovR1 RMovR1 RMovL1

%% Delays 

SDel = sprintf('%s/%s/models/SDel.txt',mgrs_root,subj);
RDel = sprintf('%s/%s/models/RDel.txt',mgrs_root,subj);

fid = fopen(SDel);
SDel1 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
fclose(fid);

fid = fopen(RDel);
RDel1 = textscan(fid,'%f %f %f %f %f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
fclose(fid);

clear SDel RDel

SDel = sort(SDel1{1,1},2);
RDel = sort(RDel1{1,1},2);

clear RDel1 SDel1

%% Trial Labels

load(sprintf('%s/ScoredData/%s_%s', mgrs_root, subj,'trialGuideMat'));

trialGuideMat(97:192,2) = trialGuideMat(97:192,2) + 8;

%Isolate saccade / reach trials 
sacc = trialGuideMat(:,4) == 0;
reach = trialGuideMat(:,4) == 1;

%shorten matrix to only sacc or reach
saccTrials = trialGuideMat(sacc,:);
reachTrials = trialGuideMat(reach,:);

%find catch trials to ignore
sacCatch = saccTrials(:,5) == 1;
reachCatch = reachTrials(:,5) == 1;

%Remove catch trials 
saccTrials = saccTrials(sacCatch,:);
reachTrials = reachTrials(reachCatch,:);

%This creates a matrix with the delay times for each trial for either all
%the reach or saccade trials, which can be input into the glm file for
%each subject. 
saccRunDelay(:,1) = round(saccTrials(:,12));
saccRunDelay(:,2) = saccTrials(:,2);

reachRunDelay(:,1) = round(reachTrials(:,12));
reachRunDelay(:,2) = reachTrials(:,2);

%% Sorting Delay into .txt files 

SDel = SDel';
RDel = RDel';

dlmwrite(fullfile(myOutputDirectory, 'SDel.txt'),SDel,'delimiter','\t');
dlmwrite(fullfile(myOutputDirectory, 'RDel.txt'),RDel,'delimiter','\t');

%Creating individual files for each trial
x = ~isnan(SDel);
SDel2 = SDel(x);

x = ~isnan(RDel);
RDel2 = RDel(x);

samp = (1:1:80)';

for ii = 1:9
    file1 = RDel2(ii);
    trial = num2str(ii);
    dlmwrite(fullfile(myOutputDirectory,(sprintf('RDelT0%s.txt', trial))),file1,'delimiter','\t');
    
    file2 = SDel2(ii);
    dlmwrite(fullfile(myOutputDirectory,(sprintf('SDelT0%s.txt', trial))),file2,'delimiter','\t');
    
end

for ii = 10:80
file1 = RDel2(ii);
    trial = num2str(ii);
    dlmwrite(fullfile(myOutputDirectory,(sprintf('RDelT%s.txt', trial))),file1,'delimiter','\t');
    
    file2 = SDel2(ii);
    dlmwrite(fullfile(myOutputDirectory,(sprintf('SDelT%s.txt', trial))),file2,'delimiter','\t');
    
end