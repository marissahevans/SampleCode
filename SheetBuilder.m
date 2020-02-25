%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TRIAL SORTING FOR  MGRS EYE AND REACH DATA %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Author: Marissa Evans; mhe229@nyu.edu
%Initalized: Oct. 17nd 2019

%Edited Feb. 5th 2020

% This script is written to sort trials eye tracking and reach mocap
% behavioral data as part of the MGRS (memory guided reach and saccade)
% project for Curtis Lab.

% Prior to running this script 'mgrs_analysis.m' must be run to create the
% files referenced below.

% We are creating a cell structure for each subject with relevant
% information from each trial which will be needed to run a MVPA analysis
% on the data.

%% INITALIZE THE WORKSPACE

clear all
close all
clc


subjParams.subj{1,1} = 's01';
subjParams.subj{2,1} = 's03';
subjParams.subj{3,1} = 's04';
subjParams.subj{4,1} = 's05';
subjParams.subj{5,1} = 's06';
subjParams.subj{6,1} = 's07';

subjParams.subj{1,2} = 's1_1_';
subjParams.subj{2,2} = 's3_1_';
subjParams.subj{3,2} = 's4_1_';
subjParams.subj{4,2} = 's5_1_';
subjParams.subj{5,2} = 's6_1_';
subjParams.subj{6,2} = 's7_1_';

subjParams.subj{1,3} = 's1_2_';
subjParams.subj{2,3} = 's3_2_';
subjParams.subj{3,3} = 's4_2_';
subjParams.subj{4,3} = 's5_2_';
subjParams.subj{5,3} = 's6_2_';
subjParams.subj{6,3} = 's7_2_';


%% Load data into the workspace

%select root file
mgrs_root = '/share/datb/mgrs';

for vv = 1:6
    
    
    %select subject to run
    subj = subjParams.subj{vv,1};
    
    %select file prefix
    filePref1 = subjParams.subj{vv,2};
    filePref2 = subjParams.subj{vv,3};
    
    
    %% Load data into the workspace
    
    %create stim structures
    sess1stim = cell(8,1);
    sess2stim = cell(8,1);
    sess1data = cell(8,1);
    sess2data = cell(8,1);
    
    numSet = (['1','2','3','4','5','6','7','8']);
    for jj = 1:8
        nn = numSet(jj);
        %Load sess 1 time series
        sess1stim{jj,:} = load(sprintf('%s/%s/stim/sess_1/timeSeries/%s_1_stim0%s', mgrs_root, subj,subj,nn));
        
    end
    
    %Filling in missing block type info for s01 session 1
    if vv == 1
        sess1stim{1,1}.task{1,1}.blocktype = 1;
        sess1stim{2,1}.task{1,1}.blocktype = 1;
        sess1stim{3,1}.task{1,1}.blocktype = 0;
        sess1stim{4,1}.task{1,1}.blocktype = 0;
        sess1stim{5,1}.task{1,1}.blocktype = 1;
        sess1stim{6,1}.task{1,1}.blocktype = 1;
        sess1stim{7,1}.task{1,1}.blocktype = 0;
        sess1stim{8,1}.task{1,1}.blocktype = 0;
    end
    
    %Fixing trial label issues for s07 session 1
    if vv == 6
        sess1stim{4,1}.task{1,1}.run = 4;
        sess1stim{5,1}.task{1,1}.run = 5;
        sess1stim{6,1}.task{1,1}.run = 6;
        sess1stim{7,1}.task{1,1}.run = 7;
        sess1stim{7,1}.task{1,1}.session = 1;
        sess1stim{8,1}.task{1,1}.run = 8;
        sess1stim{8,1}.task{1,1}.session = 1;
        
    end
    
    for jj = 1:8
        nn = numSet(jj);
        %Load sess 1 time series
        sess2stim{jj,:} = load(sprintf('%s/%s/stim/sess_2/timeSeries/%s_2_stim0%s', mgrs_root, subj,subj,nn));
        
    end
    
    for jj = 1:8
        nn = numSet(jj);
        if exist(sprintf('%s/%s/eyefiles/sess_1/%s%s_preproc.mat', mgrs_root, subj,filePref1,nn));
            %Load time series
            sess1data{jj,:} = load(sprintf('%s/%s/eyefiles/sess_1/%s%s_preproc.mat', mgrs_root, subj,filePref1,nn));
        else
            continue
        end
    end
    
    for jj = 1:8
        nn = numSet(jj);
        if exist(sprintf('%s/%s/eyefiles/sess_2/%s%s_preproc.mat', mgrs_root, subj,filePref2,nn));
            %Load time series
            sess2data{jj,:} = load(sprintf('%s/%s/eyefiles/sess_2/%s%s_preproc.mat', mgrs_root, subj,filePref2,nn));
        else
            continue
        end
    end
    
    %% ADDING TRIAL INFORMATION INTO THE SHEET
    
    trialGuideMat = nan*ones(192,17);
    
    % Columns:
    %1. SessNum - 1 or 2
    %2. RunNum - 1 through 8
    %3. TrialNum - total trial num across both sessions
    %4. RunType - 0 or 1; 0 = saccade ;  1 = reach
    %5. TrialType - 0 or 1; 0 = catch; 1= trial
    %6. TarX - x coordinate of target
    %7. TarY - y coordinate of target
    %8. TarQuadrant - 1,2,3 or 4 1 = upper right, 2=upper left, 3=lower left,
    %4 = lower right
    %9. GoCue Time Point
    %10. Fixation Time
    %11. Target Time
    %12. Delay Time
    %13. Saccade Time
    %14. Saccade Correction Time
    %15. Fixation Time
    %16. ITI Time
    %17. Total Trial Time
   
    
    % ADDING TRIAL INFORMATION
    
    %Session 1
    cc = 1;
    for ii = 1:8
        
        for jj = cc:cc+11
            trialGuideMat(jj,1) = sess1stim{ii,1}.task{1,1}.session;
            trialGuideMat(jj,2) = sess1stim{ii,1}.task{1,1}.run;
            trialGuideMat(jj,3) = jj;
            trialGuideMat(jj,4) = sess1stim{ii,1}.task{1,1}.blocktype;
            trialGuideMat(jj,10) = 1.775; %fixation
            trialGuideMat(jj,11) = .225; %target
            trialGuideMat(jj,13) = 2.5; %saccade
            trialGuideMat(jj,14) = 2.5; %correction
            trialGuideMat(jj,15) = 1; %fixation
        end
        for pp = 0:11
            trialGuideMat(jj-pp,5) = sess1stim{ii,1}.task{1,1}.randVars.trialType(12-pp);
            trialGuideMat(jj-pp,16) = sess1stim{ii,1}.task{1,1}.randVars.ITI(12-pp);    
            trialGuideMat(jj-pp,12) = sess1stim{ii,1}.task{1,1}.randVars.delayPeriod(12-pp);
        end
        cc = cc + 12;
        
    end
    
    %Session 2
    cc = 97;
    for ii = 1:8
        
        for jj = cc:cc+11
            trialGuideMat(jj,1) = sess2stim{ii,1}.task{1,1}.session;
            trialGuideMat(jj,2) = sess2stim{ii,1}.task{1,1}.run;
            trialGuideMat(jj,3) = jj;
            trialGuideMat(jj,4) = sess2stim{ii,1}.task{1,1}.blocktype;
            trialGuideMat(jj,10) = 1.775; %fixation
            trialGuideMat(jj,11) = .225; %target
            trialGuideMat(jj,13) = 2.5; %saccade
            trialGuideMat(jj,14) = 2.5; %correction
            trialGuideMat(jj,15) = 1; %fixation

        end
        for pp = 0:11
            trialGuideMat(jj-pp,5) = sess2stim{ii,1}.task{1,1}.randVars.trialType(12-pp);
            trialGuideMat(jj-pp,12) = sess2stim{ii,1}.task{1,1}.randVars.delayPeriod(12-pp);
            trialGuideMat(jj-pp,16) = sess2stim{ii,1}.task{1,1}.randVars.ITI(12-pp);
            
        end
        cc = cc + 12;
    end
    
    for tt = 1:length(trialGuideMat)
        trialGuideMat(tt,17) = sum(trialGuideMat(tt,10:16)) %trial length
        if trialGuideMat(tt,5) == 0
            trialGuideMat(tt,13:15) = 0
        else
        end
    end
    
    %Adding this missing go-cue timing for subject 1
    if vv == 1 
        trialGuideMat(1,9) = sum(trialGuideMat(1,10:12));
        for sessii = 1:2
            for runii = 1:8
                if sessii == 2
                    runii = runii + 8;
                end
                for trialii = 1:12
                    if runii > 1
                        if trialii == 1
                            trialGuideMat(((runii-1)*12)+trialii,9) = sum(trialGuideMat(((runii-1)*12)+trialii,10:12));
                        else
                            trialGuideMat(((runii-1)*12)+trialii,9) = trialGuideMat(((runii-1)*12)+(trialii-1),9) + sum(trialGuideMat(((runii-1)*12)+(trialii-1),13:16)) + sum(trialGuideMat(((runii-1)*12)+trialii,10:12));
                            
                        end
                    else
                        trialGuideMat((runii)+trialii,9) = trialGuideMat((runii)+(trialii-1),9) + sum(trialGuideMat((runii)+(trialii-1),13:16)) + sum(trialGuideMat((runii)+trialii,10:12));
                        
                    end
                end
            end
        end
    end
    %% Finding the target time points
    
    %SESSION 1
    % ADDING TAR-X AND TAR-Y w/ QUADRANTS
    
    cc = 1;
    for ii = 1:8
        targetsXFin2 = nan*ones(10,1);
        targetsYFin2 = nan*ones(10,1);
        if isempty(sess1data{ii,1})
            cc = cc + 12;
            continue
        else
            %Finding all the target points for X
            dd = 1;
            for ww = 2:length(sess1data{ii,1}.ii_data.TarX)
                
                if sess1data{ii,1}.ii_data.TarX(ww) - sess1data{ii,1}.ii_data.TarX(ww-1) ~0
                    targetsX(dd) = ww;
                    dd = dd+1;
                else
                    continue
                end
            end
            
            %Finding all the target points for Y
            
            dd = 1;
            for ww = 2:length(sess1data{ii,1}.ii_data.TarY)
                
                if ((sess1data{ii,1}.ii_data.TarY(ww)) - (sess1data{ii,1}.ii_data.TarY(ww-1))) ~0
                    targetsY(dd) = ww;
                    dd = dd+1;
                else
                    continue
                end
            end
            
            %TarX values
            targetsXFin = sess1data{ii,1}.ii_data.TarX(targetsX);
            index = find(targetsXFin ~=0);
            sample = targetsXFin(index);
            for yy = 1:length(index)
            targetsXFin2(yy,1) = sample(yy);
            end
            
            %TarY values
            targetsYFin = sess1data{ii,1}.ii_data.TarY(targetsY);
            index = find(targetsYFin ~=0);
            sample = targetsYFin(index);
            for yy = 1:length(index)
            targetsYFin2(yy,1) = sample(yy);
            end
            
            uu = 1;
            
            for jj = cc:cc+11
                if trialGuideMat(jj,5) == 1
                    if isnan(targetsXFin2(uu))
                        uu = uu+1;
                        continue
                    else
                        trialGuideMat(jj,6) = targetsXFin2(uu);
                        trialGuideMat(jj,7) = targetsYFin2(uu);
                        uu = uu+1;
                    end
                else
                    trialGuideMat(jj,6) = 0;
                    trialGuideMat(jj,7) = 0;
                end
                
                if trialGuideMat(jj,6)>0 && trialGuideMat(jj,7)>0
                    trialGuideMat(jj,8) = 1;
                elseif trialGuideMat(jj,6)<0 && trialGuideMat(jj,7)>0
                    trialGuideMat(jj,8) = 2;
                elseif trialGuideMat(jj,6)<0 && trialGuideMat(jj,7)<0
                    trialGuideMat(jj,8) = 3;
                elseif trialGuideMat(jj,6)>0 && trialGuideMat(jj,7)<0
                    trialGuideMat(jj,8) = 4;
                else
                    trialGuideMat(jj,8) = 0;
                end
                
            end
        end
        cc = cc + 12;
        clear targetsX
        clear targetsY
    end
    
    %SESSION 2
    
    % ADDING TAR-X AND TAR-Y w/ QUADRANTS
    cc = 97;
    for ii = 1:8
        if isempty(sess2data{ii,1})
            cc = cc + 12;
            continue
        else
            %Finding all the target points for X
            dd = 1;
            for ww = 2:length(sess2data{ii,1}.ii_data.TarX)
                
                if sess2data{ii,1}.ii_data.TarX(ww) - sess2data{ii,1}.ii_data.TarX(ww-1) ~0
                    targetsX(dd) = ww;
                    dd = dd+1;
                else
                    continue
                end
            end
            
            %Finding all the target points for Y
            
            dd = 1;
            for ww = 2:length(sess2data{ii,1}.ii_data.TarY)
                
                if ((sess2data{ii,1}.ii_data.TarY(ww)) - (sess2data{ii,1}.ii_data.TarY(ww-1))) ~0
                    targetsY(dd) = ww;
                    dd = dd+1;
                else
                    continue
                end
            end
            
            %TarX values
            targetsXFin = sess2data{ii,1}.ii_data.TarX(targetsX);
            index = find(targetsXFin ~=0);
            sample = targetsXFin(index);
            for yy = 1:length(index)
            targetsXFin2(yy,1) = sample(yy);
            end
            
            %TarY values
            targetsYFin = sess2data{ii,1}.ii_data.TarY(targetsY);
            index = find(targetsYFin ~=0);
            sample = targetsYFin(index);
            for yy = 1:length(index)
            targetsYFin2(yy,1) = sample(yy);
            end
            
            
            uu = 1;
            
            for jj = cc:cc+11
                if trialGuideMat(jj,5) == 1
                    if isnan(targetsXFin2(uu))
                        uu = uu+1;
                        continue
                    else
                    trialGuideMat(jj,6) = targetsXFin2(uu);
                    trialGuideMat(jj,7) = targetsYFin2(uu);
                    uu = uu+1;
                    end
                else
                    trialGuideMat(jj,6) = 0;
                    trialGuideMat(jj,7) = 0;
                end
                
                if trialGuideMat(jj,6)>0 && trialGuideMat(jj,7)>0
                    trialGuideMat(jj,8) = 1;
                elseif trialGuideMat(jj,6)<0 && trialGuideMat(jj,7)>0
                    trialGuideMat(jj,8) = 2;
                elseif trialGuideMat(jj,6)<0 && trialGuideMat(jj,7)<0
                    trialGuideMat(jj,8) = 3;
                elseif trialGuideMat(jj,6)>0 && trialGuideMat(jj,7)<0
                    trialGuideMat(jj,8) = 4;
                else
                    trialGuideMat(jj,8) = 0;
                end
                
            end
        end
        cc = cc + 12;
        clear targetsX
        clear targetsY
    end
    
    
    %% Add Time Point for Go Cue
    
    %check if time points have been entered
    if exist(sprintf('%s/%s/stim/sess_1/1SMovL.txt',mgrs_root,subj));
        
    %Session 1
    eyeGoCueLeft1 = sprintf('%s/%s/stim/sess_1/1SMovL.txt',mgrs_root,subj);
    eyeGoCueRight1 = sprintf('%s/%s/stim/sess_1/1SMovR.txt',mgrs_root,subj);
    reachGoCueLeft1 = sprintf('%s/%s/stim/sess_1/1RMovL.txt',mgrs_root,subj);
    reachGoCueRight1 = sprintf('%s/%s/stim/sess_1/1RMovR.txt',mgrs_root,subj);
    
    %Concatonate the reach trials for session 1
    fid = fopen(reachGoCueLeft1);
    s1 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
    fclose(fid);
    
    fid = fopen(reachGoCueRight1);
    s2 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
    fclose(fid);
    
    reachTime = [s1{1,1},s2{1,1}];
    for tt = 1:8
        reachTime(tt,:) = sort(reachTime(tt,:));
    end
    
    %Concatonate the eye trials for session 1
    fid = fopen(eyeGoCueLeft1);
    s1 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
    fclose(fid);
    
    fid = fopen(eyeGoCueRight1);
    s2 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
    fclose(fid);
    
    eyeTime = [s1{1,1},s2{1,1}];
    
    for tt = 1:8
        eyeTime(tt,:) = sort(eyeTime(tt,:));
    end
    
    %Put the times in the sheet
    cc = 1;
    for ii = 1:8
        
        for jj = cc:cc+11
            if trialGuideMat(jj,4) == 0
                trialGuideMat(jj,9) = eyeTime(ii,jj+1-cc);
                
            elseif trialGuideMat(jj,4) == 1
                trialGuideMat(jj,9) = reachTime(ii,jj+1-cc);
                
            end
            
        end
        cc = cc + 12;
    end
    
    % Session 2
    eyeGoCueLeft2 = sprintf('%s/%s/stim/sess_2/2SMovL.txt',mgrs_root,subj);
    eyeGoCueRight2 = sprintf('%s/%s/stim/sess_2/2SMovR.txt',mgrs_root,subj);
    reachGoCueLeft2 = sprintf('%s/%s/stim/sess_2/2RMovL.txt',mgrs_root,subj);
    reachGoCueRight2 = sprintf('%s/%s/stim/sess_2/2RMovR.txt',mgrs_root,subj);
    
    %Concatonate the reach trials for session 2
    fid = fopen(reachGoCueLeft2);
    s1 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
    fclose(fid);
    
    fid = fopen(reachGoCueRight2);
    s2 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
    fclose(fid);
    
    reachTime = [s1{1,1},s2{1,1}];
    for tt = 1:8
        reachTime(tt,:) = sort(reachTime(tt,:));
    end
    
    %Concatonate the eye trials for session 2
    fid = fopen(eyeGoCueLeft2);
    s1 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
    fclose(fid);
    
    fid = fopen(eyeGoCueRight2);
    s2 = textscan(fid,'%f %f %f %f %f %f','CollectOutput',1,'treatAsEmpty','* 	*');
    fclose(fid);
    
    eyeTime = [s1{1,1},s2{1,1}];
    
    for tt = 1:8
        eyeTime(tt,:) = sort(eyeTime(tt,:));
    end
    
    %Put the times in the sheet
    cc = 97;
    for ii = 1:8
        
        for jj = cc:cc+11
            if trialGuideMat(jj,4) == 0
                trialGuideMat(jj,9) = eyeTime(ii,jj+1-cc);
                
            elseif trialGuideMat(jj,4) == 1
                trialGuideMat(jj,9) = reachTime(ii,jj+1-cc);
                
            end
            
        end
        cc = cc + 12;
    end
    else
        %continue
    end
    %% save output
    
    save(sprintf('%s/ScoredData/%s_%s', mgrs_root, subj, 'trialGuideMat'), 'trialGuideMat')
end