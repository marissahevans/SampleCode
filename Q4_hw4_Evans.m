%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 4 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 4

% Date intialized: Oct. 31th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Analyzing and simulating experimental data. 
% An international coffee conglomerate recruits you to characterize the 
% neuropsychology underlying their customers? adoration of pumpkin spice. 
% You devise a blood-oxygen level dependent (BOLD) FMRI pilot experiment in 
% which you present one of two classes of odorants to an individual while monitoring
% the activity of three key voxels located in the amygdala, a structure known 
% to be associated with emotional responses. The file experimentData.mat contains: a (N × 3) matrix
% data, where each row is the BOLD response of the three voxels on a given trial relative to
% some baseline; and a (N × 1) vector trialConds indicating the experimental condition of
% each trial. Condition 1 are trials in which you present an odorant selected randomly from
% a library of possible control odorants, and condition 2 are trials in which the trade secret
% pumpkin spice odorant is presented.

load('experimentData.mat');
%% PART A

%% 
%* Before doing anything quantitative with your data, it is always good 
% practice to visualize it. First, determine how many trials of each trial condition were completed. 
% Display this information as a 2-bin histogram with each bin representing 
% each of the two possible trial conditions, and their heights representing their respective trial counts. 

figure
histogram(trialConds)
title(['Trial Conditions 1 and 2'])
xlabel('Condition #')
ylabel('# of Trials')
box off
set(gca, 'TickDir', 'out')
set(gca, 'XTick', [1,2])


%%
%* Next, plot a 3D scatter plot of the recorded responses, with each point color-coded according
% to its associated trial condition (use the function scatter3 in Matlab and be sure to
% label your axes). 

data2 = data;
data2(:,4) = trialConds;
data2 = sortrows(data2,4);

figure
scatter3(data2(1:300,1),data2(1:300,2),data2(1:300,3))
hold on
scatter3(data2(301:612,1),data2(301:612,2),data2(301:612,3))
title(['3D plot of data scatter'])
legend('Condition 1 - Random Scent', 'Condition 2 - Pumpkin Spice')
xlabel('voxel 1 BOLD')
ylabel('voxel 2 BOLD')
zlabel('voxel 3 BOLD')
axis equal

%%
%* Describe your data qualitatively using this figure. Is there a noticeable
% difference between the two trial conditions? What geometric shape are these ?response
% clouds?, and what distribution would you use to model them?

% There is a very noticible difference between the two trial conditions:
% Condition 2 elicits a strong response from all 3x voxels, while
% Condition 1 responses are very low across the board. All three voxels
% seem to be responding well so all dimentions are warrented. 

% The response clouds together are long and round like a baguette and the the responses
% are highly correlated, individually the clouds appear to be shaped
% similar to 3D gaussians. 

% To model the data a gaussion distribution would be the most sensible
% because the voxel activity presents with what appears to be a  multivariant 
% normal distribution. 

%% PART B

%% 
%* Quantify the response statistics of each individual trial condition. Calculate the means
% of each response cloud, as well as their respective covariance matrices. Compute the
% covariance matrices of each response cloud using matrix multiplication (remember to
% center the data first). 

D1 = data2(1:300,1:3)'; %trials in condition 1
N1 = length(D1);
D2 = data2(301:612,1:3)'; %trials in condition 2
N2 = length(D2);

%Mean
D1mean = mean(D1')
D2mean = mean(D2')

%Centering Data
D1center = D1' - D1mean;
D2center = D2' - D2mean;

%Covarianace Matrix
D1Cov = (D1center' * D1center) / [N1-1]
D2Cov = (D2center' * D2center) / [N2-1]


%%
%* Verify your calculation is correct by comparing with the output
% given by the cov function. How do the covariance matrices compare (are they similar
% at all or wildly different)?

covD1 = cov(D1center)
covD2 = cov(D2center)

%Comparing the Matlab covariance calculation to my mathmatical calculate
%they are identical. Comparing the differences between the two covariance
%matracies themselves they have a few notable elements:
%Condition 1 has higher variance in voxel 1 while condition 2 is the most
%varied in voxel 2. Voxel 3 has a farily stable variance across the
%conditions. While the total variance is similar the variance between the
%voxels seems to vary across the conditions. However in an experiement
%using thousands of voxels the variance between them would be negligible. 

%% PART C

%%
%* Next, compute the SVD of each covariance matrix. Plot the three singular 
% vectors originating from the center of each response cloud and scale their amplitude by the square
% root of the singular values. 

[U1 S1 V1] = svd(covD1);
[U2 S2 V2] = svd(covD2);

set1Center = D1mean;
set1Center(2,:) = D1mean;
set1Center(3,:) = D1mean;

set2Center = D2mean;
set2Center(2,:) = D2mean;
set2Center(3,:) = D2mean;

V1A = V1.*sqrt(S1);
V2A = V2.*sqrt(S2);

V1A = V1A + D1mean;
V2A = V2A + D2mean;

figure
P0 = set1Center;
P1 = V1A;
P2 = set2Center;
P3 = V2A;

X1 = [P0(:,1) P1(:,1)] ;
Y1 = [P0(:,2) P1(:,2)] ;
Z1 = [P0(:,3) P1(:,3)] ;

X2 = [P2(:,1) P3(:,1)] ;
Y2 = [P2(:,2) P3(:,2)] ;
Z2 = [P2(:,3) P3(:,3)] ;

X1 = X1' ;
Y1 = Y1' ;
Z1 = Z1' ;
X2 = X2' ;
Y2 = Y2' ;
Z2 = Z2' ;

scatter3(data2(1:300,1),data2(1:300,2),data2(1:300,3))
hold on
scatter3(data2(301:612,1),data2(301:612,2),data2(301:612,3))
plot3(X1(:,1),Y1(:,1),Z1(:,1),'b', 'LineWidth', 3)
plot3(X1(:,2),Y1(:,2),Z1(:,2),'b', 'LineWidth', 3,'HandleVisibility','off')
plot3(X1(:,3),Y1(:,3),Z1(:,3),'b', 'LineWidth', 3,'HandleVisibility','off')
plot3(X2(:,1),Y2(:,1),Z2(:,1),'r', 'LineWidth',3)
plot3(X2(:,2),Y2(:,2),Z2(:,2),'r', 'LineWidth',3,'HandleVisibility','off')
plot3(X2(:,3),Y2(:,3),Z2(:,3),'r', 'LineWidth',3,'HandleVisibility','off')
title('Singular Vectors of Conditions')
legend('Cond 1 - Random Scent', 'Cond 2 - Pumpkin Spice','Cond 1 Singular Vecs', 'Cond 2 Singluar Vecs')
grid on
xlabel('voxel 1 BOLD')
ylabel('voxel 2 BOLD')
zlabel('voxel 3 BOLD')
axis equal 

figure
plot3(X1,Y1,Z1,'b','LineWidth',2)
hold on
plot3(X2,Y2,Z2,'r', 'LineWidth',2)
title('Without data scatter to show variance angles')
grid on
axis equal
xlabel('voxel 1 BOLD')
ylabel('voxel 2 BOLD')
zlabel('voxel 3 BOLD')

%%
%* Relative to how similar the covariance matrices were before computing 
% their SVD, how do each trial condition?s respective set of singular values compare? 
% Describe what this tells us about the relationship between the two trial
% conditions and, more fundamentally, the relationship between the three voxels across
% conditions.

% The two conditions have very similar singular values with a difference of
% less than 1 for all dimentions. This tells us that the transformation
% applied to the data is relatively the same across conditons along the
% principle axis. Since the singluar vectors for each voxel are at differnt
% angles the variance of the data will be different when collapsed on to
% each axis individually, thus while the two conditions have similar
% variance overall the variance differs between the voxels within each
% condition. When considering a full application of this study the
% differences between the voxels themselves is negligible since thousands
% will be taken into account. The greater take away is the mean difference
% between the two conditions since ingroup variance is stable. If stacked
% on top of eachother the systems would be nearly identical with the
% variance confining and describing the entire system. 

%% PART D

%%
%* A powerful method to validate a model is by generating (i.e. simulating) new data
% matching your quantitative description of the real data, and then comparing them with
% real data. Create a function
% simResponses = odorExperiment(numTrials1,numTrials2) where numTrials1 and 
% numTrials2 are the number of trials in a simulated experiment for condition 
% 1 and 2, respectively. simResponses is a (N × 3) matrix containing simulated 
% responses of each of your 3 voxels during N = numTrials1 + numTrials2 trials. 
% [Hint: use ndRandn from the previous problem]. 

simResp = odorExperiment(300, 300);

%%
%* Plot the simulated and real 
% responses in the same figure (use subplots if you wish) to compare the two. 
% Is your simulated response data a good characterization of the real amygdala voxel responses?
figure
subplot(2,1,1)
scatter3(data2(1:300,1),data2(1:300,2),data2(1:300,3))
hold on
scatter3(data2(301:612,1),data2(301:612,2),data2(301:612,3))
title(['Real Neural Data'])
legend('Condition 1 - Random Scent', 'Condition 2 - Pumpkin Spice')
xlabel('voxel 1 BOLD')
ylabel('voxel 2 BOLD')
zlabel('voxel 3 BOLD')
axis equal 


subplot(2,1,2)
scatter3(simResp(1:300,1),simResp(1:300,2),simResp(1:300,3))
hold on
scatter3(simResp(301:600,1),simResp(301:600,2),simResp(301:600,3))
title(['Simulated Neural Data'])
legend('Condition 1 - Random Scent', 'Condition 2 - Pumpkin Spice')
xlabel('voxel 1 BOLD')
ylabel('voxel 2 BOLD')
zlabel('voxel 3 BOLD')
axis equal

% The simulated data is a good estimate of the real neural data by looking
% a the two plots visually. 

