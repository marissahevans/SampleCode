%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 5 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 3

% Date intialized: Nov. 20th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Psychopathy. 
% You are interested in causes and treatment options for Psychopathy. You
% obtained a dataset, contained in the file psychopathy.mat obtained from a prison for violent
% offenders in upstate New York (not everyone in the prison is a psychopath, but they are
% more prevalent than in the general population). Each row represents data from one prisoner.
% All study participants underwent a structural scan with a mobile, truck-mounted MRI. The
% first data column contains the estimated cortical volume of paralimbic areas, relative to the
% population median, in cm3. The second column contains the Hare Psychopathy Checklist
% (PCL-R) scores, which range from 0 to 40 (the higher the score, the more psychopathic
% traits someone exhibits). These scores are not distributed normally in the general population
% (median = 4) and definitely not normal in this subpopulation (median = 20). The third
% column indicates whether they already participated in an experimental treatment program
% known as 'decompression therapy' (0 = did not yet participate, 1 = did already participate).
% To avoid self-selection effects, everyone in this dataset agreed to the therapy, but prisoners
% were randomly assigned to an earlier and a later treatment group, so that the untreated
% prisoners could serve as a control group

load('psychopathy.mat')

%% PART A

%%
%* You want to model PCL-R scores as a function of relative volume of paralimbic areas. Use
% polynomial regression to find a model that best explains the data using cross-validation.
% What degree does it have? (Note, you can use your code from HW2.)

DATA2 = sortrows(DATA);

y = DATA2(:,2);
x = DATA2(:,1);

N = length(DATA);

%Visualizing the data to start out
figure
plot(x,y, 'ro')
axis equal
axis square
xline(0, '--');
title('PCL-R as a function of Paralimbic Area')
legend('data')
box off
set(gca, 'TickDir', 'out')
xlabel('X')
ylabel('Y')

% Using a 'leave one out' cross validation for the regression 

numPoints = size(y,1);
XX = [ones(numPoints,1), x, x.^2,x.^3,x.^4,x.^5,x.^6,x.^7];
numModels = size(XX,2);

numFolds = numPoints; 

mseXval = zeros(numFolds, numModels);

figure('Position', [10 10 1200 800]);
for ii = 1:numModels
    for jj = 1:numPoints
        
        xInd = true(numPoints,1);
        xInd(jj)=false;
        
        trainX = XX(xInd,1:ii);
        trainY = y(xInd); 
        
        testX = XX(~xInd,1:ii);
        testY = y(~xInd);
        
        
        [U,~,V] = svd(trainX);
        S_vec = svd(trainX); 
        S_inv = zeros(size(trainX));
        S_inv(1:ii,1:ii) = diag(1./S_vec);
        betaVal = V*S_inv'*U'*trainY; 
        
        fitTrain =  XX(:,1:ii)*betaVal; 
        fitTest = testX * betaVal; 
        
        mseXval(jj,ii) = mean((testY - fitTest).^2); 
        
        subplot(2,4,ii)
        plot(x,y, 'o')
        hold on
        plot(x,fitTrain, 'LineWidth', 2)
        title(['Polynomials at Order ' num2str(ii-1)])
        xlabel('X');
        ylabel('Y');
        set(gca, 'TickDir', 'out')
        box off
        hold on
    end
end

%%
figure
errorbar(0:7, mean(mseXval), std(mseXval)/sqrt(numPoints), 'LineWidth', 2)
title(['Squared error as a function of polynomial order with crossvalidation'])
legend('Error')
xlabel('Order Polynomials');
ylabel('Squared Error');
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 14)
box off
hold on

[minMSE, minInd] = min(mean(mseXval));
threshMSE = minMSE + std(mseXval(:,minInd))/sqrt(numPoints);

plot(0:7, ones(numModels,1)* threshMSE ,'k--')
legend('MSE and 95% conf. int.', 'Best Fit Thresh')
set(gca,'xtick',0:7)

%The degree of the polynomial that explains the greatest amount of error
%without overfitting the data would be order 1. Looking at the error bars
%it would not make sense to use any higher order polynomial because the
%variance (2x SEM in this case) is not outside of the error still present 
%at order 1. The fit of the higher order polynomials starts to increase in
%error because as we continue to leave data out the fit line starts
%incorportating the noise into the signal reducing it's accuracy when
%testing the point left out of the model. 


%% PART B

%%
%* Use bootstrapping methods to estimate the 95% confidence interval of the average 
% paralimbic volume of the decompression treatment group vs. the control group. 

% TREATMENT GROUP
treatGroup = DATA(:,3) == 1;
treatGroup = DATA(treatGroup,:);

empiricalTreatMean = mean(treatGroup(:,1));

numSamples = 1e6; 
resampledTreatMeans = nan(numSamples,1);
n1 = length(treatGroup);

for ii = 1:numSamples
    index = randi(n1, [n1, 1]);
    resampledTreat = treatGroup(index,1);
    resampledTreatMeans(ii,1) = mean(resampledTreat);
end

sortedTreatMeans = sort(resampledTreatMeans); 

% pick cuttoffs
CIwidth = 95;
lowerBound = (100-CIwidth)/2;
upperBound = 100-lowerBound;

%find corresponding indices
lowerBoundIndex = round(lowerBound*length(sortedTreatMeans)/100);
upperBoundIndex = round(upperBound*length(sortedTreatMeans)/100);

% reach into the sorted mean differences to find the value of of the upper
% and lower CI
lowerTreatCI = sortedTreatMeans(lowerBoundIndex);
upperTreatCI = sortedTreatMeans(upperBoundIndex);

% CONTROL GROUP
controlGroup = DATA(:,3) == 0;
controlGroup = DATA(controlGroup,:);

empiricalControlMean = mean(controlGroup(:,1));

numSamples = 1e6; 
resampledControlMeans = nan(numSamples,1);
n1 = length(controlGroup);

for ii = 1:numSamples
    index = randi(n1, [n1, 1]);
    resampledControl = controlGroup(index,1);
    resampledControlMeans(ii,1) = mean(resampledControl);
end

sortedControlMeans = sort(resampledControlMeans); 

% pick cuttoffs
CIwidth = 95;
lowerBound = (100-CIwidth)/2;
upperBound = 100-lowerBound;

%find corresponding indices
lowerBoundIndex = round(lowerBound*length(sortedControlMeans)/100);
upperBoundIndex = round(upperBound*length(sortedControlMeans)/100);

% reach into the sorted mean differences to find the value of of the upper
% and lower CI
lowerControlCI = sortedControlMeans(lowerBoundIndex);
upperControlCI = sortedControlMeans(upperBoundIndex);

%%
%* If the random assignment worked, the confidence intervals should overlap. 
% Do they? Also, do these data suggest that there is a statistically reliable 
% difference to the general population in terms of paralimbic volume?

figure
histogram(resampledTreatMeans,100)
hold on
histogram(resampledControlMeans,100)
title('Mean Paralimbic Volume of Treament & Control Groups')
xlabel('Mean Paralimbic Volume')
ylabel('# of samples')
xline(empiricalTreatMean, 'r--','LineWidth',2);
xline(empiricalControlMean, 'r--', 'LineWidth',2,'HandleVisibility','off');
xline(lowerTreatCI, 'g', 'LineWidth',2); %treatment lower CI
xline(upperTreatCI, 'g', 'LineWidth',2,'HandleVisibility','off'); %treatment Upper CI
xline(lowerControlCI, 'b', 'LineWidth',2); %control lower CI
xline(upperControlCI, 'b', 'LineWidth',2,'HandleVisibility','off'); %control Upper CI
box off
set(gca, 'TickDir', 'out')
legend('Treatment Group', 'ControlGroup', 'Empirical Means', 'Treatment 95% CI', 'Control 95% CI', 'Location', 'bestoutside')


% The confidence intervals do overlap. Since the population median (0) is
% outside both groups confidence interval we can reliably say that the prison
% population has a statistically smaller paralimbic volume than the general
% population. 

%% PART C

%%
%* Do a suitable t-test to compare the mean PCL-R score of prisoners who did and did
% not undergo decompression therapy. What is the p-value? Assuming an alpha-level of
% 0.05, is this difference significant? Can you reject the null hypothesis that decompression
% therapy is ineffective in terms of decreasing PCL-R scores?

[H,P,CI,STATS] = ttest2(controlGroup(:,2), treatGroup(:,2), 'tail', 'right')

%Calculating tstat by hand:
controlMean = mean(controlGroup(:,2));
treatMean = mean(treatGroup(:,2));
treatSTD = std(treatGroup(:,2));
controlSTD = std(controlGroup(:,2));
N = length(treatGroup);
tval = (controlMean-treatMean) / sqrt((controlSTD^2+treatSTD^2)/N)

% The two group are statistically different because 0 does not land in the
% confidence interval of their mean difference. The p-value is .0442 so it
% just makes significance at the .05 level. 

% The Null hypothesis can be rejected given this current data set, however 
% since the population is not normally distributed the t-test assumptions
% will not be met in this case meaning that this t-value may not be
% accuratly supported. This is why we are doing permutation testing below. 

%% PART D

%%
%* Do a permutation test to assess whether decompression therapy has an effect. Designate
% an appropriate test statistic and calculate its exact p value.

%Comparing difference between the ratio of control group mean/total mean
%and the treatment mean/total mean. 

empiricalM = mean(controlGroup(:,2))/mean(DATA(:,2)) - mean(treatGroup(:,2))/mean(DATA(:,2));

%Ratio difference of .1508 between control and treatment. 

%Testing the new statistic 
combinedScores = [controlGroup(:,2); treatGroup(:,2)]; %stack the ratings vertically
n1 = length(controlGroup(:,2)); %1 sample size
n2 = length(combinedScores); %Combined sample size


N = 1e6;
nullDistributionM = nan(N,1); 

for ii = 1:N
    shuffleIndices = randperm(n2);
    randomizedRatings1 = combinedScores(shuffleIndices(1:n1)); %first half
    randomizedRatings2 = combinedScores(shuffleIndices(n1+1:n2));
    
    nullDistributionM(ii,1) = (mean(randomizedRatings1)/mean(DATA(:,2)))-(mean(randomizedRatings2)/mean(DATA(:,2)));
end


nullDistributionM = sort(nullDistributionM);

% pick cuttoffs (1x tailed)
CIwidth = 95;
lowerBound = 100-CIwidth;
upperBound = 100-lowerBound;

upperBoundIndex = round(upperBound*length(nullDistributionM)/100);

% reach into the sorted mean differences to find the value of of the upper
% and lower CI
upperNullCI = nullDistributionM(upperBoundIndex);

figure
histogram(nullDistributionM,100)
title('Null distribution of Test Statistic M')
xline(empiricalM, 'r--','LineWidth',2);
xline(upperNullCI, 'g','LineWidth',2);
xlabel('Resampled M, assuming change')
ylabel('Number of Samples')
legend('Resampled M', 'Empirical M', '95% CI')
set(gca, 'TickDir', 'out')
box off


exactP = sum(nullDistributionM>empiricalM)/length(nullDistributionM)

% Given our new test statistic the treatment group is significantly lower
% in their PCL-R scores than the control group so the decompression treatment 
% was effective when using a one-tailed test (because directionality is important)
% at the .05 level. Since we only cared about the possibility of the scores
% decreasing a two tailed test was not necessary, since the PCL-R scores
% are an exsiting metric where lower scores correspond to less psychopathy
% the directionality is valid. 
