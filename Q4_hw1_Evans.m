%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 1 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 4

% Date intialized: Sept. 12th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc
%% A simple visual neuron. 

% Suppose a retinal neuron in a particular species of toad generates responses 
% that are a weighted sum of the (positive-valued) intensities of light that 
% is sensed at 6 localized regions of the retina. The weight vector is 
% [1, 3, 8, 8, 3, 1]. 

weightVec = [1, 3, 8, 8, 3, 1];

%%
%* (a) Is this system linear? If so, how do you know? 
% If not, provide a counterexample. 

disp('This is a linear system because it supports additivity and homogeneity')
disp('this is displayed by breaking a random input down to its parts then running')
disp('each part separatly though the system. The sum of these individual outcomes')
disp('is equal to the result of the complete vector going though the system.')

testVec = abs(randi(100,[1,6])); %random vector to test the system

IDmat = eye(6); %identiy matrix

for ii = 1:6
out(ii,:) = (testVec.*IDmat(ii,:)).*weightVec; %isolate each basis of the test vector and apply the system
end

testOutcome = testVec.*weightVec; %test vector applyed to the system
linearTest = sum(out) == testOutcome %test showing sum of basis vectors is the same as test vector after the system is applied
disp('The logical output shows that the sum of the weighted elements applied to')
disp('the system separatly is equal to the test vector being applied to the')
disp('system directly.') 
%%
%* (b) What unit-length stimulus vector 
% (i.e., vector of light intensities) elicits the largest response in the 
% neuron? Explain how you arrived at your answer. 

%test in 2D
testWeight = [3; 8;]; %2D test weight vector
weightlen = sqrt(testWeight'*testWeight);
unitWeight = testWeight/weightlen; %unit vector in line with the weight vector

%creating 10 random positive unit vectors 
for ii = 1:10
    
randUnit = abs(randn(2,1));
lenUnit = sqrt(randUnit'*randUnit);
unitVec = randUnit/lenUnit;

unitVecMat(:,ii) = unitVec;
end

%adding unit vectors at the maximum edges of postive values
unitVecMat(:,11) = [0, 1]; 
unitVecMat(:,12) = [1,0];
unitVecMat(:,13) = unitWeight; %adding the unit vector in line with the weight vector

dotted = unitVecMat'*testWeight; %apply the weights to all vectors

maxDot = find(dotted == max(dotted)); %find the largest response
maxUnitVec = unitVecMat(:,maxDot); %Locate the unit vector that ilicits the response
check = maxUnitVec == unitWeight;
disp('The largest value would come from a unit vector laying in the same direction as the weight vector')

minDot = find(dotted == min(dotted)); %find the smallest response
minUnitVec = unitVecMat(:,minDot); %find the unit vector that produces the smallest response
disp('The minimum value would come from a unit vector laying on the axis furthest from the weight vector, in this test case (1,0)')
 
figure
plotv(testWeight)
lineHandles = get(gca, 'children');
set(lineHandles, 'LineWidth', 2)
hold on
plotv(unitVecMat)
axis([0 1 0 1])
title('Unit Vectors and Test Weight Comparison')
xlabel('X')
ylabel('Y')
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 14)
box off

%%
%*Answer with regards to the 6D weight vector

weightTest = [1; 3; 8; 8; 3; 1;];
weightlen = sqrt(weightTest'*weightTest);
unitWeight = weightTest/weightlen;

disp('The "unitWeight" unit vector elicit the greatest response from the neuron')
disp('as it is the unit vector pointing in the direction of the weight vector')
disp('and thus has the greatest dot product.')



%%
%* 
disp('I first solved the above question this way using a testing and checking')
disp('method that more strongly depends on programming power than theory-')
disp('however I found it helpful to really show with so many possible options that the above') 
disp('reached the same ball park number (hard to be perfectly exact when randomizing). Thus Ive left it on')
disp('the script for interest sake but the formal answer is above.')

unitMat = zeros(6,10e6); %output matrix

for ii = 1:10e6
    
randUnit = abs(randn(6,1)); %random positive 6D vector
lenUnit = sqrt(randUnit'*randUnit); %length of vector
unitVec = randUnit/lenUnit; %rescaled to be a unit vector 

unitMat(:,ii) = unitVec; %added to output matrix
end

sortrows(unitMat);

weightTest = [1; 3; 8; 8; 3; 1;];

scaledUnits = unitMat'*weightTest; %scaled all unit vectors by system weights

maxOutput = max(scaledUnits) %find the largest output
maxRow = find(scaledUnits == maxOutput); %find the row of the unit vector that matches the maximum output

maxUnitVec = unitMat(:,maxRow) %unit vector that gives maximum output based on the randomized set
unitWeight %from the offical solution above

disp('Created 10e6 different positive unit vectors and applied the weights to')
disp('them all indiviually and found which of the unit vectors provided the largest output.')
disp('That unit vector is the option that is closest to the weight vector out of the randomized inputs')


%%
%* (c) What physically-realizable unit-length stimulus vector produces the 
% smallest response in this neuron? Explain. 
%[hint: think about the geometry by visualizing a simpler version of the problem, in 2 dimensions]

smallVec1 = [0; 0; 0; 0; 0; 1;];
smallVec2 = [1; 0; 0; 0; 0; 0;];

minOutput1 = smallVec1'*weightTest;
minOutput2 = smallVec2'*weightTest;

disp('The unit vector with the smallest possible output will be either one of')
disp('the two above as they are the vectors closest to the axes which are')
disp('further from the weight vector.')