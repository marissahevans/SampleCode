%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 1 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 2

% Date intialized: Sept. 12th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc
%% 2.  Inner product with a unit vector. Given unit vector u, and an arbitrary vector v, write
% (MATLAB) expressions for computing:

% (a) the component of v lying along the direction u?,
unitVec = []; %arbitrary unit vector u
vecV = []; %arbitrary vector v

dotV = vecV'*unitVec; 
compVec = dotV*unitVec; %component of v along with u

% (b) the component of v that is orthogonal (perpendicular) to u?
orthoVec = vecV - compVec;

% (c) the distance from v to the component that lies along direction u?.
lengthOrtho = sqrt(orthoVec'*orthoVec);

%% NOW WITH NUMBERS
% Now convince yourself your code is working by testing it on random vectors u and v 
% (generate these using randn, and don't forget to re-scale u so that it has unit length). 

vecV = randn(2,1); %random vector v

randUnit = randn(2,1);
lenUnit = sqrt(randUnit'*randUnit);
unitVec = randUnit/lenUnit; %random unit vector u

% (a) w/ numbers generated above
dotV = vecV'*unitVec; 
compVec = dotV*unitVec; %component vector to v

% (b) w/ numbers generated above
orthoVec = vecV - compVec; %orthoginal vector to component vector

%%
%* First, do this visually with 2-dimensional vectors, by plotting u, v, and the two components described
% in (a) and (b). (hint: execute 'axis equal' to ensure that the horizontal and vertical axes
% have the same units). 

fig1 = figure;
plotv(compVec)
hold on
plotv(orthoVec)
hold on
plotv(unitVec)
hold on
plotv(vecV)
legend('compVec', 'orthoVec', 'unitVec', 'vecV')
lineHandles = get(gca, 'children');
set(lineHandles, 'LineWidth', 2)
axis equal
title('2D Vectors a, b, u, & v')
xlabel('X')
ylabel('Y')
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 14)


%% NOW IN 4 DIMENSIONS
% Then test it numerically in higher dimensions (e.g., 4) by writing (and
% running) expressions to verify each of the following:

vecV = randn(4,1); %random vector v in 4 dimensions

randUnit = randn(4,1);
lenUnit = sqrt(randUnit'*randUnit);
unitVec = randUnit/lenUnit; %random unit vector u in 4 dimensions

dotV = vecV'*unitVec; 
compVec = dotV*unitVec; %component vector in 4 dimensions
orthoVec = vecV - compVec; %orthoginal vector in 4 dimensions

%%
%* the vector in (a) lies along the same line (i.e., points in the same or opposite direction)
% as u?. 

dotUA = compVec'*unitVec;
lengthA = sqrt(compVec'*compVec);
lengthU = sqrt(unitVec'*unitVec);


angleUA = acosd(dotUA/(lengthA*lengthU)); %angleUA will either equal 180º or 0º 
disp(['The vector in (a) lies along the same line as (u) because the angle between them is ' num2str(angleUA) ' degrees'])

%%
%* the vector in (a) is orthogonal to the vector in (b).
dotAB = compVec'*orthoVec;
lengthB = sqrt(orthoVec'*orthoVec);


angleBA = acosd(dotAB/(lengthA*lengthB)); %angleBA will equal 90º
disp(['The vector in (a) is orthogonal to the vector in (b) because the angle between them is ' num2str(angleBA) ' degrees'])


%%
%* the sum of the vectors in (a) and (b) is equal to v.

sumAB = compVec + orthoVec;

disp(['The values in vector sumAB (' num2str(sumAB') ') are equal to the values in vector vecV (' num2str(vecV') ')'])

%%
%* the sum of squared lengths of the vectors in (a) and (b) is equal to ||v||2.

sqrLenV = vecV'*vecV;

sqrLenAB = (compVec'*compVec) + (orthoVec'*orthoVec);

disp(['The length of ||v||2 is ' num2str(sqrLenV) ' which is equal to the sum of squared lengths of vectors (a) and (b) which is ' num2str(sqrLenAB)])



