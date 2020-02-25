%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 1 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 3

% Date intialized: Sept. 12th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc
%%  3. Geometry of linear transformations

%%
%* (a) Write a function 'plotVec2' that takes as an argument a matrix of height 
% 2, and plots each column vector from this matrix on 2-dimensional axes. It 
% should check that the matrix argument has height two, signaling an error if 
% not. Vectors should be plotted as a line from the origin to the vector 
% position, using circle or other symbol to denote the 'head' (see help for 
% function 'plot'). It should also draw the x and y axes, extending from -1 
% to 1. The two axes should be equal size, so that horizontal units are equal 
% to vertical units (read the help for the function 'axis').

testMat = randn(2,2);
plotVec2(testMat)

%%
%* (b) Write a second function vecLenAngle that takes two vectors as arguments
% and returns the length (magnitude, or Euclidean-norm, not dimensionality) 
% of each vector, as well as the angle between them. Decide how you would like 
% to handle cases when one (or both) vectors have zero length.
testVec1 = randn(2,1);
testVec2 = randn(2,1);

[len1, len2, angle] = vecLenAngle(testVec1, testVec2)

%%
%* (c) Generate a random 2x2 matrix M, and decompose it using the SVD, M = USV^T. 
% Now examine the action of this sequence of transformations on the two 
% 'standard basis' vectors, {e1, e2}. Specifically, use 'vecLenAngle' to 
% examine the lengths and angle between two basis vectors e'n, the two vectors 
% V^Te'n, the two vectors SV^Te'n, and the two vectors USV^Te'n. Do these 
% values change, and if so, after which transformation? Verify this is 
% consistent with their visual appearance by plotting each pair using plotVec2.

matrixM = randn(2,2);

[U, S, V] = svd(matrixM);

test = U*S*V'; %verifies that matrixM = U*S*V' (idiot check)

% e1 and e2 standard basis vectors
e = ([1 0; 0 1;]);
[Elen1, Elen2, Eangle] = vecLenAngle(e(:,1), e(:,2));
plotVec2(e)

%1st transformation V - rotates the vectors of e
EV=V'*e;
[Vlen1, Vlen2, Vangle] = vecLenAngle(EV(:,1), EV(:,2));
plotVec2(EV)

%2nd transformation SV - this is the stage in which the length and angle
%values of the vectors are subject to change. 
SVE = S*EV;
[SVlen1, SVlen2, SVangle] = vecLenAngle(SVE(:,1),SVE(:,2));
plotVec2(SVE)

%3rd transformation USV - no change occurs in magnitude or angle, only
%flips or rotation are present 
USV = U*SVE;
[USVlen1, USVlen2, USVangle] = vecLenAngle(USV(:,1),USV(:,2));
plotVec2(USV)

disp(['The first transformation, V^T rotates or flips the basis vectors,']) 
   disp(['the second transformation S stretches or shrinks the vectors, and the'])
   disp(['last transformation U flips or rotates the stretched vectors a second time'])
%%
%* (d) Generate a data matrix P with 65 columns containing 2-dimensional 
% unit-vectors u'n = [cos(?n); sin(?n)], and ?n = 2?n/64, n = 0, 1, . . . , 64. 
% [Hint: Don't use a for loop! Create a vector containing the values of ?n. ] 
%     Plot a single blue curve through these points, and a red star (asterisk)
%     at the location of the first point. 

dataMat = zeros(2,65); %output matrix

counting = [0:64]; %counting variable to increase incrementally 
thetaMath = 2*pi/64; %calculate theta
thetaVec = thetaMath'*counting %multiply by all instances of n

dataMat(1,:) = cos(thetaVec); %1st row is all cos of theta values for each n
dataMat(2,:) = sin(thetaVec); %2nd row is all sin of theta values for each n

fig1 = figure;
x = dataMat(1,:);
y = dataMat(2,:);
plot(x,y)
hold on
plot(dataMat(1,1), dataMat(2,1), '*')
lineHandles = get(gca, 'children');
set(lineHandles, 'LineWidth', 2)
xlabel('X')
ylabel('Y')
axis equal
box off
title('Cos of Theta plotted by Sin of Theta')
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 14)

%% 
%* Consider the action of the matrix M from the previous problem on this 
% set of points. In particular, apply the SVD transformations one at a time 
% to full set of points (again, think of a way to do this without using a 
% for loop!), plot them, and describe what geometric changes you see (and why).

%Original plot
fig2 = figure;
x = dataMat(1,:);
y = dataMat(2,:);
plot(x,y, 'b')
hold on
plot(dataMat(1,1), dataMat(2,1), '*','HandleVisibility','off')
hold on
lineHandles = get(gca, 'children');
set(lineHandles, 'LineWidth', 2)
xlabel('X')
ylabel('Y')
axis equal
box off
title('Transformations of a Circle by SVD of MatrixM')
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 14)

%1st transformation V - original plot has been rotated based on the values
%in V'
dataV=V'*dataMat;
[Vlen1, Vlen2, Vangle] = vecLenAngle(dataV(:,1), dataV(:,2));
x = dataV(1,:);
y = dataV(2,:);
plot(x,y, 'g')
lineHandles2 = get(gca, 'children');
set(lineHandles2, 'LineWidth', 2)
hold on
plot(dataV(1,1), dataV(2,1), '*','HandleVisibility','off')
hold on


%2nd transformation S - This stretches or shrinks the circle based on the
%values in S- the circle now becomes an oval. 
SVd = S*dataV;
[SVlen1, SVlen2, SVangle] = vecLenAngle(SVd(:,1),SVd(:,2));
x = SVd(1,:);
y = SVd(2,:);
plot(x,y, 'y')
lineHandles3 = get(gca, 'children');
set(lineHandles3, 'LineWidth', 2)
hold on
plot(SVd(1,1), SVd(2,1), '*','HandleVisibility','off')
hold on


%3rd transformation U - no change occurs in magnitude or angle, only
%flips or rotation of the oval are present 
USV = U*SVd;
[USVlen1, USVlen2, USVangle] = vecLenAngle(USV(:,1),USV(:,2));
x = USV(1,:);
y = USV(2,:);
plot(x,y, 'm')
lineHandles4 = get(gca, 'children');
set(lineHandles4, 'LineWidth', 2)
hold on
plot(USV(1,1), USV(2,1), '*','HandleVisibility','off')
legend('Original', 'Rotated by V^T', 'Stretched by S', 'Rotated by U', 'Location', 'NorthEastOutside')


disp('The circle is rotated after the transformation by V^T but remains the')
 disp('same size and shape. Transforming by S makes the circle into an ellipse')
   disp('and then the U transform rotates that ellipse.')



