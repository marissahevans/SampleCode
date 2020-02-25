%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 2 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 3

% Date intialized: Sept. 27th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Constrained Least Squares Optimization. 
% Load the file constrainedLS.mat into 
% MATLAB. This contains an N × 2 data matrix, data, whose columns correspond to horizontal
% and vertical coordinates of a set of 2D data points, ~dn. It also contains a 2-vector w. Consider
% a constrained optimization problem (SEE PRINT OUT):
% Thus, the constraint on ~v is that it must lie on a line, perpendicular to ~w, whose perpendicular
% distance from the origin is 1/|| ~w||.

load('constrainedLS.mat');


%% PART 1

%%
%*(a) Rewrite the optimization problem in matrix form. 
% Then rewrite the problem in terms
% of a new optimization variable, v? (a linear transformation of ~v), such that the quantity
% to be minimized is now ||v?||2. Note: you must also rewrite the constraint in terms of v?.

disp('Problem rewritten in Matrix form')
disp('min(v)||Dv||^2 -> SVD(D) -> ||USV"v||^2 -> ||Sv*|| ->||v~||^2')
disp('Moving v into the transformed space (~)')

disp('v~ = SV"v ')
disp('v = VS#v~ - showing the components that make up v using v~')

disp('Constraint rewritten in terms of v~')
disp('(VS#v~)"*w = 1')

%%
%* (b) The transformed problem is one that you should be able to solve. In particular, you
% must find the shortest vector v that lies on the constraint line. Compute the solution for
% v~, and plot the solution vector, the constraint line and the transformed data points.

D = data;

%Run the SVD on the data set
[U S V] = svd(D);

%Calculate the inverse of S
Sinv(1,1) = 1/S(1,1);
Sinv(2,2) = 1/S(2,2);

%Find the square matrix for S
sVal = find(diag(S));
S2 = S(sVal,:);

%Transform w into the v~ space
wSS = Sinv'*V'*w; % w~

%Find 1/||w|| in the v~ space
lenWSS = sqrt(wSS'*wSS); %||w~||
distOrigSS = 1/lenWSS; %1/||w~||

%Find vector points with length 1/||w~|| to origin
x = (wSS(1,1)*distOrigSS)/lenWSS; %find x value using similar right triangles law
y = sqrt(distOrigSS^2 - x^2); %find y using pythagorean theorem 
vPoint = [x; -y]; %y is negative due to axis points
lenVTil = sqrt(vPoint'*vPoint); %length of vector v~
test = distOrigSS == lenVTil; %check that  ||v~|| is the same as 1/||w~||

% Plotting the constraint and w in the v~ space
slope = vPoint(1)/vPoint(2); %find slop of v~ (same as w~)
perpSlope = -1/slope; %find the inverse slope
x2 = -50:50;
x1 = vPoint(1); 
y1 = vPoint(2); 
y2 = perpSlope*(x2 - y1) + x1;

%Transform the data into the v~ space
DSS = D*V*Sinv;
disp('Transformed data is circular (although very small when plotted with the constraint), it"s dot product is the idenity')
IDMAT = DSS'*DSS %shows that the dot product is the identity thus data is circular


transFig = figure;
plot(DSS(:,1), DSS(:,2), '.')
hold on
plot(y2, x2, 'LineWidth', 2)
plot([0,vPoint(1,1)],[0, vPoint(2,1)], 'LineWidth', 2)
xline(0);
yline(0);
axis equal
axis([-20 20 -20 20])
title('Constrain line, vecV and Data in v~ space')
legend('data v~','constrain v~','vec v~','Location','bestoutside')
box off
set(gca, 'TickDir', 'out')
xlabel('X')
ylabel('Y')


%%
%* (c) Transform the solution back into the original space (i.e., solve for ~v). Plot ~v, the original
% constraint line, and the original data points. Is the optimal vector ~v perpendicular to
% the constraint line? 

%rotate back to the origninal space
vecV = V*Sinv*vPoint;

%Constraint line in original space (using w)
lenW = sqrt(w'*w); %||w||
unitW = w/lenW; 
distOrig = 1/lenW; %1/||w||

xOrig = (w(1,1)*distOrig)/lenW;
yOrig = sqrt(distOrig^2 - xOrig^2);
vOPoint = [xOrig; yOrig]; 

lenVTil = sqrt(vOPoint'*vOPoint);
test2 = distOrig == lenVTil;
slope2 = vOPoint(1)/vOPoint(2);
perpSlope2 = -slope2;
x3 = -5:5;
x4 = vOPoint(1); 
y4 = vOPoint(2); 
y3 = perpSlope2*(x3 - x4) + y4;

% Plot the constrain line, v, w and data in the original space
origFig = figure;
plot(data(:,1),data(:,2),'.')
hold on
plot (x3,y3, 'LineWidth', 2)
plot([0,vecV(1,1)],[0, vecV(2,1)],'LineWidth',2)
plot([0,w(1)],[0,w(2)], 'LineWidth', 2)
xline(0);
yline(0);
axis equal
axis([-5 5 -5 5])
title('Constrain line, vecV, W and Data in Original space')
legend('data', 'constrain', 'vecV','w','Location','bestoutside')
box off
set(gca, 'TickDir', 'out')
xlabel('X')
ylabel('Y')

disp('The optimal vecV is not perpindicular to the constraint line in the')
disp('original space.')

%% PART 2
%%
%* On the same graph, plot the total least squares solution (i.e., the
% vector that minimizes the same objective function, but that is constrained to be a unit
% vector). Are the two solutions the same?


disp('LEAST SQUARES SOLUTION')
disp(' ||Dû||^2 s.t. ||û||^2 = 1 (is a unit vector)')
disp('||Dû||^2 = ||USV"û||^2 = ||SV"û||^2 = ||Sû*||^2 = ||u~||^2')

%Create a unit circle to represnt the unit length constraint options
unitCircle = [cos(linspace(0,2*pi));sin(linspace(0,2*pi))];

%transform circle to the v~ space
uSS = S2*V'*unitCircle;

%find the point closest to the origin
minYval = min(uSS(2,:));
minU = [0; minYval];

%rotate that point back to the original space
originalUnit = V*Sinv*minU;

% Plotting the unit circle in the original space, including the best choice
% unit vector based on the transformation (in addition to the original
% constraint line, v and w)
figure
plot(data(:,1),data(:,2),'.')
hold on
plot (x3,y3, 'LineWidth', 2)
plot([0,vecV(1,1)],[0, vecV(2,1)],'LineWidth',2)
plot(unitCircle(1,:), unitCircle(2,:), 'LineWidth', 1)
plot([0,originalUnit(1)], [0,originalUnit(2)],'LineWidth', 2)
plot(originalUnit(1),originalUnit(2), '*')
plot([0,w(1)],[0,w(2)], 'LineWidth', 2)
xline(0);
yline(0);
axis equal
axis([-5 5 -5 5])
title('Unit circle in original space')
legend('data', 'constrain', 'vecV','unit circle', 'min unit vector', 'min point','w','Location','bestoutside')
box off
set(gca, 'TickDir', 'out')
xlabel('X')
ylabel('Y')

%For comparison I'm also plotting the transformed unit circle in the v~ space
figure
plot(DSS(:,1), DSS(:,2), '.')
hold on
plot(y2, x2, 'LineWidth', 2)
plot([0,vPoint(1,1)],[0, vPoint(2,1)], 'LineWidth', 2)
plot(uSS(1,:), uSS(2,:), 'LineWidth', 1)
plot(minU(1),minU(2),'*')
xline(0);
yline(0);
axis equal
axis([-20 20 -20 20])
title('Unit circle transformed into v~ space')
legend('data v~','constrain v~','vec v~', 'unit circle in v~ space', 'min unit point','Location','bestoutside')
box off
set(gca, 'TickDir', 'out')
xlabel('X')
ylabel('Y')

disp('The two solutions are not the same because they are depending on different')
disp('ways to constrain the error, although they appear in similar space on the graph')
disp('the first is constrained by a line 1/||w|| away from the origin while the second')
disp('is constrained by a line 1 away, and these values are not the same especially')
disp('when the slope of the line is taken into account')
