%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MATHTOOLS HOMEWORK 3 - MARISSA EVANS - mhe229@nyu.edu %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION 2

% Date intialized: Oct. 15th 2019
% Written by: Marissa Evans

%Start with a fresh work space
clear all 
close all
clc

%% Convolution in matlab. 

% Suppose you have a linear shift-invariant system with impulse response:
r = [4 2 1 0.5];
% Because it is LSI, the response of this system to any input vector in can be computed as a
% convolution.

%% PART A

%%
%*Compute responses to eight 8-dimensional impulse vectors, using MATLAB's conv
% function: out = conv(in, r). How do these compare to what you'd expect from
% the convolution formula given in class, [SEE ORIGINAL PROMPT FOR FORMULA] 
% Specifically, compute the matrix that represents the linear system. 
% What is the size, and organization of this matrix?

impulseVecs = eye(8);

for ii = 1:8
out(ii,:) = conv(impulseVecs(ii,:), r);
end 

disp(out)

% The formula from class gives an output with zero padding at the beginning
% as the first combindation of the first impulse vector and the values of r
% flipped will give you what is the 4th column of the matlab conv output.
% The matrix you end up with will be 8x8 due to the stretching from the
% impulse vector, however it will not be as large as the 8x11 matlab matrix
% due to the additional padding used there. The conv function flips the
% kernel while the formula flips the input. 


% The linear system matrix is larger than the inputs used. The output is expanded
% to 8x11 instead of 8x8 as the values are let to run to completion on each
% row, with padding used on the shorter rows. The values do not wrap around
% so the signal in this case is not periodic. The values of r are represented
% in the columns starting with the last values of the vector and they are
% also all seen across the rows in the same orientation as the original
% impulse repsonse. 

%% PART B

%%
%* How does MATLAB's conv function handle boundaries?

% When in 'full' the conv function adds zeros onto the ends of all the
% vectors, creating the largest matrix possible to include all values from
% the orginal signal. this is the automatic setting for conv. 

% When 'valid' is selected the output will only include the center colums
% which contain a full range of values, clipping both the edges where zeros
% would be added, this makes the output smaller than the input but does not
% change the shape of the data. 

% When 'same' is selected the edges are slighly shortened to make the row
% length the same as the input vector. Some padding with zeros is included
% however not as much as the full option. 

%% PART C

%%
%* Using conv, compute the response to an input vector of length 48 
% containing a singlecycle cosine. Is this a single-cycle sinusoid? Why or why not? 

k = 1;
n = 48;
x = 0:n-1;

input = cos(2*pi*k/n * x);

output = conv(input, r);

figure
plot(input, '--k', 'Linewidth',2) 
hold on 
plot(output, 'r','Linewidth',2)
xlabel('Sample #')
ylabel('Amplitude')
title('Input Vector Convolved with r (full)')
legend ('Input Wave', 'Output Function')
box off
set(gca, 'TickDir', 'out')

% After the convolution the output is no longer a single cycle sinousoid.
% This is because of the padding of zeros on the edges of the function
% because they falsely decrease the values toward the ends of the vector. 

%%
%* If not, what modification is necessary to the conv function to ensure that it 
% will behave according to the 'sinein, sine-out' behavior we expect of LSI systems?

output = conv(input, r, 'valid');

figure
plot(input, '--k','Linewidth',2) 
hold on 
plot(output, 'r','Linewidth',2)
legend ('Input Wave', 'Output Function')
xlabel('Sample #')
ylabel('Amplitude')
title('Input Vector Convolved with r (valid)')
box off
set(gca, 'TickDir', 'out')

%By using the 'valid' function this does not add any extra zero padding so
%it holds the shape constant from the input to output. This does however
%clip a small amount of the wave. To insure no data will be lost, adhoc
%choices would have to be made as to how to design the boundaries, whether
%they are periodic or some other form of data repitition based on the
%structure and magnitude of the data. For instance the padding could be the
%original values mirrored to keep the continuing wave form, or an
%extrapolation of the edge values to keep the trend going in it's current
%direction. 