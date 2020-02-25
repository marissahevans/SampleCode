function [] = plotVec2(plotmat)
% This function takes as an argument a matrix of height 2,
% and plots each column vector from this matrix on 2-dimensional axes. It 
% checks that the matrix argument has height two, signaling an error if 
% not. Vectors are plotted as a line from the origin to the vector 
% position, using circle or other symbol to denote the 'head' 
% It also draws the x and y axes, extending from -1 
% to 1. The two axes are equal size, so that horizontal units are equal 
% to vertical units. 
sizeMat = size(plotmat);
figure;
title('2D plotVec2 Plot')
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 14)
yline(0,'r:');
xline(0,'r:');
hold on;
axis equal
xlabel('X')
ylabel('Y')
axis([-1 1 -1 1])


if sizeMat(1) > 2
    error('Matrix may not exceed 2 rows')
else
    for ii = 1:sizeMat(2)
        plotv(plotmat(:,ii));
        hold on
        plotv(plotmat(:,ii),'kd');
        hold on
    end
end 
lineHandles = get(gca, 'children');
set(lineHandles, 'LineWidth', 2);
end

