function [mseXval] = betaMathTools2(x, y, numPoly)
%This function calculates the beta value and error value for a given x and
%y data set. It also plots the polynomial fit to the data. 
figure
plot(x,y, 'o')
axis equal
axis square
box off
xline(0,'--','HandleVisibility','off');
yline(0, '--','HandleVisibility','off');
hold on

for jj = 1:length(x)
    x2 = x;
    y2 = y;
    
    testValx = x2(jj);
    testValy = y2(jj);
    
    x2(jj) = nan;
    y2(jj) = nan;
    x2 = x2(~isnan(x2));
    y2 = y2(~isnan(y2));

m = length(x2);
pp = 1:numPoly;

XX = ones(m,numPoly);

for ii = 0:length(pp)-1
XX(:,numPoly-ii) = XX(:,numPoly-ii).*(x2.^((numPoly-1)-ii));

end

[U, ~, V] = svd(XX);
sVec = svd(XX);
sInv = zeros(size(XX));

for ii = 1:numPoly
sInv(ii,ii) = 1./sVec(ii);
end

betaVal = V*sInv'*U'*y2;

errorVal = XX*betaVal;

fit_test = testValx * betaVal;

mseXval(jj,1) = mean((testValy - fit_test).^2); 
                

plot(x2,XX*betaVal, 'LineWidth', 2)
title(['Plot with Polynomials at Order ' num2str(numPoly-1)])
legend('Data', 'PolyFit Line')
xlabel('X');
ylabel('Y');
set(gca, 'TickDir', 'out')
box off
set(gca, 'FontSize', 12)
hold on
end