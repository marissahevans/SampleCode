function Y = randp(p, num, plot)

% 3rd argument: Enter 1 to get a plot, 0 for no plot

% Uses a clipped Poisson vector p as input. The function generates num samples from the probability 
% distribution function (PDF) specified by p. [Hint: use the rand 
% function, which generates real values over the interval [0...1], and 
% partition this interval into portions proportional in size to the probabilities in p]. 
randNums = rand(num,1);
k = 0:length(p)-1;

bin(1) = 0;
for ii = 1:length(p)-1
    bin(ii+1) = bin(ii) + p(ii);
end

Y = discretize(randNums,bin);
Y = Y-1;

if plot ==1 
figure
histogram(Y,k)
title(['Poisson Distribution with ' num2str(num) ' Samples'])
xlabel('Spikes/Interval')
ylabel('# of Samples')
box off
set(gca, 'TickDir', 'out')

else
end

end