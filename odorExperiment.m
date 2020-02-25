function [simResponses] = odorExperiment(numTrials1,numTrials2)
%simResponses = odorExperiment(numTrials1,numTrials2) where numTrials1 and 
% numTrials2 are the number of trials in a simulated experiment for condition 
% 1 and 2, respectively. simResponses is a (N × 3) matrix containing simulated 
% responses of each of your 3 voxels during N = numTrials1 + numTrials2 trials. 
% [Hint: use ndRandn from the previous problem]. 

D1mean = [2.9771 4.2488 4.9744];
D2mean = [9.9819 12.0838 11.0741];

D1Cov =  [26.0444 9.7876 4.8923; 9.7876 15.6932 3.3961; 4.8923 3.3961 3.3716;];
D2Cov = [12.8625 8.0714 2.9013; 8.0714 26.8591 6.2990; 2.9013 6.2990 4.4109;];

simData1 =  ndRandn(D1mean, D1Cov, numTrials1);
simData2 = ndRandn(D2mean, D2Cov, numTrials2);
simResponses = [simData1; simData2];

end

