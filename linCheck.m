function [sanityCheck1, sanityCheck2, sanityCheck3] = linCheckSum(numCheck)
%This function checks the linearity of unknown systems 64 rows in length.
%Inputs are first which system you'd like to test, followed by how many
%random points at the system you'd like to compare (10 points leads to 9
%comparisons). 

sampleVec = randi(64,1,numCheck);

for ii = 1:numCheck-1
    
    testSoundA = zeros(64,1);
    testSoundB = zeros(64,1);
    testSoundC = zeros(64,1);
    sample1 = sampleVec(ii);
    sample2 = sampleVec(ii+1);
    
    testSoundA(sample1) = 1;
    testSoundA(sample2) = 1;
    
    testSoundB(sample1) = 1;
    testSoundC(sample2) = 1;
    
    outputA = unknownSystem1(testSoundA);
    outputB = unknownSystem1(testSoundB);
    outputC = unknownSystem1(testSoundC);
    
    test1 = sum(outputA);
    test2 = sum(outputB + outputC);
    
    sumCheck1A = round(test1, 6);
    sumCheck1B = round(test2, 6);
    
    sanityCheck1(ii,1) = sumCheck1A == sumCheck1B;
    
end
for ii = 1:numCheck-1
    clear test1 test2 outputA outputB outputC;
    testSoundA = zeros(64,1);
    testSoundB = zeros(64,1);
    testSoundC = zeros(64,1);
    sample1 = sampleVec(ii);
    sample2 = sampleVec(ii+1);
    
    testSoundA(sample1) = 1;
    testSoundA(sample2) = 1;
    
    testSoundB(sample1) = 1;
    testSoundC(sample2) = 1;
    
    outputA = unknownSystem2(testSoundA);
    outputB = unknownSystem2(testSoundB);
    outputC = unknownSystem2(testSoundC);
    
    test1 = sum(outputA);
    test2 = sum(outputB + outputC);
    
    sumCheck2A = round(test1, 6);
    sumCheck2B = round(test2, 6);
    
    sanityCheck2(ii,1) = sumCheck2A == sumCheck2B;
    
end
for ii = 1:numCheck-1
    clear test1 test2 outputA outputB outputC;
    testSoundA = zeros(64,1);
    testSoundB = zeros(64,1);
    testSoundC = zeros(64,1);
    sample1 = sampleVec(ii);
    sample2 = sampleVec(ii+1);
    
    testSoundA(sample1) = 1;
    testSoundA(sample2) = 1;
    
    testSoundB(sample1) = 1;
    testSoundC(sample2) = 1;
    
    outputA = unknownSystem3(testSoundA);
    outputB = unknownSystem3(testSoundB);
    outputC = unknownSystem3(testSoundC);
    
    test1 = sum(outputA);
    test2 = sum(outputB+outputC);
    
    sumCheck3A = round(test1, 6);
    sumCheck3B = round(test2, 6);
    
    sanityCheck3(ii,1) = sumCheck3A == sumCheck3B;
    
end
end

