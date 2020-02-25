function [gramMat] = gramSchmidt(N)
% This function takes a single argument, N, specifying 
% the dimensionality of the basis. It should then generate an N × N matrix 
% whose columns contain a set of orthogonal normalized unit vectors.
% First, one generates an arbitrary unit vector 
% (e.g., by normalizing a vector created with randn). Each subsequent basis 
% vector is created by generating another arbitrary vector, subtracting off 
% the projections of that vector along each of the previously created basis 
% vectors, and normalizing the remaining vector.

vecMat = randn(N);

unit1 = vecMat(:,1)/sqrt(vecMat(:,1)'*vecMat(:,1));
e1 = unit1/sqrt(unit1'*unit1);
gramMat(:,1) = e1;

 for ii = 2:N
        vec = vecMat(:,ii);
        for jj = 1: ii-1
            vec = vec - proj(vecMat(:,ii),gramMat(:,jj));
       
        end
        gramMat(:,ii) = vec / sqrt(vec'*vec);
 end
           
orthoMat = gramMat'*gramMat; %check to verify the resulting matrix is orthogonal 



end 





