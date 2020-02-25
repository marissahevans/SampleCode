function [lenVec1,lenVec2, angleVecs] = vecLenAngle(vec1,vec2)
%This function takes 2x vectors and returns their lengths and the angle between them


vecDot = vec1'*vec2;

lenVec1 = sqrt(vec1'*vec1);
lenVec2 = sqrt(vec2'*vec2);

if lenVec1 && lenVec2 == 0
    error('Neither of the vectors is greater than zero, no angle can be calculated')
    
elseif lenVec2 == 0
    error('One of the vectors has a length of zero, cannot compute angle')
    
elseif lenVec1 == 0
    error('One of the vectors has a length of zero, cannot compute angle')

else
angleVecs = acosd(vecDot/(lenVec1*lenVec2));

end 
end

