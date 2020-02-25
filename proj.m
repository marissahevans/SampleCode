function [p] = proj(vec1,vec2)
%This function projects the vector vec1 orthogonally onto the line spanned by vector vec2. 
p = (vec1' * vec2) / (vec2' * vec2) .* vec2;

end

