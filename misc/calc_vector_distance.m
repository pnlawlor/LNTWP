function [ angle ] = calc_vector_distance( vec1, vec2)

% Add a very small number to avoid imaginary angles
vec2 = vec2 + 10*eps;

cos_angle = dot(vec1,vec2)/(norm(vec1)*norm(vec2));
angle = acos(cos_angle);

end

