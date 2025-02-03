function normal = Normal3points(central,left,right)
dir = left  - right;

% candidates (all possible normal vectors)
n_vec1 = [-dir(2),  dir(1)];  
n_vec2 = [ dir(2), -dir(1)]; 


to_center = central-right;

% Test which normal points away from the center
if dot(n_vec1, to_center) < 0
    normal = n_vec1; % n1 points to the opposite side
else
    normal = n_vec2; % n2 points to the opposite side
end

normal = normal' / norm(normal); % normalization