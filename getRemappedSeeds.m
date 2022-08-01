function X = getRemappedSeeds(bx,X,per_x,per_y)

% Funtion to get copies of seeds in the fundamental domain, which is
% symmetric about the origin in the directions in which it is periodic.
    
    p = [per_x,per_y];
    
    bxDims = [bx(3)-bx(1),bx(4)-bx(2)]; % domain dimensions
    Binv   = diag(p./bxDims);           % 2x2 matrix to normalise seeds in directions in which domain is periodic
    
    % Get n x 2 matrix
    % k = argmin_{l}|Z + l*diag(bxDims)|,
    % where the minimum is taken over all n x 2 matrices with integer
    % entries and n is the number of seeds
    k = floor(-X*Binv+0.5*ones(1,2));
    
    X = X + k*diag(bxDims);

end