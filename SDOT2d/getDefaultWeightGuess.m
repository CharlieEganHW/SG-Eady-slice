function w = getDefaultWeightGuess(bx,Z,per_x,per_y)

    if (per_x == false) && (per_y == false)
        error("Default weight guess not valid when domain is not periodic in either direction.");
    elseif (per_x == true) && (per_y == false)
        % Let w(i) be the square of the x-periodic distance from seed i to the source domain (the c-transform of 0 at seed i)
        w = (- (bx(2)/2)*(sign(Z(:,2)-bx(2))-1)...
               + Z(:,2).*(1/2).*(sign(Z(:,2)-bx(2)) + sign(Z(:,2)-bx(4)))...
               - (bx(4)/2)*(sign(Z(:,2)-bx(4))+1)).^2;
    elseif (per_x == false) && (per_y == true)
        % Let w(i) be the square of the y-periodic distance from seed i to the source domain (the c-transform of 0 at seed i)
        w = (- (bx(1)/2)*(sign(Z(:,1)-bx(1))-1)...
               + Z(:,1).*(1/2).*(sign(Z(:,1)-bx(1)) + sign(Z(:,1)-bx(3)))...
               - (bx(3)/2)*(sign(Z(:,1)-bx(3))+1)).^2;
    elseif (per_x == true) && (per_y == true)
        w = zeros(N,1);
    end

end