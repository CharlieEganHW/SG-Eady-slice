function k = kper_bx(x,y,bx,p)
% This caculates k0 such that
% \| x - y + k0 \| = min_k \| x - y + k \|
% where the minimum is taken over p \mathbb{Z}^N

% Inputs are x  - a 1xN vector
%            y  - a 1xN vector
%            bx - a 1xN vector of box side lengths
%            p  - a 1xN vector of 0s and 1s    

[~,nx]=size(x);

%[mx,nx]=size(x);
%     if(mx~=1)
%         error('x should be a 1xN vector');
%     end
%     
%     [my,ny]=size(y);
%     if(my~=1)
%         error('y should be a 1xN vector');
%     end
%    
%     if(ny~=nx)
%         error('x and y should both be 1xN vectors');
%     end
%     
%     [mbx,nbx]=size(bx);
%     if(mbx~=1)
%         error('bx should be a 1X2N vector');
%     end
%     
%     if(nbx~=2*nx)
%         error('bx should be a 1x2N vector');
%     end
    
    % Diagonal matrix made of bx (box sizes)
    boxDims = [bx(3)-bx(1),bx(4)-bx(2)];
    Binv    = diag(p./boxDims);
    k       = floor((y-x)*Binv+0.5*ones(1,nx));
    
end 
