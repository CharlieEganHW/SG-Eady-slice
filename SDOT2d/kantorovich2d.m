function [g,Dg,H,actual_areas] = kantorovich2d(bx,X,target_areas,w,per_x,per_y)

% 21/12/20
% [g,Dg,H,actual_vols] = get_g(w,X,target_vols,bx,periodic)
%
% Input arguments
%         w            - is the weights an Nx1 array
%         X            - is the positions, three Nx3 array
%         target_vols  - target volumes an Nx1 array
%         bx           - box size, a 1x3 array
%         per_x        - periodic flag (a boolean true/false to indicate periodicity in x-direction)
%         per_y        - periodic flag (a boolean true/false to indicate periodicity in y-direction)
%
% Return arguments
%         g            - function g(w;x)
%         Dg           - gradient of g wrt w
%         H            - the Hessian d^2g, a sparse matrix
%         actual_vols  - the actual volumes of the Laguerre diagram with seeds X and weights w
%  Demo_layers = 5.4 and data fitting = 5.3 
    
    %% Catch errors
    [Nw,~]=size(w);  
    %{
    [Nw,Mw]=size(w);    
    
    if(Mw~=1)
        error('w should be an N x 1 array where N is the number of cells');
    end
    
    [NX,MX]=size(X);
    
    if(MX~=2)
        error('X should be an N x 2 array where N is the number of cells');
    end

    if(NX~=Nw)
        error('The number of cells represented by w and X disagree, w: %d and X %d',Nw,NX);
    end
    
    [Ntv,Mtv]=size(target_areas);
    
    if(Mtv~=1)
        error('target_vols should be an N x 1 array where N is the number of cells');
    end

    if(Ntv~=Nw)
        error('The number of cells represented by w and target_areas disagree, w: %d and target_areas: %d',Nw,Ntv);
    end
    
    total_area=sum(target_areas);
    bx_area=(bx(3)-bx(1))*(bx(4)-bx(2));
    
    if(abs(total_area-bx_area)/bx_area>1e-10)
        error('The sum of the target volumes is different to the volume of the box');
    end
    %}
    %% Computations
    
    % Use mexPD to calculate the Laguerre diagram
    [actual_areas,transport_costs,~,vfn]=mexPDall_2d(bx,X,w,per_x,per_y);
    
    % Gradient of g
    Dg = actual_areas-target_areas;
    
    % Definition of g, a convex function of the weights 
    g = dot(Dg,w)-sum(transport_costs); 
    
    % Vectors to store the values and indices for the sparse Hessian matrix
    H_spvals_i=zeros(1,Nw^2);
    H_spvals_j=zeros(1,Nw^2);
    H_spvals_val=zeros(1,Nw^2);
    
    % Index to keep track of the number of entries in the sparse matrix
    r=1;
    
    % we now start computing the Hessian values
    for i = 1:Nw
        
        % Coordinates of the seed location of the ith cell
        xi = X(i,:);
        
        % The neighbours of the ith cell 
        N_i = vfn{i,2};
        
        % The vertices for the ith cell
        vertices_cell_i = vfn{i,1};
        
        edge_Lengths_cell_i = get_edgeLengths(vertices_cell_i);
        edge_midpoints_cell_i = 0.5*(vertices_cell_i + circshift(vertices_cell_i,-1));
        
        % The number of neighbours - boundaries have negative indices
        Num_N_i = size(N_i);
        
        % calculate H_ij for j in N_i 
        for j=1:Num_N_i
            
            k = N_i(j); % k is the j-th neighbour of cell i
            
            % we only need to calculate if it is a true neighbour, i.e.,
            if (k>i)  % upper triangular elements only
                
                del_k = zeros(1,length(xi));
                xk = X(k,:); % coordinates of seed k
                %if(per_x || per_y)
                    % 1. Calculate Centroid of the face
                    point_on_edge = edge_midpoints_cell_i(j,:);
                    
                    % 2. Calculate K(x,y) and K(x,z) from x = centroid
                    k_xk = kper_bx(point_on_edge,xk,bx,[per_x,per_y]);
                    k_xi = kper_bx(point_on_edge,xi,bx,[per_x,per_y]);
                    
                    % 3. del_k = K(x,y) - K(x_z)
                    del_k = k_xk - k_xi;
                %end
                % length of the edge 
                l_ik = edge_Lengths_cell_i(j);
                
                % Distance between the seeds
                % if periodic del_k = del_k and if not del_k = 0    
                dist_ik = norm(xi-xk+del_k*[bx(3)-bx(1),0;0,bx(4)-bx(2)]);
                
                % Value of the Hessian
                H_spvals_i(r)=i;
                H_spvals_j(r)=k;
                H_spvals_val(r)=-0.5*l_ik/dist_ik;
                r=r+1;
            end            
        end
        
    end
    
    % Discard extra entries of vectors with preallocated size
    num_entries = find(H_spvals_i==0,1,'first')-1;
    
    % Now use H_spvals to construct the sparse matrix
    H=sparse(H_spvals_i(1:num_entries),H_spvals_j(1:num_entries),H_spvals_val(1:num_entries),Nw,Nw);
    
    % Having calculated the upper triangular part we can find all off-diagonal entries by adding the transpose
    H=H+H';

    % Now having all the off-diagonal entries we can add the diagonal entries
    H=H-spdiags(sum(H)',0,Nw,Nw);
    
end
                
               
                
                
               
               
                