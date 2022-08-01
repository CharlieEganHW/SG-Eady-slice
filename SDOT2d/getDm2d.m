function [areas,H,Hz1,Hz2] = getDm2d(bx,X,w,per_x,per_y)

% [areas,H,Hz1,Hz2] = get_g(w,X,bx,per)
%
% Input arguments
%         w            - is the weights an Nx1 array
%         X            - is the positions, Nx2 array
%         bx           - box size, a 1x4 array
%         per_x        - periodic flag (a boolean true/false to indicate periodic in x or not)
%         per_y        - periodic flag (a boolean true/false to indicate periodic in y or not)
%
% Return arguments

%         H            - the Hessian d^2g, a sparse matrix, dm/dw
%         Hz1          - sparse matrix dm/dz, 1st component
%         Hz2          - sparse matrix dm/dz, 2nd component
    
    %% Catch errors
    [Nw,~]=size(w);    

    %[Nw,Mw]=size(w);    
    %{
    if(Mw~=1)
        error('w should be an N x 1 array where N is the number of cells');
    end
    %}
    %[NX,MX]=size(X);
    %{
    if(MX~=2)
        error('X should be an N x 2 array where N is the number of cells');
    end
    %}
    %{
    if(NX~=Nw)
        error('The number of cells represented by w and X disagree, w: %d and X %d',Nw,NX);
    end
    %}
    Lx=bx(3)-bx(1);
    Ly=bx(4)-bx(2);
    
    % Get copies of seeds in the fundamental domain
    X = getRemappedSeeds(bx,X,per_x,per_y);
    
    %% Computations
    
    % Use mexPD to calculate the Laguerre diagram
    [areas,~,~,vfn]=mexPDall_2d(bx,X,w,per_x,per_y);
    
    % Matrices to store the values and indices for the sparse Hessian matrix
    H_spvals_i=zeros(1,Nw^2);
    H_spvals_j=zeros(1,Nw^2);
    H_spvals_val=zeros(1,Nw^2);
    
    Hz1_spvals_i=zeros(1,Nw^2);
    Hz1_spvals_j=zeros(1,Nw^2);
    Hz1_spvals_val=zeros(1,Nw^2);
    
    Hz2_spvals_i=zeros(1,Nw^2);
    Hz2_spvals_j=zeros(1,Nw^2);
    Hz2_spvals_val=zeros(1,Nw^2);
    
    % Index to keep track of the number of entries in the sparse matrix
    r=1;rz=1;
    
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
        edge_centroids_cell_i = 0.5*(vertices_cell_i + circshift(vertices_cell_i,-1));
        
        % The number of neighbours - boundaries have negative indices
        Num_N_i = size(N_i);
        
        % calculate H_ij for j in N_i 
        for j=1:Num_N_i
            
            k = N_i(j); % k is the j-th neighbour of cell i
            
            if(k>0)
                
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
                
                % centroid of the edge
                x_ik = edge_centroids_cell_i(j,:);
                
                % Distance between the seeds
                % if periodic del_k = del_k and if not del_k = 0    
                
                xi_minus_xk=xi-xk+del_k*[Lx,0;0,Ly];
                dist_ik = norm(xi_minus_xk);
                
                % Values of the Hessian
                
                H_spvals_i(r)=i;
                H_spvals_j(r)=k;
                H_spvals_val(r)=-0.5*l_ik/dist_ik;
                r=r+1;
                
                % Hz
                
                entry=-l_ik*((x_ik-xi)+xi_minus_xk)/dist_ik;
                diag_entry=l_ik*(x_ik-xi)/dist_ik;
                
                % First component
                
                Hz1_spvals_i(rz)=i;
                Hz1_spvals_j(rz)=k;
                Hz1_spvals_val(rz)=entry(1);
                
                Hz1_spvals_i(rz+1)=i;
                Hz1_spvals_j(rz+1)=i;
                Hz1_spvals_val(rz+1)=diag_entry(1);
                
                % Second component
                
                Hz2_spvals_i(rz)=i;
                Hz2_spvals_j(rz)=k;
                Hz2_spvals_val(rz)=entry(2);
                
                Hz2_spvals_i(rz+1)=i;
                Hz2_spvals_j(rz+1)=i;
                Hz2_spvals_val(rz+1)=diag_entry(2); 
                
                rz=rz+2;
            end            
        end
    end
    
    % Discard extra entries of vectors with preallocated size
    num_entries_H  = find(H_spvals_i==0,1,'first')-1;
    num_entries_Hz = find(Hz1_spvals_j==0,1,'first')-1;

    % Now use H_spvals to construct the sparse matrix
    H=sparse(H_spvals_i(1:num_entries_H),H_spvals_j(1:num_entries_H),H_spvals_val(1:num_entries_H),Nw,Nw);
    
    Hz1=sparse(Hz1_spvals_i(1:num_entries_Hz),Hz1_spvals_j(1:num_entries_Hz),Hz1_spvals_val(1:num_entries_Hz),Nw,Nw);
    Hz2=sparse(Hz2_spvals_i(1:num_entries_Hz),Hz2_spvals_j(1:num_entries_Hz),Hz2_spvals_val(1:num_entries_Hz),Nw,Nw);
    
    % Now having all the off-diagonal entries we can add the diagonal entries
    H=H-spdiags(sum(H)',0,Nw,Nw);
    
end





