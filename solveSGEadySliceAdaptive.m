function [Z,w,t,tolerances,halvings,runTime] = solveSGEadySliceAdaptive(g,s,f,th0,L,H,Z0,M,T,h,K,eta)

% Function for solving the Lagrangian SG Eady slice equations with periodic 
% boundary conditions in the horizontal using the Geometric method of 
% Cullen & Purser (1984)
%
% Give initial 'seed' positions Z0 in 'geostrophic space', the ODEs for the
% seeds are solved using the two-step Adams-Bashforth method with adaptive 
% time-stepping. At each function evaluation a semi-dicrete optimal transport problem is 
% solved numerically using a damped Newton algorithm (see Merigot, Thibert (2019)). 
% To ensure convergence of this algorithm, the initial guess for the
% weights must not generate any zero-area cells. In practice, due to rounding 
% errors, it is necessary to ensure that all cells have sufficiently large area. 
% At each step, we approximate the optimal weights to first order. 
% If the subsequent cell areas are not all sufficiently large, the time-step size is halved. 
% To avoid excessively small time-steps, the time-step size reverts to a 
% default size after four unsuccesful halvings.
%
% Data is saved at every K-th time step for K specified in line 50.

% Input
%{
Format: variable - class; size; description.

Dimensional parameters

    g     - double; 1 x 1; acceleration due to gravity
    s     - double; 1 x 1; latitudinal temperature gradient
    f     - double; 1 x 1; coriolis parameter
    th0   - double; 1 x 1; reference potential temperature
    L     - double; 1 x 1; half-domain length
    H     - double; 1 x 1; domain height

Initial condition (n := number of seeds)

    Z0    - double; n x 2; array of initial seed positions
    w0    - double; n x 1; initial guess for the weights
    M     - double; n x 1; vector of cell masses

Simulation parameters

    T     - double; 1 x 1; final time
    h     - double; 1 x 1; default time step size
    K     - double; 1 x 1; integer such that data is recorded at every K-th time point
    eta   - double; 1 x 1; percentage mass tolerance for optimal transport solver

%}

% Output
%{

numStepsRecorded := number of steps at which data was recorded

Format: variable - class; size; description.

    Z          - double; n x 2 x numStepsRecorded; seed positions
    w          - double; n x numStepsRecorded; optimal weight vectors
    t          - double; numStepsRecorded x 1; time points
    tolerances - double; numStepsRecorded x 1; tolerances achieved
    halvings   - double; numStepsRecorded x 1; number of times the time step size was halved at each step
    runTime    - double; 1 x 1; time taken for code to run
%}

tic

%% Set up
bx        = [-L,-H/2,L,H/2]; % Define fundamental fluid domain

% Logical variables specifying domain periodicity
perL = true;  % specifies periodicity in longitudinal direction
perV = false; % specifies periodicity in vertical

timescale = g*s/f/th0;       % Define timescale

% Define matrices for storing data
numSteps = 2*ceil(T/h);              % Determine the approximate number of timesteps needed, taking into account that the time-step size may be reduced
numStepsRecorded = ceil(numSteps/K); % Determine the approximate number of timesteps at which to record data

n = length(M);                       % Number of seeds

Z          = nan(n,2,numStepsRecorded);  % 3D array to store seed positions
w          = nan(n,numStepsRecorded);    % 2D array to store weights 
tolerances = nan(1,numStepsRecorded);    % Vector to store tolerances achieved by the optimal transport solver
t          = [0,nan(1,numStepsRecorded-1)];    % Vector to store time points
halvings   = nan(1,numStepsRecorded);    % Vector to store number of times the default time-step is halved at each time-step

% Define the minimum desired cell area
% NB: This is used to avoid initialising the damped Newton algorithm with a weight vector that generates very small cells
% NB: (bx(3)-bx(1))*(bx(4)-bx(2))/n is the average target cell area
areaThreshold  = max([(1e-14)*2*L*H/n,1e-14]);

%% STEP k=1: apply 1st order Adams-Bashforth method (i.e., forward Euler) to compute Z(:,:,2)
k=1;

% Evaluate the right-hand side of the ODE
w0 = initialiseDampedNewton2d(bx,Z0,M,perL,perV,areaThreshold);                    % initial guess for weights
[vCurr,wk,tolAchieved] = getEadySliceVelocity(bx,Z0,M,w0,eta,timescale,perL,perV); % instantaneous velocity, remapped seeds, optimal weight vector, percentage mass tolerance achieved by optimal transport solver
tolerances(k)          = tolAchieved;
halvings(k)            = 0;

% Record remapped seed positions and optimal weights
Zk       = getRemappedSeeds(bx,Z0,perL,perV); % remap seeds into fundamental domain using function defined at the end of this script
Z(:,:,k) = Zk;
wk       = wk - wk(end); % use convention that n-th weight is zero
w(:,k)   = wk;

% Calculate increments for the seeds...
ZInc = h*vCurr;

% ...and weights such that the masses of the cells don't change to first 
% order in the magnitude of the velocities of the seeds
[~,DmDw,DmDz1,DmDz2] = getDm2d(bx,Zk,wk,perL,perV);   % calculate derivatives of the mass map with respect to the weights and the two seed components
DmDzTimesZInc        = DmDz1*ZInc(:,1) + DmDz2*ZInc(:,2);   % multiply the seed increment vector ZInc by the derivative of the mass map with respect to the seed locations
DmDwMod              = DmDw(1:n-1,1:n-1);                   % since we use the convention that w(end)=0, this (n-1) x (n-1) matrix is invertible
wInc                 = zeros(n,1);
wInc(1:end-1)        = -DmDwMod\DmDzTimesZInc(1:end-1);

% Calculate new seed positions and set the new guess for the weights
Zk = Zk + ZInc;
wk = wk + wInc;

h0 = h; % previous step size

% Check that this guess for the weights doesn't generate any cells with
% area less than the area threshold. If it does, use the default weight
% guess.
areas = mexPDall_2d(bx,Zk,wk,perL,perV);

if min(areas)<=areaThreshold
    wk = initialiseDampedNewton2d(bx,Zk,M,perL,perV,areaThreshold); % initial guess for weights
end

%% STEPS k=2:K : apply 2nd order Adams-Bashforth method to compute Z(:,:,k)
while max(t) <= T
    k = k+1;
    
    % Index at which to record data
    % NB: the sequence of indices thus generated is 
    % (2,...,2,3,...,3,4,...,4,...). Data is therefore overwritten K-many 
    % times before being recorded and retained, as desired.
    ind = ceil((k-1)/K)+1;
    
    % Store the k-th timepoint and reset the proposed step size
    t(ind) = max(t) + h0;
    h1     = h;
    
    % Update velocities
    vPrev = vCurr;
    
    %Solve OT problem to find geostrophic velocity
    [vCurr,wk,tolAchieved] = getEadySliceVelocity(bx,Zk,M,wk,eta,timescale,perL,perV);
    tolerances(ind)      = tolAchieved;
    
    % Record remapped seed positions, optimal weights and optimal centroids
    Zk         = getRemappedSeeds(bx,Zk,perL,perV); % remap seeds into fundamental domain using function defined at the end of this script
    Z(:,:,ind) = Zk;
    wk         = wk - wk(end); % use convention that Nth weight is zero
    w(:,ind)   = wk;
    
    % Calculate increments for the seeds...
    cPrev = -(1/2)*(h1^2/h0);
    cCurr = (1/2)*((h1+h0)^2/h0 - h0);
    ZInc  = cPrev*vPrev + cCurr*vCurr;

    % ...and weights such that the masses of the cells don't change to first 
    % order in the magnitude of the velocities of the seeds
    [~,DmDw,DmDz1,DmDz2] = getDm2d(bx,Zk,wk,perL,perV);
    DmDzTimesZInc       = DmDz1*ZInc(:,1) + DmDz2*ZInc(:,2);
    DmDwMod         = DmDw(1:n-1,1:n-1);                 
    wInc(1:end-1)  = -DmDwMod\DmDzTimesZInc(1:end-1);
    
    % Calculate new seed positions and set the new guess for the weights
    ZNew = Zk + ZInc;
    wNew = wk + wInc;
    
    % Find appropriate step size: if minimum area is too small, halve the
    % step size and recompute the velocity
    areas = mexPDall_2d(bx,ZNew,wNew,perL,perV);
    
    halving = 0;
    
    while min(areas) <= areaThreshold && h1 >= h/2^7
        % halve the proposed step-size
        h1      = h1/2;
        halving = halving + 1;
        
        % Calculate increments for seeds...
        cPrev = -(1/2)*(h1^2/h0);
        cCurr = (1/2)*((h1+h0)^2/h0 - h0);
        ZInc  = cPrev*vPrev + cCurr*vCurr;

        % ...and weights such that the masses of the cells don't change to first 
        % order in the magnitude of the velocities of the seeds
        [~,DmDw,DmDz1,DmDz2] = getDm2d(bx,Zk,wk,perL,perV);
        DmDzTimesZInc        = DmDz1*ZInc(:,1) + DmDz2*ZInc(:,2);
        DmDwMod              = DmDw(1:n-1,1:n-1);                 
        wInc(1:end-1)        = -DmDwMod\DmDzTimesZInc(1:end-1);
    
        % Calculate new seed positions and set the new guess for the weights
        ZNew = Zk + ZInc;
        wNew = wk + wInc;
    
        % Calculate areas generated by Znew and wnew
        areas = mexPDall_2d(bx,ZNew,wNew,perL,perV);
    end
    
    % If, after seven halvings of the original time step (i.e., if h1 =
    % h/2^7), wnew still generates cells whose area is less than the area
    % tolerance, reset the time step size to the default h.
    % NB: the performance of the damped Newton algorithm when using the
    % first order approximation of the optimal weights is not
    % sufficiently improved to justify halving the time step size further.
    if min(areas) <= areaThreshold
        % Reset the time step size
        h1 = h;
        
        % Calculate increments for the seeds
        cPrev = -(1/2)*(h1^2/h0);
        cCurr = (1/2)*((h1+h0)^2/h0 - h0);
        ZInc  = cPrev*vPrev + cCurr*vCurr;
    
        % Calculate new seed positions and set the new guess for the
        % weights using initialiseDampedNewton
        ZNew = Zk + ZInc;
        wNew = initialiseDampedNewton2d(bx,ZNew,M,perL,perV,areaThreshold); % initial guess for weights
        
    end
    
    % Record the number of times the time step size was halved
    halvings(ind) = halving;
    
    h0 = h1;
    Zk = ZNew;
    wk = wNew;
    
end

%% Final step: Compute the optimal weights and remapped seeds at the final timestep
k   = k+1;
ind = ceil((k-1)/K)+1;

[w(:,ind),tolAchieved] = dampedNewton2d(bx,Zk,M,wk,eta,perL,perV);
tolerances(ind)        = tolAchieved;

Zk         = getRemappedSeeds(bx,Zk,perL,perV); % remap seeds into fundamental domain using function defined at the end of this script
Z(:,:,ind) = Zk;

runTime = toc;

%% Delete any unnecessary entries of outputs
stepsRecorded = ~isnan(t);
Z             = Z(:,:,stepsRecorded);  
w             = w(:,stepsRecorded);  
tolerances    = tolerances(stepsRecorded);  
t             = t(stepsRecorded);    
halvings      = halvings(stepsRecorded);  

%% Save relevant variables to a Version 7.3 MAT-File
clear vfn tolAchieved ind k wk Zk h0 wNew ZNew ZInc cCurr cPrev h1 areas ...
      wInc DmDw DmDwMod DmDzTimesZInc DmDz2 DmDz1 DmDz halving vCurr vPrev ...
      numStepsRecorded numSteps K stepsRecorded
%save('EadySliceSimulationResults','-v7.3')
end

function [F,w,tolAchieved] = getEadySliceVelocity(bx,Z,M,w0,eta,timescale,perL,perV)
%%
% 01/02/22
% Function to evaluate the right hand side of the ODE

% Input
%{
n := number of seeds

Format: variable - class; size; description

    bx        - double; 1 x 4; [xmin ymin xmax ymax] specifying rectangular domain
    Z         - double; 1 x 1; seed locations
    M         - double; 1 x 1; target masses
    w0        - double; 1 x 1; guess for optimal weights
    eta       - double; 1 x 1; percentage mass tolerance for optimal transport solver
    timescale - double; 1 x 1; timescale
    perL      - logical (either true or false); 1 x 1; specifies periodicity in longitudinal direction
    perV      - logical (either true or false); 1 x 1; specifies periodicity in vertical

%}
% Output
%{
n := number of seeds

Format: variable - class; size; description

   F     - double; n x 2; for i in {1,...,n}, [F(i,1),F(i,2)] is the time derivative of the i-th seed
   w     - double; n x 1; optimal weight vector
   etaAchieved - dobule; 1 x 1; percentage mass tolerance achieved by optimal transport solver
%}

n = length(Z(:,1));  % number of Dirac masses

% Solve the semi-discrete OT problem to find the weights using the damped
% Newton algorithm of Kitagawa, Merigot and Thibert (2019)
[w,tolAchieved] = dampedNewton2d(bx,Z,M,w0,eta,perL,perV);

% Find the centroids of the Laguerre cells
[~,~,centroids] = mexPDall_2d(bx,Z,w,perL,perV);

% Compute the right hand side of the ODE
F      = zeros(n,2);
F(:,1) = -timescale*centroids(:,2);         
F(:,2) = timescale*(centroids(:,1)-Z(:,1));

end