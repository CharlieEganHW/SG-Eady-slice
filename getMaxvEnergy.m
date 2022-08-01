function [maxv,totalEnergy,energyConservationErrors] = getMaxvEnergy(f,L,H,N,perL,perV,t,Z,w)
% 02/02/2022
%
% Function for computing the time evolution of the maximum meridional velocity
% (over all cell centroids), total energy and energy conservation relative
% error of a numerical solution of the Lagrangian SG Eady slice equations
% obtained using the geometric method.

% Input
%{
Format: variable - class; size; description.

Dimensional parameters

    f - double; 1 x 1; coriolis parameter
    L - double; 1 x 1; half-domain length
    H - double; 1 x 1; domain height
    N - double; 1 x 1; bouyancy frequency

Logical variables specifying domain periodicity

    perL - logical (either 'true' or 'false'); 1 x 1; specifies periodicity in longitudinal direction
    perV - logical (either 'true' or 'false'); 1 x 1; specifies periodicity in vertical

Data from numerical solution (n := number of seeds, numStepsRecorded := approximate number of time steps specified in solveSGEadySliceAdaptive.m)
    
    t - dobule; numStepsRecorded x 1; time points
    Z - double; n x 2 x numStepsRecorded; seed locations
    w - double; n x numStepsRecorded; optimal weights

%}

% Output
%{
%}

% Example 1
%{
% Unstable normal mode from Williams (1967) with parameters used in Egan
% et. al (2022) generated using getICForSteadyShearPerturbation.m
% Total run time is approximately 10 seconds

% Choice of perturbation

    perturbation = 'unstable'; % choose unstable normal mode from Williams (1967)

% Dimensional parameters

    g   = 10;     % acceleration due to gravity
    s   = -3e-6;  % latitudinal temperature gradient
    f   = 1e-4;   % coriolis parameter
    th0 = 3e+2;   % reference potential temperature
    L   = 1e+6;   % half-domain length
    H   = 1e+4;   % domain height (will be redefined to maximise growth rate of the perturbation)
    N   = 5e-3;   % bouyancy frequency
    a   = -7.5;   % amplitude of perturbation

% Simulation parameters

    numCols = 6;    % desired number of columns of seeds
    T       = 100;  % final time
    h       = 1;    % default time step size
    K       = 10;   % integer such that data is recorded at every K-th time point
    eta     = 1;    % percentage mass tolerance for optimal transport solver

% Get seeds masses and optimal weights

    [Z0,M,w0,H] = getICForSteadyShearPerturbation(perturbation,g,s,f,th0,L,H,N,a,numCols);

% Logical variables specifying domain periodicity

    perL = true;  % specifies periodicity in longitudinal direction
    perV = false; % specifies periodicity in vertical

% Solve ODE

    [Z,w,t,tolerances,halvings,runTime] = solveSGEadySliceAdaptive(g,s,f,th0,L,H,perL,perV,Z0,w0,M,T,h,K,eta);


% Now compute the RMSv, total energy and
% maximum v (at cell centroids) of the numerical solution

    [maxv,totalEnergy,energyConservationErrors] = getEnergyMaxv(f,L,H,N,perL,perV,t,Z,w);

%}

    % In solveSGEadySliceAdaptive.m, arrays have preallocated
    % sizes for efficiency. We first discard any unused entries.
    [~,numSteps] = max(t);
    
    % Define arrays to store results
    totalEnergy     = zeros(numSteps,1);
    maxv            = zeros(numSteps,1);
    
    % Rectangular fluid domain
    bx = [-L,-H/2,L,H/2];
    
    % Calculate the potential energy due to the linear steady temperature profile
    % By convention, in the SG model, the steady temperature is combined
    % with the perturbation to give the full temperature.
    % In previous works (in particular Yamazaki, Shipton,
    % Cullen, Marshall, Cotter (2017)), the steady temperature is not
    % included in the calculation of the potential energy.
    % For comparison with previous works, this should therefore be subtracted 
    % from the total energy computed using the full temperature from the SG
    % simulation.
    backgroundPotentialEnergy = -L*N^2*H^3/6;

    % Compute total energy and maximum v at each time point
    for k=1:numSteps
        
        Zk        = Z(:,:,k);
        [M,tc,xc] = mexPDall_2d(bx,Z(:,:,k),w(:,k),perL,perV);
        Zk        = getRemappedSeeds(bx,Zk,perL,perV);
        
        totalEnergyG       = f^2*sum(tc(:,1))/2 - f^2*sum(M.*Zk(:,2).^2)/2 - f^2*L*H^3/12;
        totalEnergy(k)     = totalEnergyG - backgroundPotentialEnergy;
        
        % maximum of v over all cell centroids
        maxv(k) = f*max(Zk(:,1)-xc(:,1));
        
    end
    
    % Compute energy conservation relative errors
    energyConservationErrors = zeros(numSteps,1);
    
    meanTotalEnergy = mean(totalEnergy);
    
    for k=1:numSteps
        
        energyConservationErrors(k) = (meanTotalEnergy-totalEnergy(k))/(meanTotalEnergy);
        
    end

end

