function stillsDomainCellsEadySlice(Z,w,perL,perV,g,f,th0,L,H,N,t,steps,colouring)
% Function to produce figures of fluid domain coloured by potential temperature, 
% meridional velocity or perturbation of potential temperature from shear
% flow steady state

% Input:
%{
    Z              - N x 2 x (K+1) 3d array containing seed positions
    w              - N x (K+1) matrix of weights
    perL          - periodicity in x-direction
    perV          - periodicity in y-direction
    g              - acceleration due to gravity (meters per second^2)
    f              - Coriolis parameter (seconds^{-1})
    th0            - background temperature (Kelvin)
    L              - half domain length (meters)
    H              - domain height (meters)
    N              - Buoyancy frequency (seconds^{-1})
    t              - 1 x (K+1) row vector of times
    steps          - time steps to be plotted
    frames_per_day - desired number of frames per day
    colouring      - one of three 4-letter strings, 'mvel', 'temp' or 'thPe', specifying colouring of cells
%}

    
    %% Make steps a row vector so that it can be used for indexing
    steps = steps(:)';
    
    num_steps = length(steps(:));
    
    %% Define the fluid domain
    bx = [-L,-H/2,L,H/2];
    
    %% Make figure
    if prod(colouring == 'mvel') == 1
        %% Colour by meridional velocity
        
        % Set colours of cells and variables defining the colour bar
        % NB: this invloves computing Laguerre tessellations so we save the
        % data to pass to the plotting function later
        colour_label           = 'Meridional velocity $\mathrm{ms}^{-1}$';
        meriodional_velocities = zeros(size(Z,1),num_steps);
        vfn_all                = cell(num_steps,1);
    
        for k=steps
            Zk               = Z(:,:,k);
            wk               = w(:,k);
            [~,~,xck,vfnk]   = mexPDall_2d(bx,Zk,wk,perL,perV);
            mer_vels_k       = f*(Zk(:,1)-xck(:,1));
        
            meridional_velocities(:,k) = mer_vels_k;
            vfn_all{k}                 = vfnk;
        end
        
        max_vel = max(meridional_velocities,[],'all'); % upper limit for colour bar
        min_vel = min(meridional_velocities,[],'all'); % lower limit for colour bar

        denom = 50; % distance between ticks on the colourbar

        colours    = cell(3,1);
        colours{1} = colour_label;
        colours{2} = [denom*floor(min_vel/denom),denom*ceil(max_vel/denom),denom];

        % Plot coloured cells at each selected step, save figure and export
        % graphic
        for k=steps
            time       = t(k)/60/60/24;
            colours{3} = meridional_velocities(:,k);
            vfn        = vfn_all{k};
            
            clf
            
            cellPlot = tiledlayout(1,1,'Padding','none');
            cellPlot.Units = 'inches';
            cellPlot.OuterPosition = [0.25 0.25 2.5 2];
            nexttile;

            plotCellsEadySlice(bx,vfn,colours,time);
        
            %saveas(gcf,strcat('merVelStep',num2str(k),'.fig'));
            %exportgraphics(gcf,strcat('merVelStep',num2str(k),'.png'),'Resolution',300);
            
        end
        
    elseif prod(colouring == 'temp') == 1
        %% Colour by temperature
        colour_label = 'Deviation from reference temperature (Kelvin)';
    
        % Set colours of cells and variables defining the colour bar
        c = (th0*f^2)/g;
        max_temp = c*max(Z(:,2,:),[],'all'); % upper limit for colour bar
        min_temp = c*min(Z(:,2,:),[],'all'); % lower limit for colour bar
        
        denom    = 10; % distance between ticks on the colourbar

        colours    = cell(3,1);
        colours{1} = colour_label;
        colours{2} = [min_temp,max_temp,denom];
        
        % Plot coloured cells at each selected step, save figure and export
        % graphic
        for k=steps
            time        = t(k)/60/60/24;
            Zk          = Z(:,:,k);
            wk          = w(:,k);
            [~,~,~,vfn] = mexPDall_2d(bx,Zk,wk,perL,perV);
            colours{3}  = c*Zk(:,2);
            
            clf
            
            cellPlot = tiledlayout(1,1,'Padding','none');
            cellPlot.Units = 'inches';
            cellPlot.OuterPosition = [0.25 0.25 4 3];
            nexttile;

            plotCellsEadySlice(bx,vfn,colours,time);
            
            %saveas(gcf,strcat('tempStep',num2str(k),'.fig'));
            %exportgraphics(gcf,strcat('tempStep',num2str(k),'.png'),'Resolution',300);
        end
    elseif prod(colouring == 'thPe') == 1
        %% Colour by temperature perturbation from steady state
        
        %colour_label = 'Deviation from steady temperature (Kelvin)';
        colour_label = 'Kelvin';
        
        % Set colours of cells and variables defining the colour bar
        % find the maximum and minimum temperature perturbations at given
        % time-points
        
        k          = steps(1);
        Zk         = Z(:,:,k);
        wk         = w(:,k);
        [~,~,xc,~] = mexPDall_2d(bx,Zk,wk,perL,perV);
        temp_perts = th0*f^2/g*(Zk(:,2)-N^2/f^2*(xc(:,2)+H/2));

        min_temp_pert = min(temp_perts);
        max_temp_pert = max(temp_perts);
        
        for k=steps(2:end)
            Zk         = Z(:,:,k);
            wk         = w(:,k);
            [~,~,xc,~] = mexPDall_2d(bx,Zk,wk,perL,perV);
            temp_perts = th0*f^2/g*(Zk(:,2)-N^2/f^2*(xc(:,2)+H/2));
            
            min_temp_pert = min([min_temp_pert;temp_perts]);
            max_temp_pert = max([max_temp_pert;temp_perts]);
        end
        
        denom    = 2; % distance between ticks on the colourbar

        colours    = cell(3,1);
        colours{1} = colour_label;
        colours{2} = [min_temp_pert,max_temp_pert,denom];
        
        % Plot coloured cells at each selected time-point and write to video
        % file
        %for i=17:20 % use this to plot days 4-5 with frames_per_day = 4
        for k=steps
            time        = t(k)/60/60/24;
            Zk          = Z(:,:,k);
            wk          = w(:,k);
            [~,~,xc,vfn] = mexPDall_2d(bx,Zk,wk,perL,perV);
            colours{3}  = th0*f^2/g*(Zk(:,2)-N^2/f^2*(xc(:,2)+H/2));
            
            cellPlot = tiledlayout(1,1,'Padding','none');
            cellPlot.Units = 'inches';
            cellPlot.OuterPosition = [0.25 0.25 3 3];
            nexttile;
            
            plotCellsEadySlice(bx,vfn,colours,time);
            
            %saveas(gcf,strcat('tempPertStep',num2str(k),'.fig'));
            %exportgraphics(gcf,strcat('tempPertStep',num2str(k),'.png'),'Resolution',300);
            
        end
    else
        error('Unrecognised colour code. Use strings "temp", "thPe", or "mvel" to specify colouring.')
    end
   
end