function animateDomainCellsEadySlice(Z,w,perL,perV,g,f,th0,L,H,N,t,frames_per_day,colouring)
% Function to produce animation of fluid domain coloured by potential temperature, 
% meridional velocity or perturbation of potential temperature from shear
% flow steady state.

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
    frames_per_day - desired number of frames per day
    colouring      - one of three 4-letter strings, 'mvel', 'temp' or 'thPe', specifying colouring of cells
%}
    %% Suppress output of figures at each step
    set(0,'DefaultFigureVisible','off');
    
    %% Set the maximum number of days to animate (comment out if animating full results)
    %num_days_to_animate = 10;
    %final_index         = find(t>=num_days_to_animate*24*60*60,1,'first');
    %t                   = t(1:final_index);
    
    %% Define the fluid domain
    bx = [-L,-H/2,L,H/2];
    
    %% Make animation
    if prod(colouring == 'mvel') == 1
        %% Colour by meridional velocity
        % Set up video-writer
        videoName   = ['cellsMerVels',num2str(frames_per_day),'FramesPerDay.avi']; % define the name of the video
        v           = VideoWriter(videoName);
        frames      = ceil(frames_per_day*max(t)/(24*60*60))-1; % number of frames at frames_per_day frames per day
        v.FrameRate = frames_per_day;  % 1 second in video is 1 day of siulation

        open(v);
        
        % Select time-points to plot
        indices = zeros(frames,1);
        indices(1)=1;
        for i=2:frames
            indices(i) = find(t>=(i-1)*24*60*60/frames_per_day,1,'first');
        end
        
        % Set colours of cells and variables defining the colour bar
        % NB: this invloves computing Laguerre tessellations so we save the
        % data to pass to the plotting function later
        colour_label           = 'Meridional velocity $\mathrm{ms}^{-1}$';
        meridional_velocities = zeros(size(Z,1),frames);
        vfn_all                = cell(frames,1);
    
        for i=1:frames
            k                = indices(i);
            Zk               = Z(:,:,k);
            wk               = w(:,k);
            [~,~,xck,vfnk]   = mexPDall_2d(bx,Zk,wk,perL,perV);
            mer_vels_k       = f*(Zk(:,1)-xck(:,1));
        
            meridional_velocities(:,i) = mer_vels_k;
            vfn_all{i}                 = vfnk;
        end
        
        max_vel = max(meridional_velocities,[],'all'); % upper limit for colour bar
        min_vel = min(meridional_velocities,[],'all'); % lower limit for colour bar

        denom = 20; % distance between ticks on the colourbar

        colours    = cell(3,1);
        colours{1} = colour_label;
        colours{2} = [denom*floor(min_vel/denom),denom*ceil(max_vel/denom),denom];

        % Plot coloured cells at each selected time-point and write to video
        % file
        for i=1:frames
            i
            time        = (i-1)/frames_per_day;
            colours{3} = meridional_velocities(:,i);
            vfn        = vfn_all{i};
            
            clf

            plotCellsEadySlice(bx,vfn,colours,time);
                       
            % write to video
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        
    elseif prod(colouring == 'temp') == 1
        %% Colour by temperature
        % Set up video-writer
        videoName   = ['cellsTemp',num2str(frames_per_day),'FramesPerDay.avi']; % define the name of the video
        v           = VideoWriter(videoName);
        frames      = ceil(frames_per_day*max(t)/(24*60*60))-1; % number of frames at frames_per_day frames per day
        v.FrameRate = frames_per_day;                           % 1 second in video is 1 day of siulation

        open(v);
        
        % Select time-points to plot
        indices = zeros(frames,1);
        indices(1)=1;
        for i=2:frames
            indices(i) = find(t>=i*24*60*60/frames_per_day,1,'first');
        end
        colour_label = 'Deviation from reference temperature (Kelvin)';
    
        % Set colours of cells and variables defining the colour bar
        c = (th0*f^2)/g;
        max_temp = c*max(Z(:,2,:),[],'all'); % upper limit for colour bar
        min_temp = c*min(Z(:,2,:),[],'all'); % lower limit for colour bar
        
        denom    = 10; % distance between ticks on the colourbar

        colours    = cell(3,1);
        colours{1} = colour_label;
        colours{2} = [min_temp,max_temp,denom];
        
        % Plot coloured cells at each selected time-point and write to video
        % file
        for i=1:frames
            i
            time        = (i-1)/frames_per_day;
            k           = indices(i);
            Zk          = Z(:,:,k);
            wk          = w(:,k);
            [~,~,~,vfn] = mexPDall_2d(bx,Zk,wk,perL,perV);
            colours{3}  = c*Zk(:,2);

            plotCellsEadySlice(bx,vfn,colours,time);

            frame = getframe(gcf);
            writeVideo(v,frame);
        end
    elseif prod(colouring == 'thPe') == 1
        %% Colour by temperature perturbation from steady state
        % Set up video-writer
        videoName   = ['cellsTempPert',num2str(frames_per_day),'FramesPerDay.avi']; % define the name of the video
        v           = VideoWriter(videoName);
        frames      = ceil(frames_per_day*max(t)/(24*60*60))-1; % number of frames at frames_per_day frames per day
        v.FrameRate = frames_per_day;                           % 1 second in video is 1 day of siulation

        open(v);
        
        % Select time-points to plot
        indices = zeros(frames,1);
        indices(1)=1;
        for i=2:frames
            indices(i) = find(t>=i*24*60*60/frames_per_day,1,'first');
        end
        colour_label = 'Deviation from steady temperature (Kelvin)';
        
        % Set colours of cells and variables defining the colour bar
        
        % find the maximum and minimum temperature perturbations at given
        % time-points
        
        k          = indices(1);
        Zk         = Z(:,:,k);
        wk         = w(:,k);
        [~,~,xc,~] = mexPDall_2d(bx,Zk,wk,perL,perV);
        temp_perts = th0*f^2/g*(Zk(:,2)-N^2/f^2*(xc(:,2)+H/2));

        min_temp_pert = min(temp_perts);
        max_temp_pert = max(temp_perts);
        
        for i=2:frames
            k          = indices(i);
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
        for i=1:frames
            i
            time        = (i-1)/frames_per_day;
            k           = indices(i);
            Zk          = Z(:,:,k);
            wk          = w(:,k);
            [~,~,xc,vfn] = mexPDall_2d(bx,Zk,wk,perL,perV);
            colours{3}  = th0*f^2/g*(Zk(:,2)-N^2/f^2*(xc(:,2)+H/2));

            plotCellsEadySlice(bx,vfn,colours,time);
            
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
    else
        error('Unrecognised colour code. Use strings "temp", "thPe", or "mvel" to specify colouring.')
    end
   
end