function animateSeedsEadySlice(Z,w,perL,perV,g,f,th0,L,H,t,frames_per_day,colouring)
% Function to produce animation of seeds coloured either by potential temperature
% or meridional velocity.

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
    colouring      - one of two 4-letter strings, 'mvel' or 'temp', specifying colouring of cells
%}
 
    %% Define the fluid domain
    bx = [-L,-H/2,L,H/2];
    
    %% Get limits for geostrophic domain
    ymin = min(Z(:,2,:),[],'all');
    ymax = max(Z(:,2,:),[],'all');
    
    %% Suppress output of figures at each step
    set(0,'DefaultFigureVisible','off');
    
    %% Create animation
    if prod(colouring == 'mvel') == 1
        %% Colour by meridional velocity
        % Set up video-writer
        v           = VideoWriter(strcat('Meridional_velocity_seeds_',num2str(frames_per_day),'_frames_per_day','.avi'));
        frames      = ceil(frames_per_day*max(t)/(24*60*60))-1; % number of frames at frames_per_day frames per day
        v.FrameRate = frames_per_day;  % 1 second in video is 1 day of siulation

        open(v);
        
        % Select time-points to plot
        indices = zeros(frames,1);
        indices(1)=1;
        for i=2:frames
            indices(i) = find(t>=i*24*60*60/frames_per_day,1,'first');
        end
        
        % Set colours of seeds and variables defining the colour bar
        % NB: this invloves computing Laguerre tessellations so we save the
        % data to pass to the plotting function later
        colour_label           = 'Meriodional velocity ms^{-1}';
        meriodional_velocities = zeros(size(Z,1),frames);
    
        for i=1:frames
            k           = indices(i);
            Zk          = Z(:,:,k);
            wk          = w(:,k);
            [~,~,xck,~] = mexPDall_2d(bx,Zk,wk,perL,perV);
            mer_vels_k  = f*(Zk(:,1)-xck(:,1));
        
            meridional_velocities(:,i) = mer_vels_k;
        end
        
        max_vel = max(meridional_velocities,[],'all'); % upper limit for colour bar
        min_vel = min(meridional_velocities,[],'all'); % lower limit for colour bar

        denom = 40; % distance between ticks on the colourbar

        colours    = cell(3,1);
        colours{1} = colour_label;
        colours{2} = [min_vel,max_vel,denom];

        % Plot coloured seeds at each selected time-point and write to video
        % file
        for i=1:frames
            figure;
            clf;
            i
            
            % Plot the coloured seeds
            k          = indices(i);
            Zk         = Z(:,:,k);
            colours{3} = meridional_velocities(:,i);
            
            scatter(Zk(:,1),Zk(:,2),3.2,colours{3},'filled')
            
            % Set the title
            day  = ceil(i/frames_per_day);
            title(['Day',' ',num2str(day)],'Interpreter','latex')
            
            % Set the axis limits
            xlim([bx(1),bx(3)]);
            ylim([ymin,ymax]);
            
            % set up colourbar with ticks at intervals specified by colours{2}(3) 
            % and spanning range specified by colours{2}(1:2)
            denom = colours{2}(3);
            ticks = denom*(ceil(colours{2}(1)/denom):floor(colours{2}(2)/denom));
            cbr   = colorbar('Ticks',ticks);
            caxis([ticks(1),ticks(end)]);
            cbr.Label.String = colours{1};
            
            % Write to video file
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        
    elseif prod(colouring == 'temp') == 1
        %% Colour by temperature
        % Set up video-writer
        v           = VideoWriter(strcat('Temperature_seeds_',num2str(frames_per_day),'_frames_per_day','.avi'));
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
    
        % Set colours of seeds and variables defining the colour bar
        c = (th0*f^2)/g;
        max_temp = c*max(Z(:,2,:),[],'all'); % upper limit for colour bar
        min_temp = c*min(Z(:,2,:),[],'all'); % lower limit for colour bar
        
        denom    = 10; % distance between ticks on the colourbar

        colours    = cell(3,1);
        colours{1} = colour_label;
        colours{2} = [min_temp,max_temp,denom];
        
        % Plot coloured seeds at each selected time-point and write to video
        % file

        for i=1:frames
            i
            figure;
            clf;
            
            % Plot the coloured seeds
            k          = indices(i);
            Zk         = Z(:,:,k);
            colours{3} = c*Zk(:,2);
            
            scatter(Zk(:,1),Zk(:,2),[],colours{3},'filled')
            
            % Set the title
            day  = ceil(i/frames_per_day);
            title(['Day',' ',num2str(day)])
            
            % Set the axis limits
            xlim([bx(1),bx(3)]);
            ylim([ymin,ymax]);
            
            % set up colourbar with ticks at intervals specified by colours{2}(3) 
            % and spanning range specified by colours{2}(1:2)
            denom = colours{2}(3);
            ticks = denom*(ceil(colours{2}(1)/denom):floor(colours{2}(2)/denom));
            cbr   = colorbar('Ticks',ticks);
            caxis([ticks(1),ticks(end)]);
            cbr.Label.String = colours{1};
            
            % Write to video file
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
    else
        error('Unrecognised 4-digit colour code. Use strings "temp" or "mvel" to specify colouring.')
    end
       
end