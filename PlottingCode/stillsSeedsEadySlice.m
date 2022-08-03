function stillsSeedsEadySlice(Z,w,perL,perV,g,f,th0,L,H,steps,colouring)
% Function to produce figures of seeds coloured by potential temperature, 
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
    steps          - time steps to be plotted
    frames_per_day - desired number of frames per day
    colouring      - one of two 4-letter strings, 'mvel' or 'temp', specifying colouring of cells
%}
 
    %% Define the fluid domain
    bx = [-L,-H/2,L,H/2];
    
    %% maximum and minimum vertical seed coordinates
    ymin = min(Z(:,2,:),[],'all');
    ymax = max(Z(:,2,:),[],'all');
    
    %% Make steps a row vector so that it can be used for indexing
    steps = steps(:)';
    
    num_steps = length(steps(:));
    
    %% Make figure
    if prod(colouring == 'mvel') == 1
         %% Colour by meridional velocity
        
        % Set colours of cells and variables defining the colour bar
        % NB: this invloves computing Laguerre tessellations so we save the
        % data to pass to the plotting function later
        colour_label          = 'Meridional velocity $\mathrm{ms}^{-1}$';
        meridional_velocities = zeros(size(Z,1),num_steps);
        
        for k=steps
            Zk               = Z(:,:,k);
            wk               = w(:,k);
            [~,~,xck,~]   = mexPDall_2d(bx,Zk,wk,perL,perV);
            mer_vels_k       = f*(Zk(:,1)-xck(:,1));
        
            meridional_velocities(:,k) = mer_vels_k;
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
   
            colours{3} = meridional_velocities(:,k);
            
            clf
            
            seedPlot = tiledlayout(1,1,'Padding','none');
            seedPlot.Units = 'inches';
            seedPlot.OuterPosition = [0.25 0.25 2.5 2];
            nexttile;
            
            Zk = Z(:,:,k);
            scatter(Zk(:,1),Zk(:,2),3.2,colours{3},'filled')
            
            xlim([bx(1),bx(3)]);
            ylim([ymin,ymax]);
            
            set(gca,'FontSize',9,'TickLabelInterpreter','latex')
            
            % set up colourbar with ticks at intervals specified by colours{2}(3) 
            % and spanning range specified by colours{2}(1:2)    
            denom = colours{2}(3);
            ticks = denom*(ceil(colours{2}(1)/denom):floor(colours{2}(2)/denom));
            caxis([ticks(1),ticks(end)]);

            cbr = colorbar('Ticks',ticks);
            cbr.Label.Interpreter = 'latex';
            cbr.Label.FontSize = 10;
            cbr.Label.String = colours{1};
            cbr.TickLabelInterpreter = 'latex';
        
            %saveas(gcf,strcat('seedsMerVelStep',num2str(k),'.fig'));
            %exportgraphics(gcf,strcat('seedsMerVelStep',num2str(k),'.png'),'Resolution',300);
            
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
            Zk          = Z(:,:,k);
            colours{3}  = c*Zk(:,2);
            
            clf
            
            cellPlot = tiledlayout(1,1,'Padding','none');
            cellPlot.Units = 'inches';
            cellPlot.OuterPosition = [0.25 0.25 2.5 2];
            nexttile;

            scatter(Zk(:,1),Zk(:,2),3.2,colours{3},'filled')
            
            xlim([bx(1),bx(3)]);
            ylim([ymin,ymax]);
            
            % set up colourbar with ticks at intervals specified by colours{2}(3) 
            % and spanning range specified by colours{2}(1:2)    
            denom = colours{2}(3);
            ticks = denom*(ceil(colours{2}(1)/denom):floor(colours{2}(2)/denom));
            caxis([ticks(1),ticks(end)]);

            cbr = colorbar('Ticks',ticks);
            cbr.Label.Interpreter = 'latex';
            cbr.Label.String = colours{1};
            cbr.TickLabelInterpreter = 'latex';
            
            %saveas(gcf,strcat('tempStep',num2str(k),'.fig'));
            %exportgraphics(gcf,strcat('tempStep',num2str(k),'.png'),'Resolution',300);
        end
    else
        error('Unrecognised 4-digit colour code. Use strings "temp" or "mvel" to specify colouring.')
    end
   
end