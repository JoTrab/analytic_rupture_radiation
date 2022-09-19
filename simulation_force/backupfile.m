%make a nice source structur for several sources. the same holds  for X0
%and Source ind and so on. Then run each on overwriting the old ones but
%saving the result in wave so that you can add them all up in the end. 





userInput.maxDepth              = 145;   % Maximum depth of image, measured in wavelengths (typ. 192)
userInput.PRF                   = 3000;  % The ultrafast frame rate (pulse repetition frequency) in Hz 
userInput.NbrFrame              = 3000;    % Number of frames
userInput.NbrAngles             = 1;     % Number of angles for compound imaging
userInput.timeBetweenCompound   = 160;   % Time (µs) between frames for compound imaging
userInput.TrigDelay             = 1; 

[f,func_as,func_ps] = spect(dple(:,:,:),userInput.PRF ,'Dim',3);
func_m = mean(mean(func_as));
figure
plot(f,squeeze(func_m))
figure
plot(squeeze(dple(50,100,:)))
hold on
plot(squeeze(dple(50,50,:)),'o')
plot(squeeze(dple(50,10,:)),'r')


[f,func_as,func_ps] = spect(dpl(:,:,:),3000,'Dim',3);func_m = mean(mean(func_as));
hold on; figure;plot(f,squeeze(func_m))


%% grid for simulation
dpl = zeros(128,139,3000); dple = dpl; %initialize dpl if no data for comparison

Trans.frequency = 5.208; % nominal frequency in MHz
         % Vantage:  5.208 is closest supported frequency to 5 MHz
Trans.Bandwidth = [4, 7]
Trans.numelements = 128;
Trans.elementWidth = .250; % width in mm
Trans.spacingMm = .298;   % Spacing between elements in mm.
Trans.ElementPos = zeros(Trans.numelements,4);
Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));


% grid = zeros(size(dpl,2),size(dpl,1)); 
depth_meter = (0.*(1:size(dpl,2)))' ;
% put the origin in upperleft corner
x_vec = ((Trans.ElementPos(:,1)-Trans.ElementPos(1,1)) ./ 1000)';
x_mat = repmat(x_vec,[length(depth_meter), 1]);
depth_meter(:) = ((1:length(depth_meter)).*(1540/(Trans.frequency*1e6)));  %    WRONG!!! due to stelus speckle tracking?
depth_mat = repmat(depth_meter,[1 length(x_vec)]);
%source point



%% Parameters
addpath('/home/johnny/work/functions/my_functions/simulation')
beta = 4;
density = 1000
 %only use positive non-values of the source as source vector
observation_dir = 'horizontal' ;   %'vertical' 'horizontal'
timesteps = 500 ;  %timesteps to give in calculation
% time = kgrid.t_array ;%0:(1/3000):1;
total_timesteps = 3000 ;
t_end = 0.2 ;  
delta_t  = t_end/total_timesteps ; %time of one timestep in seconds

time = 0:(delta_t):t_end;

%time sampling equal to rupture time of 1 point to the next:

rup_speed = 10 ;
rup_time = Trans.spacingMm/1000 / rup_speed ; % time needed for rupture to advance one x-point in m
delta_t = rup_time ; 
time= 0:delta_t:t_end;
total_timesteps = length(time)

%%  structures for different sources - this is the same in the efficient script
nr_sources = 10
source_ind = [139 50]
source_cycle = 1
source_freq = 250
source_mag = 1
bw = 0.6 % bandwidth of gausspuls

%ricker
% [rw,trick] = ricker(source_freq,total_timesteps,delta_t,0);
%         source_ux = rw(1:find(abs(rw) >= 0.0001 ,1,'last'));
        
        
% %  gaussian shaped sinusoidal source 
tc = gauspuls('cutoff',source_freq,bw,[],-40);  %last truncation of db , 60 % bandwidth
t = -tc : delta_t : tc;  %sample the pulse at 3 khz
% time = [-fliplr(time(1:round(tc/delta_t))) time]  % add the negative timesteps 
% source_ux = gauspuls(time(1:2*round(tc/delta_t)),source_freq,0.6);
source_ux = gauspuls(t,source_freq,0.6);
plot(t,source_ux)
%         
        
        
        
        
        

r = repmat(struct('source_distance',zeros(length(depth_meter), length(x_vec))),1,nr_sources);


Source = repmat(struct('Point',source_ind, ...
                       'Cycles', source_cycle, ...
                       'Frequency_Hz',source_freq, ...
                       'Magnitude', source_mag, ...
                       'Function', source_ux),1,nr_sources);


for i = 1:nr_sources
    Source(i).Point(1,2) = source_ind(2) +(i-1)
    r(i).source_distance(:,:)=sqrt( ( x_mat(:,:) - x_mat(Source(i).Point(1),Source(i).Point(2) )).^2 + ( abs(depth_mat(:,:) - depth_mat(Source(i).Point(1),Source(i).Point(2) ))).^2)  ;
    source_ux = [0, source_ux]
    
    % this  is for a mixed phase wavelet:
%     [rw,trick] = ricker(source_freq,total_timesteps,delta_t,(i-1)*delta_t);
%         source_ux = rw(1:find(abs(rw) >= 0.0001 ,1,'last'));
if i == nr_sources-1
    break
end
    Source(i+1).Function = source_ux ;
end


%get r from random source point
% imagesc(r(30).source_distance)



% X0(Source_ind(1), Source_ind(2), 1) = 100;
  

%% assign the source function dependance

X0 = source_function(time, r(i).source_distance, beta, Source(i).Function)  ; 
wave = GSourceConvolution(observation_dir,beta,X0,timesteps,r(i).source_distance,Source(i).Point,depth_mat,density);

%% several sources
clear dpl dple
% wave_str = repmat( struct( 'wave', zeros([size(r(1).source_distance) timesteps])),1, nr_sources)
wave_tot = zeros([size(r(1).source_distance) timesteps]);
%     h = waitbar(0,'Please wait...')
for i = 2:nr_sources

X0 = source_function(time(1:timesteps), r(i).source_distance, beta, Source(i).Function)  ; 
wave = GSourceConvolution(observation_dir,beta,X0,timesteps,r(i).source_distance,Source(i).Point,depth_mat,density);
wave_tot = wave_tot+wave;
% wave_str(i).wave = wave;
%               waitbar(i/(nr_sources),h)
i
end
% close(h)

%delay the source function


%%

%sinus source
% 
% source.ux    =  -source_mag*sin(2*pi*source_freq*time( ... 
%                     1:floor( (length(time) / t_end) * source_cycles*(1/source_freq)))); 
%    impulssource  
% source.ux = zeros(size(time));
% source.ux(1) = 1

% %  gaussian shaped sinusoidal source 
% % yi = gauspuls(t,fc,bw)
% tc = gauspuls('cutoff',source_freq,0.6,[],-40);
% t = -tc : 1/3000 : tc;
% yi = gauspuls(t,source_freq,0.6);
% plot(t,yi)
% 

% %ricker
[rw,trick] = ricker(250,3000,1/3000,0);
source.ux = rw(1:find(abs(rw) >= 0.0001 ,1,'last'));



source_input = source.ux ;
wave_Str = {}
% for i = 1:3
X0 = source_function(time, r, beta, source.ux)  ; 
wave = GSourceConvolution(observation_dir,beta,X0,timesteps,r,Source_ind,depth_mat,density);
% wave_Str{1,i} = wave
% end






% rupture speed (take 10 m/s)
rup_speed = 10 ;
rup_time = Trans.spacingMm/1000 / rup_speed  % time needed for rupture to advance one x-point
%following sources
% source.ux(2:size(dpl,1)-Source_ind(2)) = 
% source.ux(2) =
% source.ux(3) =
Source_ind = Source_ind'
Source_ind(2,:) = [Source_ind(1,1) 62]
[rw,t] = ricker(250,3000,1/3000,12);    %delay time comes into play here
source.ux(2,:) = rw(1:find(abs(rw) >= 0.0001 ,1,'last'));
Source_ind(3,:) = [Source_ind(1,1) 74]
[rw2,t] = ricker(250,3000,1/3000,24);    %delay time comes into play here
source.ux(3,:) = rw(1:find(abs(rw) >= 0.0001 ,1,'last'));
source.ux = source_input(i,:);















%% simple source


for t = 1:length(time)

        X_test = zeros(size(r,1), size(r,2));
        
            test(:,:,t) =  time(t).*ones(size(r)) - r(:,:)./beta;
            testplus = test(:,:,t);
      
            if t > 1
                testminus = test(:,:,t-1);
                
                
                Sign_1 = sign(testminus) ;
                Sign1 = sign(testplus) ;
                [A] = find(Sign_1 == -1);
                [B] = find(Sign1 == 1);
                [C] = find(Sign_1 == 0);
                [IND] =  intersect(A,B);
                [I,J] = ind2sub([139,128],IND)
                
                X_test(IND) = 100;
                
                if isempty(C) == 0
                    [I,J] = ind2sub([139,128],C)
                    X_test(C) = 100 ;
                end
                %         imagesc(squeeze(X0(:,:,t-1)))
                %         pause
              X0(:,:,t) = X_test;
              imagesc(squeeze(X0(:,:,t-1)))
                        pause
            end
            
      
end



    
%% k wave make grid  ---- length is right but position is shifted by half an elemnt width compared to verasonics acquisition
Trans.frequency = 5.208; % nominal frequency in MHz
         % Vantage:  5.208 is closest supported frequency to 5 MHz
Trans.Bandwidth = [4, 7]
Trans.numelements = 128;
Trans.elementWidth = .250; % width in mm
Trans.spacingMm = .298;   % Spacing between elements in mm.
Trans.ElementPos = zeros(Trans.numelements,4);
Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));


samplesPerWavelength = 4;
Nx = Trans.numelements; dx =  Trans.spacingMm / 1000; Ny = (size(dpl,2))  ; dy = (1540/(Trans.frequency*1e6)) / samplesPerWavelength;
Nz = 30; dz = dx;
% create the computational grid
PML_size = 10;

kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);


% define the properties of the upper layer of the propagation medium
medium.sound_speed_compression = 1500*ones(Nx, Ny, Nz); % [m/s]
medium.sound_speed_shear       = 4*ones(Nx, Ny, Nz);     % [m/s]
medium.density                 = 1000*ones(Nx, Ny, Nz); % [kg/m^3]

% define the properties of the lower layer of the propagation medium
% medium.sound_speed_compression(Nx/2:end, :, :) = 2000;  % [m/s]
% medium.sound_speed_shear(Nx/2:end, :, :)       = 800;   % [m/s]
% medium.density(Nx/2:end, :, :)                 = 1200;  % [kg/m^3]

% define the absorption properties
medium.alpha_coeff_compression = 0.1; % [dB/(MHz^2 cm)]
medium.alpha_coeff_shear       = 0.0; % [dB/(MHz^2 cm)]

% create the time array
cfl = 0.5;  % courant friedrichi criterion
t_end = 30/3000;
kgrid.t_array = makeTime(kgrid, max(medium.sound_speed_shear(:)), cfl, t_end);


source_freq = 250;      % [Hz]
source_cycles = 1;
source_mag = 1e-6;

% define source mask to be a square piston
source_x_pos = 50;      % [grid points]
% source_radius = 2;     % [grid points]
source.u_mask = zeros(Nx, Ny, Nz);
% sinus source function adapted to the time of acquisition and the t_array
source.ux    =  -source_mag*sin(2*pi*source_freq*kgrid.t_array( ... 
                    1:floor( (length(kgrid.t_array) / t_end) * source_cycles*(1/source_freq))));     
source.uy     =   0   
source.uz     =    0   
% source.u_mask(source_x_pos, Ny/2 - source_radius + 1:Ny/2 + source_radius, Nz/2 - source_radius + 1:Nz/2 + source_radius) = 1;
source.u_mask(Nx/4, Ny(1):Ny(end) , Nz/2) = 1;
% source.p_mask(Nx/4, Ny/2 - source_radius:Ny/2 + source_radius, Nz/2 - source_radius:Nz/2 + source_radius) = 1;
% define source to be a velocity source

% source.ux = source_mag*toneBurst(1/kgrid.dt, source_freq, source_cycles);

% set source focus
% source.ux = focus(kgrid, source.ux, source.u_mask, [0, 0, 0], 1500);

% % define a single source point
% source.u_mask = zeros(Nx, Ny);
% source.u_mask(end - Nx/4, Ny/2) = 1;
% 
% % define a time varying sinusoidal velocity source in the x-direction
% source_freq = 0.25e6;
% source_mag = 2/(medium.sound_speed*medium.density);
% source.ux = -source_mag*sin(2*pi*source_freq*kgrid.t_array);
% 
% % filter the source to remove high frequencies not supported by the grid
% source.ux = filterTimeSeries(kgrid, medium, source.ux);

% define sensor mask in x-y plane using cuboid corners, where a rectangular
% mask is defined using the xyz coordinates of two opposing corners in the
% form [x1, y1, z1, x2, y2, z2].'

rect1_x_start = 15;
rect1_y_start = 15;
rect1_x_end   = 118;
rect1_y_end   = 129 ;
rect1_z_start   = 14;
rect1_z_end   = 15;
sensor.mask = [rect1_x_start, rect1_y_start, rect1_z_start, rect1_x_end, rect1_y_end, rect1_z_end].';
% sensor.mask = [1 + PML_size, 1 + PML_size, Nz/2, Nx - PML_size, Ny - PML_size, Nz/2].';

% record the maximum pressure in the plane
sensor.record = {'p_max','u_final'};

% define input arguments
input_args = {'PlotScale', [-2, 2, -0.1, 0.1], 'DataCast', 'single',...
    'PMLSize', PML_size, 'DisplayMask', source.u_mask};


% run the simulation with PML inside
sensor_data = pstdElastic3D(kgrid, medium, source, sensor, input_args{:});

%%

% define a series of Cartesian points to collect the data
x = repmat(kgrid.x_vec',[1 139]);          % [m]
y = repmat(ones(size(kgrid.x_vec')),[1 139]);     % [m]
for i = 1:length(kgrid.y_vec); y(i*length(kgrid.x_vec)-(length(kgrid.x_vec)-1):i*length(kgrid.x_vec)) = kgrid.y_vec(i); end
z =ones(size(x))*0          % [m]
sensor.mask = [x; y; z];

% input arguments
input_args = {'PlotLayout', true, 'PlotPML', false, ...
    'DataCast', 'single', 'CartInterp', 'nearest'};

% run the simulation
% sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
%%
% =========================================================================
% VISUALISATION
% =========================================================================

% % plot the simulation layout using voxelplot
% voxelPlot(double(source.p0 | cart2grid(kgrid, sensor.mask)));
% view([50, 20]);

% plot the simulated sensor data
figure;
imagesc(sensor_data, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;





for t = 1:50
    SOL(:,:,t) = (1./ (4*pi.*r(:,:))).*1000.*(((t-1)/3000)-(r./(8)));
    imagesc(SOL(:,:,t))
    pause
end

1*exp(i*(2*pi*Trans.frequency*1e6).*r(:,:) );








































%% The FUNctions to be called 

%% complicated source

 h = waitbar(0,'Please wait...');

for t = 1:length(time)-1
    if t <= length(source.ux)-1
        X_test = zeros(size(r,1), size(r,2),t);
        
        for s_iter = 1:t
        X_iter =   zeros(size(r,1), size(r,2)); 
            
            
            
            
            test(:,:,t) =  time(t-(s_iter-1)).*ones(size(r)) - r(:,:)./beta;
            testplus = test(:,:,t);
            %     imagesc(test(:,:,t))
            %     colorbar
            %     pause
            if t > 1 & s_iter  < t
                testminus = test(:,:,t-(s_iter-1)-1);
                
                
                Sign_1 = sign(testminus) ;
                Sign1 = sign(testplus) ;
                [A] = find(Sign_1 == -1);
                [B] = find(Sign1 == 1);
                [C] = find(Sign_1 == 0);
                [IND] =  intersect(A,B);
                X_iter(IND) = source.ux(s_iter);
                
                if isempty(C) == 0
                    X_iter(C) = source.ux(s_iter);
                end
                %         imagesc(squeeze(X0(:,:,t-1)))
                %         pause
            end
            X_test(:,:,s_iter) = X_iter;
        end
        
  
%     end
% end
    
    else
         X_test = zeros(size(r,1), size(r,2),length(source.ux)-1);
        
        for s_iter = 1:length(source.ux)-1;
           
            
             
        X_iter =   zeros(size(r,1), size(r,2)); 
             
            test(:,:,t) =  time(t-(s_iter-1)).*ones(size(r)) - r(:,:)./beta;
            testplus = test(:,:,t);
            %     imagesc(test(:,:,t))
            %     colorbar
            %     pause
            if t > 1 & s_iter  < t
                testminus = test(:,:,t-(s_iter-1)-1);
                
                
                Sign_1 = sign(testminus) ;
                Sign1 = sign(testplus) ;
                [A] = find(Sign_1 == -1);
                [B] = find(Sign1 == 1);
                [C] = find(Sign_1 == 0);
                [IND] =  intersect(A,B);
                X_iter(IND) = source.ux(s_iter);
                
                if isempty(C) == 0
                    X_iter(C) = source.ux(s_iter);
                end
                %         imagesc(squeeze(X0(:,:,t-1)))
                %         pause
            end
            X_test(:,:,s_iter) = X_iter;
        end
            
    end
              X0(:,:,t) = mean(X_test,3);
              waitbar(i/(length(time)-1),h)
%    imagesc(X0(:,:,t),[min(source.ux) max(source.ux)]); title(t);pause(0.000001);
end
% end


close(h)

    
    

    
 
%% include tne wave function  x/r and 1 (or zero) for directivity
%measure against the vertical  ith component of displacement at x t  by dirac force
%at x,t = 0 in direction xj
% obs_dir = 0 ; obs_dir = deg2rad(obs_dir) ;
% Source direction is always horizontal 
% observation can be chosen between vertical and horizontal 
observation_dir = 'vertical'    ; % 'horizontal' 'vertical'  obs_dir
beta  = 4;
% source_dir = 90 ; source_dir = deg2rad(source_dir) ;  %of the source against the horizontal
wave = zeros(size(X0,1),size(X0,2), 200);





r_0 = r; r_0(Source_ind(1), Source_ind(2)) = 0.5*r_0(Source_ind(1), Source_ind(2)+1);     %no division by zero
xi = flipud(depth_mat)-depth_mat(1,1);
x_distance = r(139,:) ; x_distance = repmat(x_distance,size(r,1),1) ;
x_distance(:,1:Source_ind(2)-1) = -1* x_distance(:,1:Source_ind(2)-1) ;
switch(observation_dir)
    case{'vertical'}
        r_0(:,1:Source_ind(2)-1) = -1* r_0(:,1:Source_ind(2)-1) ;

        yi = (real(xi(:,:)./r_0(:,:)));   %NEEDED FOR VALUES CLOSE TO ZERO  this is the cosine of the angel of S-R  to the vertical
        'vertical'
        obs_dir = deg2rad(0)
                source_dir = 90 ; source_dir = deg2rad(source_dir) ;

    case{'horizontal'}
                r_0(:,1:Source_ind(2)-1) = -1* r_0(:,1:Source_ind(2)-1) ;

        yi =  (real(x_distance(:,:)./r_0(:,:)));   %NEEDED FOR VALUES CLOSE TO ZERO  this is the cosine of the angel of S-R  to the hprizontal
%         yi = cos(deg2rad((90))) - abs(yi);
%                 yi = (real(xi(:,:)./r_0(:,:)));  
        source_dir = 0 ; source_dir = deg2rad(source_dir) ;
                obs_dir = deg2rad(0)
        'horizontal'
    otherwise
        yi = (real((cos(obs_dir).*x_distance(:,:) )./r_0(:,:)));
        yi = (real((cos(obs_dir).*xi(:,:) )./r_0(:,:)));
%         yi = (real((xi(:,:) )./r_0(:,:)));
%         source_dir = abs(obs_dir - source_dir) ;
        'angle'
end

yij = real(cos((source_dir - acos(yi)))); %cos of angle between force and source-receiver 
delta_angle = ((source_dir - obs_dir)) ;
% delta = zeros(size(r));
delta = 0
if floor(delta_angle) == 0   
%     IND = find(yij >= 0.999);  % find cos =1 therefore angle = 0
%     delta(IND) = 1;
%     imagesc(delta);
delta = 1
end
% [I, J] = ind2sub(size(yij),IND);
% I(find(J == 50,1, 'first'):end) J(J >= Source_ind(2)) 
% IND = sub2ind( size(delta) ,I(J >= Source_ind(2)),J(J >= Source_ind(2)));


% yi = abs(yi); yj = abs(yj) ; yij = abs(yij) ;

% if source_dir == 0 | source_dir == 180
%     delta = 1
% else
%     delta = 0
% end
for  t = 1:200
%     directivity

%     yj = real(cos(abs(source_dir - acos(yi))));
    
%     delta(:,:)
    
    wave(:,:,t) = 1./(4*pi*1000*beta^2) .*  (delta - yij.*yi)  .*  1./(r(:,:)) .*  X0(:,:,t)  ;%.* (1 - [0 20])
%     imagesc(wave(:,:,t));
%     pause
%     ytick = wave(:,:,t) 8     ;
end
