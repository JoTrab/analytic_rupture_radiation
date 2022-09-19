if isequal(Params.observation_plane,'x-z')
    %no division by zero
    r_0 = r; r_0(Source.Point(1), Source.Point(2)) = 0.5*r_0(Source.Point(1), Source.Point(2)+1);     
else
    %r_0 = r; r_0(Source.Point(1), Source.Point(2)) = 0.5*r_0(Source.Point(1), Source.Point(2)+1);     
    r_0 = r ;
end

xi = (depth_mat)-depth_mat(Source.Point(1,1),Source.Point(1,2));
x_distance =  r(Source.Point(1),:) ; x_distance = repmat(x_distance,size(r,1),1) ;
x_distance(:,1:Source.Point(2)-1) = -1* x_distance(:,1:Source.Point(2)-1) ;

% th = atan(xi./sqrt(x_distance.^2 + 0.^2)) ;     %elevation in matlab notation
%elevation in shearer/igel/... notation ---> acos(xi./r_0)



%careful with r for 3d case and for 2d case.
if isequal(Params.observation_plane,'x-z')
  r_0(Source.Point(1), Source.Point(2)) = nan  ;  %avoid complex results    
  xi_r = xi./r_0;
  xi_r(xi_r > 1) = 1 ;
  th = acos(xi_r);
  phi = atan2(0,x_distance);
  %zeros(size(th));  %atan(0./x_distance) ;  aber x_distance contains zeros
  %atan2 returns 0 for 0,0 atan nan for 0/0    %+ azimuth---> phi is 0 in 2d plane
  %note atan2(x,y) should be the same as atan(x/y) but taking into account
  %quadrants ---> phi should be either pi or zero
elseif isequal(Params.observation_plane,'x-y')
  xi_r = Source.Point(3)./r_0;
%   xi_r(xi_r > 1) = 1 ; %NOT NEEDED since no division by zero?
  th = acos(xi_r); %winkel zur vertikalen, haengt nur von r ab
  phi = atan2(abs(depth_mat),x_distance); %zeros(size(th));  %atan(0./x_distance) ;  aber
else
  error('Observation plane not defined, define observation plane as x-z or x-y' )  
%   error('us 2 entry source point for z-x simulation with y plane = 0 or 3 entry source point for y-x simulation with z plane = height in meter (third entry) ')
end






switch(wave_type)
    case{'s-wave' }
        
        
        if isequal(observation_dir,'vertical')
            [~,~,G_fct_s] = Green_function(beta,Source.X0,timesteps,r,phi,th,density,wave_type,[]);
        elseif isequal(observation_dir,'horizontal')
            [G_fct_s,~,~] = Green_function(beta,Source.X0,timesteps,r,phi,th,density,wave_type,[]);
        elseif isequal(observation_dir,'side')
            [~,G_fct_s,~] = Green_function(beta,Source.X0,timesteps,r,phi,th,density,wave_type,[]);
        end
        %%convolve
        
    case{ 'p-wave' }
        
        %P-field
        
        if isequal(observation_dir,'vertical')
            [~,~,G_fct_p] = Green_function(p_speed,Source.X1,timesteps,r,phi,th,density,wave_type,[]);
        elseif isequal(observation_dir,'horizontal')
            [G_fct_p,~,~] = Green_function(p_speed,Source.X1,timesteps,r,phi,th,density,wave_type,[]);
        elseif isequal(observation_dir,'side')
            [~,G_fct_p,~] = Green_function(p_speed,Source.X1,timesteps,r,phi,th,density,wave_type,[]);          
        end
        
    case{ 'inter_s' }
        if isequal(observation_dir,'vertical')
            [~,~,G_fct_I_s] = Green_function(beta,Source.X0,timesteps,r,phi,th,density,wave_type,[]);
        elseif isequal(observation_dir,'horizontal')
            [G_fct_I_s,~,~] = Green_function(beta,Source.X0,timesteps,r,phi,th,density,wave_type,[]);
        elseif isequal(observation_dir,'side')
            [~,G_fct_I_s,~] = Green_function(beta,Source.X0,timesteps,r,phi,th,density,wave_type,[]);  
        end
        
    case{ 'inter_p' }
        if isequal(observation_dir,'vertical')
            [~,~,G_fct_I_p] = Green_function(p_speed,Source.X1,timesteps,r,phi,th,density,wave_type,[]);
        elseif isequal(observation_dir,'horizontal')
            [G_fct_I_p,~,~] = Green_function(p_speed,Source.X1,timesteps,r,phi,th,density,wave_type,[]);
        elseif isequal(observation_dir,'side')
            [~,G_fct_I_p,~] = Green_function(p_speed,Source.X1,timesteps,r,phi,th,density,wave_type,[]);     
        end
        
    case{ 'near-field' }
        
        
        %         Source.X2 = zeros(size(r,1), size(r,2), length(time));
        %         for t = 1:length(time)-1
        %             X2 = zeros(size(r,1), size(r,2)) ;
        %             testplus_s = zeros(size(r,1), size(r,2));
        %             testplus_p = zeros(size(r,1), size(r,2));
        %
        %             testplus_s(:,:) = time(t).*ones(size(r)) - r(:,:)./beta;
        %             testplus_p(:,:) = time(t).*ones(size(r)) - r(:,:)./p_speed;
        %             Signp = sign(testplus_p) ;
        %             [B] = find(Signp == 1);
        %             Signs = sign(testplus_s) ;
        %             [A] = find(Signs == -1) ;
        %
        %             [C] = intersect(A,B) ;
        %             X2(C) = 1 ;
        %             Source.X2(:,:,t) = X2 ;
        %         end
        
        if isequal(observation_dir,'vertical')
            [~,~,G_fct_near] = Green_function(beta,Source.X2,timesteps,r,phi,th,density,wave_type,time,Params);
        elseif isequal(observation_dir,'horizontal')
            [G_fct_near,~,~] = Green_function(beta,Source.X2,timesteps,r,phi,th,density,wave_type,time,Params);
        elseif isequal(observation_dir,'side')
            [~,G_fct_near,~] = Green_function(beta,Source.X2,timesteps,r,phi,th,density,wave_type,time,Params);                 
        end
    case{ 'static' }
        
        if isequal(observation_dir,'vertical')
            [~,~,G_fct_static] = Green_function(beta,Source.X2,timesteps,r,phi,th,density,wave_type,time,Params);
        elseif isequal(observation_dir,'horizontal')
            [G_fct_static,~,~] = Green_function(beta,Source.X2,timesteps,r,phi,th,density,wave_type,time,Params);
         elseif isequal(observation_dir,'side')
            [~,G_fct_static,~] = Green_function(beta,Source.X2,timesteps,r,phi,th,density,wave_type,time,Params);                  
        end
        
end
