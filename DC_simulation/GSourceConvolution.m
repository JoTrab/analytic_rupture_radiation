function wave = GSourceConvolution(observation_dir,beta,X0,timesteps,r,Source_ind,depth_mat,density,wave_type)
% not needed anymore (????)
% SOURCE Green Convolution for far field s , gives the
% resulting field at observation point in time
% USAGE:
%       X0 = source_function(time, r, beta, source)  X0(t) is point force in source direction
%                                                    X0(t-r/beta)
% INPUTS:
%       observation_dir = 'vertical' or 'horizontal' for position of US
%       probe
%
%       timesteps   - timesteps to calculate [s]
%       r  - distance vector field from source (source - receiver vector) [m]
%       beta  - wave speed  [m/s]
%       Source_ind  - Source position in the receiver array (shifted in
%                      between the fiven and the next grid point in this function to avoid
%                      division by 0
%       depth_mat
%       density    - density of the homogenous material (1000 for gels) [kg/m³]







r_0 = r; r_0(Source_ind(1), Source_ind(2)) = 0.5*r_0(Source_ind(1), Source_ind(2)+1);     %no division by zero
xi = flipud(depth_mat)-depth_mat(1,1);
x_distance = r(size(r,1),:) ; x_distance = repmat(x_distance,size(r,1),1) ;
x_distance(:,1:Source_ind(2)-1) = -1* x_distance(:,1:Source_ind(2)-1) ;
yi = zeros(size(xi)) ;
yij = zeros(size(xi)) ;


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
wave = zeros([size(r) timesteps]);
for  t = 1:timesteps
    %     directivity
    
    %     yj = real(cos(abs(source_dir - acos(yi))));
    
    %     delta(:,:)
    switch(wave_type)
        case{'s-wave'}
        wave(:,:,t) = 1./(4*pi*density*beta^2) .*  (delta - yij.*yi)  .*  (1./r(:,:)) .*  X0(:,:,t)  ;%.* (1 - [0 20])
        case{'p-wave'}
            %wave(:,:,t) = 1./(4*pi*density*beta^2.*r(:,:)) .*  ((yij.*yi)./(r(:,:).^2))  .*  X0(:,:,t)  ;
	    wave(:,:,t) = 1./(4*pi*density*beta^2) .*  (yij.*yi) .* (1./r(:,:))  .*  X0(:,:,t)  ;   
        
        
        wave(:,:,t) = 1./(4*pi*density) .*  (3.* (yij.*yi) - delta)  .*  (1./(r(:,:).^3)) .* ...
            X0(:,:,t)  ;%.* (1 - [0 20])

        

 end
    
    %     imagesc(wave(:,:,t));
    %     pause
    %     ytick = wave(:,:,t) 8     ;
end
