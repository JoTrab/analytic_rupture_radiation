function wave = Green_function(observation_dir,observation_plane,beta,...
    X0,timesteps,r,Source_ind,depth_mat,density,wave_type,time,Source_Direction)
% source direction is given relative to observation direction 
% SOURCE Green function for far field s , gives the
% resulting field at observation point in time for a delta
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
source_dir = Source_Direction;
source_dir = deg2rad(source_dir) ;

if isequal(observation_plane,'x-z')
    %no division by zero
    r_0 = r; r_0(Source_ind(1), Source_ind(2)) = 0.5*r_0(Source_ind(1), Source_ind(2)+1);     
    r_0(Source_ind(1), Source_ind(2)) = nan;  
else
    r_0 = r ;
end
xi = (depth_mat)-depth_mat(Source_ind(1,1),Source_ind(1,2));
x_distance = r(Source_ind(1),:) ; x_distance = repmat(x_distance,size(r,1),1) ;
x_distance(:,1:Source_ind(2)-1) = -1* x_distance(:,1:Source_ind(2)-1) ;
yi = zeros(size(xi)) ;
yij = zeros(size(xi)) ;



switch(observation_plane)
    case{'x-y'}
        switch(observation_dir)
            case{'vertical'}
                yi = (real(Source_ind(3)./r_0)); %ankathete/hypo =cos
                %this is the cosine
                yi(yi > 1) = 1 ;
                yi(yi < -1) = -1 ;
                yi_angle =  acos(yi);
                %define if
%                 yi_angle(:,1:Source_ind(2)-1) = -yi_angle(:,1:Source_ind(2)-1) ; %angle to vertical   
%                 or
%                 yi_angle(Source_ind(1)+1:end,:) = -yi_angle(Source_ind(1)+1:end) ; %angle to vertical   
%                   depending on source. 
                error('define proper angles')
                'vertical'
                obs_dir = deg2rad(0) ;  %neasured to vertical  
            case{'horizontal'}
                yi =  (real(x_distance./r_0));    %xi?
                yi(yi > 1) = 1 ; 
                yi_angle =  acos(yi);  
                %given towards the axis!!
                yi_angle(Source_ind(1)+1:end,:) = -yi_angle(Source_ind(1)+1:end,:) ;
%                 yi_angle(:,1:Source_ind(2)) = -yi_angle(:,1:Source_ind(1)) ;

                obs_dir = deg2rad(0) ; %always 0 because source is relative to observation direction. 
            case{'side'}
                yi =  (real(xi./r_0));    %xi?
                yi(yi > 1) = 1 ; 
                yi_angle =  acos(yi);        
                yi_angle(:,1:Source_ind(2)-1) = -yi_angle(:,1:Source_ind(2)-1) ; %angle to vertical   
                
            otherwise
                error('define for angles') ;
        end
    case{'x-z'}
        switch(observation_dir)
            case{'vertical'}
                %this is the cosine
                yi = (real(xi./r_0));   
                %NEEDED FOR VALUES CLOSE TO ZERO  this is the cosine of the angel of S-R  to the vertical
                yi(yi > 1) = 1 ;
                yi_angle =  acos(yi);
                %set negative: neede for the force receiver angle (to get full 180
                %degree in positive and negative direction)
                yi_angle(:,1:Source_ind(2)-1) = -yi_angle(:,1:Source_ind(2)-1) ; %angle to vertical   
                'vertical'
                obs_dir = deg2rad(0) ;  %neasured to vertical  
                %both angles are to vertical
        %         yij = real(cos((source_dir - (yi_angle)))); %yij (yj in aki is the cos of angle between force direction and source-receiver angle (angle to observation direction)
            case{'horizontal'}
                yi =  (real(x_distance./r_0));   
                %NEEDED FOR VALUES CLOSE TO ZERO  this is the cosine of the angel of S-R  to the horizontal
                yi(yi > 1) = 1 ; 
                yi_angle =  acos(yi);        
                yi_angle(Source_ind(1)+1:end,:) = -yi_angle(Source_ind(1)+1:end,:) ;
%                 yi_angle(1:Source_ind(1),:) = -yi_angle(1:Source_ind(1),:) ;

                obs_dir = deg2rad(0) ; %always 0 because source is relative to observation direction. 
            otherwise
                %to do
                %create a distance vector in the observation direction
        %         yi = (real((cos(obs_dir).*x_distance(:,:) )./r_0(:,:)));
        %         yi = (real((cos(obs_dir).*xi(:,:) )./r_0(:,:)));
                %         yi = (real((xi(:,:) )./r_0(:,:)));
                %         source_dir = abs(obs_dir - source_dir) ;
                'angle' ;
        end
end

%if the angle is larger 90 source direction is opposite to observatio
%direction

yij = real(cos((source_dir - (yi_angle)))); %yij 
%(yj in aki is the cos of angle between force direction and uP 
%(for the P wave this is source-receiver angle)

%Sourrce direction is mainly needed for the kronecker delta
%test if source and observation are in one line
delta_angle = ((abs(cos(source_dir)) - abs(cos(obs_dir)))) ;
% delta = zeros(size(r));
kroneck = 0 ;
if abs(0-(delta_angle)) <= 0.001
    %     IND = find(yij >= 0.999);  % find cos =1 therefore angle = 0
    %     delta(IND) = 1;
    %     imagesc(delta);
    kroneck = 1 ;
    %kroneck is needed for the s-wave to bring in orthogonality to P
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
wave = single(zeros([size(r) timesteps]));
for  t = 1:timesteps
    %     directivity

    %     yj = real(cos(abs(Del_dir - acos(yi))));

    %     delta(:,:)
    switch(wave_type)
        case{'s-wave'}  %the kroneck in the s-wave brnigs in the orthogonality  og yij compared to  p -wave
            wave(:,:,t) = -1 .* (1/(4*pi*density*(beta^2))) .*  (yij.*yi - kroneck)  .*  (1./r) .*  X0(:,:,t)  ;%.* (1 - [0 20])
        case{'p-wave'}
            %wave(:,:,t) = 1./(4*pi*density*beta^2.*r(:,:)) .*  ((yij.*yi)./(r(:,:).^2))  .*  X0(:,:,t)  ;
            wave(:,:,t) = (1/(4*pi*density*beta^2)) .*  (yij.*yi) .* (1./r)  .*  X0(:,:,t)  ;

        case{'near-field'}
            wave(:,:,t) = (1/(4*pi*density)) .*  (3.*(yij.*yi) - kroneck)  .*  (1./(r.^3)) .* time(t) .* X0(:,:,t)  ;%.* (1 - [0 20])  time(t) is the tau



    end

    %     imagesc(wave(:,:,t));
    %     pause
    %     ytick = wave(:,:,t) 8     ;
end
