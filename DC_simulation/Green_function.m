function [Gx,Gy,Gz] = Green_function(beta,X0,timesteps,r,phi,th,...
    density,wave_type,time,Params)
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

% x_distance is x - axis
% xi is z axis
% y axis is 0
%----> theta phi und r definiert wie in igel/aki/shearer
test = [] ;


wave = single(zeros([size(r) 3 timesteps]));
Gr = single(zeros([size(r) timesteps]));
Gth = single(zeros([size(r) timesteps]));
Gphi = single(zeros([size(r) timesteps]));
Gx =single(zeros([size(r) timesteps]));
Gy = single(zeros([size(r) timesteps]));
Gz = single(zeros([size(r) timesteps]));

sin_phi = sin(phi) ;
sin_phi(-0.00000002 > sin_phi < 0.00000002) = 0;

switch(wave_type)
    case{'s-wave'}
        AFS = single(zeros([size(th), 3]));
        AFS(:,:,1) = 0; 
        AFS(:,:,2) = cos(2*th).*cos(phi); 
        AFS(:,:,3) = -cos(th).*sin_phi  ;
    case{'p-wave'}
        AFP = single(zeros([size(th), 3]));
        AFP(:,:,1) = sin(2*th).*cos(phi); 
        AFP(:,:,2) = 0; 
        AFP(:,:,3) = 0   ;

    case{'inter_s'}
        AIS = single(zeros([size(th), 3]));
        
        AIS(:,:,1) = -3*sin(2*th).*cos(phi); 
        AIS(:,:,2) = +3*cos(2*th).*cos(phi); % me
        AIS(:,:,3) = +1*cos(th).*sin_phi  ;
%         AIS(:,:,2) = -3*cos(2*th).*cos(phi); %aki
%         AIS(:,:,3) = +1*cos(th).*sin_phi  ;
    case{'inter_p'}
        AIP = single(zeros([size(th), 3]));
        
        AIP(:,:,1) = 4*sin(2*th).*cos(phi); 
        AIP(:,:,2) =  -2*cos(2*th).*cos(phi);
        AIP(:,:,3) = -cos(th) .* sin_phi  ;
    case{'near-field'}

        AN = single(zeros([size(th), 3]));
        AN(:,:,1) = (9*sin(2*th).*cos(phi)) ;
        AN(:,:,2) = -6*cos(2*th).*cos(phi) ;
        AN(:,:,3) = -1*cos(th).*sin_phi ;  %(aki + I minus)
        %test with 69 being x of the source
% wave_type = 'near-field' ;
% Green_type_switch ;
% imagesc(squeeze(G_fct_near(:,69,1:60)))  vergleichen in gleichem abstand
% nan, minimalst unterschied liegen eventuell in der diskretisierung
    case{'static'}
%         AIS = zeros([size(th), 3]);
%         AIS(:,:,1) = -3*sin(2*th).*cos(phi); AIS(:,:,2) = 3*cos(2*th).*cos(phi);
%         AIS(:,:,3) = -3*cos(th).*sin_phi  ;
%         AN = zeros([size(th), 3]);
%         AN(:,:,1) = (9*sin(2*th).*cos(phi)) ;
%         AN(:,:,2) = -6*cos(2*th).*cos(phi) ;
%         AN(:,:,3) = 6*cos(th).*sin_phi ;
%         AIP = zeros([size(th), 3]);
%         AIP(:,:,1) = 4*sin(2*th).*cos(phi); AIP(:,:,2) =  -2*cos(2*th).*cos(phi);
%         AIP(:,:,3) = 2*cos(th) .* sin_phi  ;

        wave = single(zeros([size(r) 3]));
%         alpha = evalin('base','Params.p_speed') ;
%         beta = evalin('base','beta') ;
        alpha = Params.p_speed% evalin('base','Params.p_speed') ;
        beta = Params.s_speed ;%Sevalin('base','beta') ;
        dt = abs(time(2) - time(1)) ; %AIS and AIO are scaled by dt (delta function INTEGRAL =1)
%         wave = bsxfun(@times,(  (1/(4*pi*density)) .* (1./r.^2) ) , ...
%            ((1/dt).*AN.*(1/(2*(beta^2)) - 1/(2*(alpha^2)))) + (1/dt) .* (AIP./(alpha^2)) +  (1/dt) .* (AIS./(beta^2) )) ;
%

       %i_s only
%              wave = bsxfun(@times,(  (1/(4*pi*density)) .* (1./(r.^2)) ) , ...
%           (1/dt) .* (AIS./(beta^2) )) ;
      %near only
%       wave = bsxfun(@times,(  (1/(4*pi*density)) .* (1./r.^2) ) , ...
%           ((1/dt).*AN.*(1/(2*(beta^2)) - 1/(2*(alpha^2))))  ) ;

%whole expression
                A_static = single(zeros([size(th),3]));
                A_static(:,:,1) = 0.5*(3/beta^2 -1/alpha^2) .* sin(2*th).*cos(phi) ;
%             (1/alpha^2 ).* cos(2*th).*cos(phi) , (-1/alpha^2) .* cos(th).*sin(phi)])
                A_static(:,:,2) =  (1/alpha^2 ).* cos(2*th).*cos(phi) ;
                 A_static(:,:,3) =   (-1/alpha^2) .* cos(th).*sin_phi;

  wave = bsxfun(@times, ((1/(4*pi*density)) .* (1./r.^2) ) , (1/dt).*A_static) ;

        Gx = zeros(size(r));
        Gy = zeros(size(r));
        Gz = zeros(size(r));
        [Gx,Gy,Gz] = sph2cart_vecProj(wave(:,:,1), wave(:,:,2), wave(:,:,3),th,phi)   ;

end


% for i = 1:3   %loop over the component ( theta/phi and r)
if isequal(wave_type, 'static') == 0
    for  t = 1:timesteps
        %     directivity

        %     yj = real(cos(abs(Del_dir - acos(yi))));

        %     delta(:,:)
        switch(wave_type)
            case{'s-wave'}   %bsxfun does elementwise product of matrix --> like outer product but for matrix
                wave(:,:,:,t) = bsxfun(@times,((1/(4*pi*density*(beta^3))) .*   (1./r)),  bsxfun(@times,AFS, X0(:,:,t))) ; %.* (1 - [0 20])
            case{'p-wave'}
                %wave(:,:,t) = 1./(4*pi*density*beta^2.*r(:,:)) .*  ((yij.*yi)./(r(:,:).^2))  .*  X0(:,:,t)  ;
                wave(:,:,:,t) = bsxfun(@times,((1/(4*pi*density*(beta^3))) .*  (1./r)), bsxfun(@times,AFP , X0(:,:,t)))  ;

            case{'inter_s'}
                wave(:,:,:,t) =   bsxfun(@times,((1/(4*pi*density*(beta^2))) .*  (1./(r.^2))), bsxfun(@times, AIS  , X0(:,:,t)));%.* (1 - [0 20])
            case{'inter_p'}
                wave(:,:,:,t) =   bsxfun(@times,((1/(4*pi*density*(beta^2))) .*  (1./(r.^2))) , bsxfun(@times,AIP , X0(:,:,t))) ;%.* (1 - [0 20])
            case{'near-field'}

                wave(:,:,:,t) = bsxfun(@times,((1/(4*pi*density)) .* (1./(r.^4))  .* time(t) ), bsxfun(@times,AN,   X0(:,:,t)))  ;%.* (1 - [0 20])  time(t) is the tau

        end

        %     imagesc(wave(:,:,t));
        %     pause
        %     ytick = wave(:,:,t) 8     ;
        % end

        Gr(:,:,t) = squeeze(wave(:,:,1,t)); % spherical componets of the field
        Gth(:,:,t) = squeeze(wave(:,:,2,t));
        Gphi(:,:,t) = squeeze(wave(:,:,3,t));
        [Gx(:,:,t),Gy(:,:,t),Gz(:,:,t)] = sph2cart_vecProj(Gr(:,:,t), Gth(:,:,t), Gphi(:,:,t),th,phi)   ;

        %     if t == 1000
        %         max(max(max(Gx)))
        %         max(max(max(Gy)))
        %         max(max(max(Gz)))
        %         max(max(max(Gr)))
        %         max(max(max(Gth)))
        %         max(max(max(phi)))
        %         pause
        %     end
    end
end

% NOT use [Gx, Gy, Gz] = sph2cart( Gphi ,Gth, Gr );    % spherical to cartesian coordinates

% %G is the Amplitude of the Green vector at that point. th and ph the
% %coordinates neede for the unit vector conversion
%     function [x,y,z] = sph2cart_vecProj(Gr, Gth, Gphi, th, phi)
%         %
%         % Transform spherical coordinates unit vectors to cartesian components
%         % meaning project the spherical components to cartesian components
%
%         %initialize matrixes of the right size (dim 3 because each grid
%         %point is a vector)
%
%         sin_phi = sin(phi) ;
%         sin_phi(-0.00000002 > sin_phi < 0.00000002) = 0;
%         r_unit = zeros([size(th),3]) ;th_unit = zeros(size(r_unit)); phi_unit = zeros(size(r_unit));
%         r_unit(:,:,1) = sin(th) .* cos(phi);
%         r_unit(:,:,2) =  sin(th) .* sin_phi;
%         r_unit(:,:,3) =   cos(th);
%         th_unit(:,:,1) = cos(th) .* cos(phi);
%         th_unit(:,:,2) =  cos(th) .* sin_phi;
%         th_unit(:,:,3) =   -sin(th);
%         phi_unit(:,:,1) = -sin_phi;
%         phi_unit(:,:,2) =  cos(phi);
%         phi_unit(:,:,3) =    zeros(size(phi));
%
%         xyz_r = bsxfun(@times,Gr, r_unit) ;
%         xyz_th = bsxfun(@times,Gth,th_unit) ;
%         xyz_phi =  bsxfun(@times,Gphi,phi_unit) ;
%
%         x = xyz_r(:,:,1) + xyz_th(:,:,1) + xyz_phi(:,:,1) ;
%         y = xyz_r(:,:,2) + xyz_th(:,:,2) + xyz_phi(:,:,2)  ;
%         z = xyz_r(:,:,3) + xyz_th(:,:,3) + xyz_phi(:,:,3) ;
%     end

end
