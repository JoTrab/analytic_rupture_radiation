%G is the Amplitude of the Green vector at that point. th and ph the
%coordinates neede for the unit vector conversion
    function [x,y,z] = sph2cart_vecProj(Gr, Gth, Gphi, th, phi)
        %
        % Transform spherical coordinates unit vectors to cartesian components
        % meaning project the spherical components to cartesian components

        %initialize matrixes of the right size (dim 3 because each grid
        %point is a vector)

        sin_phi = sin(phi) ;
        sin_phi(-0.00000002 > sin_phi < 0.00000002) = 0;
        r_unit = zeros([size(th),3]) ;th_unit = zeros(size(r_unit)); phi_unit = zeros(size(r_unit));
        r_unit(:,:,1) = sin(th) .* cos(phi);
        r_unit(:,:,2) =  sin(th) .* sin_phi;
        r_unit(:,:,3) =   cos(th);
        th_unit(:,:,1) = cos(th) .* cos(phi);
        th_unit(:,:,2) =  cos(th) .* sin_phi;
        th_unit(:,:,3) =   -sin(th);
        phi_unit(:,:,1) = -sin_phi;
        phi_unit(:,:,2) =  cos(phi);
        phi_unit(:,:,3) =    zeros(size(phi));

        xyz_r = bsxfun(@times,Gr, r_unit) ;
        xyz_th = bsxfun(@times,Gth,th_unit) ;
        xyz_phi =  bsxfun(@times,Gphi,phi_unit) ;

        x = xyz_r(:,:,1) + xyz_th(:,:,1) + xyz_phi(:,:,1) ;
        y = xyz_r(:,:,2) + xyz_th(:,:,2) + xyz_phi(:,:,2)  ;
        z = xyz_r(:,:,3) + xyz_th(:,:,3) + xyz_phi(:,:,3) ;
    end