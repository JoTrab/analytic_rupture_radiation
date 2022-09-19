clear all
% This code reproduce 

% new x axis: Trans.numelements_true*Trans.spacingMm - Trans.spacingMm*128
% 22.6480 : 38.1440

addpath(genpath('analytic_rupture_simulation/simulation_force/'))
set(0,'DefaultFigureWindowStyle','Docked')
Params.Vp_Dpl = 'Dpl'   % decide if you want displacement or Vp as primary output
% Trans.numelements_true = 202; %90 to have a third 'outside imaging region'
Trans.numelements_true = 204; %90 to have a third 'outside imaging region'

Trans.size_factor = 10; %10 reduce size to save memory and increase distance in z, memory saving
Trans.Depthpoints = 1680/Trans.size_factor;  %delta depth is 1/8 of ...
        % wavelength ultrasound scaled by size factor
Trans.frequency = 5.208; % nominal frequency in MHz
Trans.Bandwidth = [4, 7] ;
Trans.elementWidth = .250; % width in mms
Trans.spacingMm = .298;   % Spacing between elements in mm.
Params.Fs = 3000;
% Params.returnSwitch = 'Field' ; % 'Green', 'Field' ,'none'
Params.Field_only = 'off'
Params.s_speed = 6.0;   %1 m/s error well possible at 1/3000 PRF and the measured distance of 14.40 mm



Rupture.rup_speed_start_factor =   4 %4.0462 %3.8116%2.7  %sqrt(2)+0.1  %1.0 ;% (0.01/0.0497)*52.5000 +((0.027+(10.9*0.001))/0.0497)*10.6
Rupture.rup_speed_end_factor =  1.3%10/6.9    %sqrt(2)+0.1    % 0.1  ;%sqrt(2)  % sqrt(2) ;  %sqrt(2)
%Give a maximum allowed rupture time
Rupture.max_rup_time = 0.03  ;  %total time from one end of rupture to the other
%Give a minimum allowed rupture time
Rupture.min_rup_time = 0.000001  ;
%add a constant speed vector at start/end
Rupture.final_speed = 0.045+0.0052 - 0.016 -0.0089 ;%+0.0059%+50.005 % 0.027+(10.9*0.001) ; %add number of sources / points/ timesteps in meter that work at final speed
Rupture.start_speed = 0.012   %0.01
Rupture.length_total = Rupture.final_speed+Rupture.start_speed+0.005+0.01 

Source.ReleaseTogether = 0%ignore; round(8/Trans.spacingMm) %round(12/Trans.spacingMm) ;
Source.ReleaseTogetherFactor = 0
Source.increase_length = 1  ;%5 ; %increases number of  timesteps since it
    %decreases delta_t to allow for stretching
    %becomes important if the source length becomes too small due too many
    %shifted sources. d_shift min defines minimal delta_t before scalitg by this factor
Params.observation_dir = 'vertical' ;   %'vertical' 'horizontal'  Gx - horizontal  Gy - side Gz - vertical
                                        % vertical is z and horizontal is x  , side is y    % defines the
                                         % position of the probe and thus which vp is retrieved, vx,vy,vz

                                                         
                                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SOURCE DIRECTION OBSERVATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                      
Source.Direction = 90  
Rupture.rup_direction = 'pos' ; %'pos' 'neg'  pos is right or up, neg is down or left

%define relative to observation direction. 90 means right for vertical ,
%  0 means right for horizontal; define always positive
                                                   
Params.observation_plane = 'x-z' % 'x-y'  ,'x-z' %defines the plane of observation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Params.p_speed = 1540;
Params.density = 1000;
Params.t_end = 0.016    ;% 0.025 ; %final simulation time  in seconds  - careful that not too short!!
Params.calc_type = 'discrete'  ;% % smooth should be obsolete!! 'smooth' 'discrete'
Params.video = 'off'  ;  %'on'  %%'off'
                                %Source point also automatically defines the position
                                %2 source points: above or on the side top to bottom x-z plane
                                %3 source points: x-y plane


%%%%%%%%%%%%%% SOURCE POINT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%Source.Point = [floor(0.5*Trans.Depthpoints) floor(0.3*Trans.numelements_true) 0.003];

Source.Point = [Trans.Depthpoints 1];
                                    %---> source length is calculated from this.
                                    %only use positive non-values of the source as source vector
                                    %time sampling equal to rupture time of 1 point to the next:
                                    %Define a more or less centred SourcePoint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Source.Sign = 1   ; %Direction of the Source.

Source.Magnitude = 1;  %M0 for the double couple
Source.type_switch = 'gausshape'   ;
Source.figure_switch = 'save'  ;  %'display', 'no' ,'save'
Source.displaytype = 'speed' ; 
% Source.Trise = 0.5 ;  %defines tris in ramp. make it 0.5 for gaussian
Source.Trise = 0.1 % new way of source
Source.GaussShift = 0.5 ; 



% the smaller, the longer Source.signal_length
Source.d_shift_min = 0.022*Source.Trise ; % start  or end depends if accel or decel
Source.d_shift_min = 0.04*Source.Trise ; % start  or end depends if accel or decel
Source.d_shift_min = 0.012*0.5 %Source.Trise ; % start  or end depends if accel or decel

% Source.d_shift_start = 2*Source.Trise ; %
Params.filename = strcat('Force_', num2str(Rupture.rup_speed_start_factor), '_', ... 
    num2str(Rupture.rup_speed_end_factor), '_', Params.observation_dir, '_',  ... 
    Params.observation_plane, '_', Params.Vp_Dpl , '_', num2str(Params.calc_type),...
        '_', num2str(Source.Sign), '_', Rupture.rup_direction, '_', num2str(Source.d_shift_min)  )  ;
Params.filename = strrep(Params.filename,'.','-');

Params.testrun = 'off' ;%'off'
Params.saveSwitch = 'off'
Params.Field_only = 'off'
%%
if isequal(Params.testrun,'on') == 1
    [~,~,~,Source, Params, Rupture] =  wrapper_function(Trans, Source, Params,  Rupture) ;
else
    if isequal(Params.Field_only,'on') 
        [ ~,Fieldterms,~,Source, Params, Rupture] = wrapper_function(...
        Trans, Source, Params,  Rupture)   ;
    else
        [Ut, Fieldterms, vt, Source, Params, Rupture] = wrapper_function(...
        Trans, Source, Params,  Rupture)   ;
    end
end


return

%%
shiftnr = -40
XposVector =Params.padnumber.space+Source.Point(2)-shiftnr:Params.padnumber.space...
    +Source.Point(2)-shiftnr+127;

x2= Params.x_vec(XposVector-Params.padnumber.space)-Params.x_vec(XposVector(1)-Params.padnumber.space)  
y2= flipud(Params.depth_meter)

%show snapshots
figure; nexttile; imagesc(x2,y2, normalised_diff(-1*Ut(:,XposVector,10),99),[-1 1])
axis image

function norm_value = normalised_diff( data,PrcTile )
% Normalise values of an array to be between -1 and 1
% original sign of the array values is maintained.
%if Prctile is given not equal to 100, all values above percentile will be
%cut


if PrcTile <100
    Cut = prctile(data(:),PrcTile); 
    if Cut < 0
          max_range_value = abs(Cut(:));
          min_range_value = Cut;
    else
          max_range_value = abs(Cut(:));
          min_range_value = -Cut;
    end
    data(data<min_range_value)  = min_range_value;
    data(data>max_range_value)  = max_range_value;

    norm_value = 2 .* data ./ (max_range_value - min_range_value);
else
    

if abs(min(data(:))) > max(data(:))
      max_range_value = abs(min(data(:)));
      min_range_value = min(data(:));
  else
      max_range_value = max(data(:));
      min_range_value = -max(data(:));
  end
norm_value = 2 .* data ./ (max_range_value - min_range_value);
end

end

