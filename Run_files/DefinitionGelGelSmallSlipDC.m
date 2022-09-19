%% DEfinitions
clear all
addpath(genpath('/analytic_rupture_simulation/DC_simulation/'))
set(0,'DefaultFigureWindowStyle','Docked')
% check that the convoluted source signal in the near field has a plateau
% ---> source signal long enough
Params.Vp_Dpl = 'Vp'   % decide if you want displacement or Vp as primary output
Trans.numelements_true = 256; %90 to have a third 'outside imaging region'
Trans.size_factor = 2; %reduce size and increase distance in z
Trans.Depthpoints = 900/Trans.size_factor;  %delta depth is 1/8 of ...
        % wavelength ultrasound scaled by size factor
Trans.frequency = 5.208; % nominal frequency in MHz
Trans.Bandwidth = [4, 7] ;
Trans.elementWidth = .250; % width in mms
Trans.spacingMm = .298/2;   % Spacing between elements in mm.
Params.Fs = 3000;
% Params.returnSwitch = 'Field' ; % 'Green', 'Field' ,'none'
Params.Field_only = 'off'
Params.s_speed = 6.0;   %1 m/s error well possible at 1/3000 PRF and the measured distance of 14.40 mm


Rupture.rup_speed_start_factor = 0.6  %1.0 ;
Rupture.rup_speed_end_factor = 0.6    % 0.1  ;%sqrt(2)  % sqrt(2) ;  %sqrt(2)
%Give a maximum allowed rupture time
Rupture.max_rup_time = 0.03  ;  %total time from one end of rupture to the other
%Give a minimum allowed rupture time
Rupture.min_rup_time = 0.000001  ;
% Rupture.max_rup_steps = 250; %maximum time_points wished for calculation

%add a constant speed vector at start/end
Rupture.final_speed = 0.01 ; %add number of sources / points/ timesteps in meter that work at final speed
Rupture.start_speed = 0.0000001
Source.ReleaseTogether = 0% round(8/Trans.spacingMm) %round(12/Trans.spacingMm) ;
Source.ReleaseTogetherFactor = 0
Source.increase_length = 5  ;%5 ; %increases number of  timesteps since it
Rupture.length_total = Rupture.final_speed+Rupture.start_speed  ;% in m accelerating and steady length
    %decreases delta_t to allow for stretching
    %becomes important if the source length becomes too small due too many
    %shifted sources. d_shift min defines minimal delta_t before scalintg by this factor
Params.observation_dir = 'vertical' ;   %'vertical' 'horizontal'  Gx - horizontal  Gy - side Gz - vertical
                                                % vertical is z and horizontal is x  , side is y    % defines the
                                                         % position of the probe and thus which vp is retrievec, vx,vy,vz

                                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SOURCE DIRECTION OBSERVATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                      
Source.Direction = 90  
Rupture.rup_direction = 'pos' ; %'pos' 'neg'  pos is right or up, neg is down or left

%define relative to observation direction 90 means right for vertical ,
%  0 means right for horizontal; define always positive
                                                   
Params.observation_plane = 'x-z' % 'x-y'  ,'x-z' %defines the plane of observation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Params.s_speed = 6;
Params.s_speed = 3.725;
Params.p_speed = 1540;
Params.density = 1000;
Params.t_end = 0.014    ;% 0.025 ; %final simulation time  in seconds  - careful that not too short!!
Params.calc_type = 'discrete'  ;% % smooth should be obsolete!! 'smooth' 'discrete'
Params.video = 'off'  ;  %'on'  %%'off'
                                %Source point also automatically defines the position
                                %2 source points: above or on the side top to bottom x-z plane
                                %3 source points: x-y plane

                                %%Source.Point = [floor(0.5*Trans.Depthpoints) floor(0.3*Trans.numelements_true) 0.003];
                      %%Source.Point = [floor(0.5*Trans.Depthpoints) floor(0.3*Trans.numelements_true) 0.003];
%%%%%%%%%%%%%% SOURCE POINT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Source.Point = floor([0.5*Trans.Depthpoints 0.6*Trans.numelements_true+22]);
                                    %---> source length is calculated from this.
                                    %only use positive non-values of the source as source vector
                                    %time sampling equal to rupture time of 1 point to the next:
                                    %Define a more or less centred SourcePoint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Source.Sign = 1   ; %Direction of the Source.
% Source.Cycles = 1;
% Source.Frequency_Hz = 250;
Source.Magnitude = 1  %0.001;  %M0 for the double couple

% Source.Bandwidth = 0.6; % bandwidth of gausspuls
%for the case of a non-accelerating source ; for accelerating source: adapt
Source.type_switch = 'GaussDC'   ;
Source.figure_switch = 'save'  ;  %'display', 'no' ,'save'
Source.displaytype = 'speed' ; ;
Source.figure_switch = 'save'  ;  %'display', 'no' ,'save'
Source.Trise_old = 0.5 ;  %defines tris in ramp. make it 0.5 for gaussian
Source.Trise = 0.1 ;  %defines tris in ramp. make it 0.5 for gaussian

Source.d_shift_min =0.02*Source.Trise_old ; % start  or end depends if accel or decel
Source.d_shift_min = 0.5*Source.Trise_old ;


% Source.d_shift_start = 2*Source.Trise ; %
Params.filename = strcat('DC_', num2str(Rupture.rup_speed_start_factor), '_', ... 
    num2str(Rupture.rup_speed_end_factor), '_', Params.observation_dir, '_',  ... 
    Params.observation_plane, '_', Params.Vp_Dpl , '_', num2str(Params.calc_type),...
        '_', num2str(Source.Sign), '_', Rupture.rup_direction, '_', num2str(Source.d_shift_min)  )  ;
Params.filename = strrep(Params.filename,'.','-');
% obsolete ? Params.Field_only = 'off' %stop before shifting the sources
Params.testrun = 'off'  ;  %'on', 'off'
%%
Params.Field_only = 'off'
Params.returnSwitch = 0
if isequal(Params.testrun,'on') == 1
    [~,~,~,Source, Params, Rupture] =  wrapper_function(Trans, Source, Params,  Rupture) ;
else
    if isequal(Params.Field_only,'on') 
        [ ~,~,Fieldterms,Source, Params, Rupture] = wrapper_function(...
        Trans, Source, Params,  Rupture)   ;
    else
        [Ut, vt, Fieldterms ,  Source, Params, Rupture] = wrapper_function(...
        Trans, Source, Params,  Rupture)   ;
    end
end
%%

% new_y_vec = y_vec - y_vec(1) ;
% cutUpDown = 0.5*(y(end)-new_y_vec(end)) ; 
% CutUp = find(y >= cutUpDown,1,'first') ;
% CutDown = find(y >= (y(end)-cutUpDown),1,'first') ;
% y_sim_new = y(CutUp:CutDown)-y(CutUp)
% 
% % cutDepth = find(new_y_vec >= y(end),1,'first')
% % y2 = y2(1:cutDepth) ;
% % adapt to length of simulation
% Sim_extract = resize(Sim_extract(CutUp:CutDown,:,:) , [length(new_y_vec) size(p4,2) size(p4,3)]);
% y2 = resize(y_sim_new,length(new_y_vec))
% x2 = resize(x,length(x_vec))


%%
x2= [0:Trans.spacingMm:Trans.spacingMm.*length(...
    Params.padnumber.space+Params.Trans.numelements_true:size(vt,2))]./1000
y2= flipud(Params.depth_meter)

%show snapshots
figure; nexttile; imagesc(x2,y2, normalised_diff(-1*vt(...
    :,Params.padnumber.space+Params.Trans.numelements_true:end,10),99),[-1 1])
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

